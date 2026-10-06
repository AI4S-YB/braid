//! Directed graph used by the splice graph and the path graph.
//!
//! Source and sink occupy ids 0 and 1 and are not counted by `len`.
//! Successor lists stay sorted by id. Depth-first topological order then
//! walks the highest id first, which is reproducible across runs.

use std::collections::HashMap;
use std::hash::Hash;

use crate::util::insert_sorted;

#[derive(Clone, Debug)]
pub(crate) struct Graph<T> {
    map: HashMap<T, usize>,
    nodes: Vec<Option<T>>,
    pub succs: Vec<Vec<usize>>,
    pub preds: Vec<Vec<usize>>,
    n: usize,
}

impl<T: Clone + Eq + Hash> Graph<T> {
    pub(crate) fn new(source: T, sink: T) -> Self {
        Self {
            map: HashMap::from([(source.clone(), 0), (sink.clone(), 1)]),
            nodes: vec![Some(source), Some(sink)],
            succs: vec![Vec::new(), Vec::new()],
            preds: vec![Vec::new(), Vec::new()],
            n: 0,
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.n
    }

    pub(crate) fn node_ids(&self) -> impl Iterator<Item = usize> + '_ {
        (2..self.nodes.len()).filter(|&id| self.nodes[id].is_some())
    }

    pub(crate) fn get(&self, id: usize) -> Option<&T> {
        self.nodes.get(id)?.as_ref()
    }

    pub(crate) fn id_of(&self, node: &T) -> Option<usize> {
        self.map.get(node).copied()
    }

    pub(crate) fn add_node(&mut self, node: T) -> usize {
        if let Some(&id) = self.map.get(&node) {
            return id;
        }
        let id = self.nodes.len();
        self.map.insert(node.clone(), id);
        self.nodes.push(Some(node));
        self.succs.push(Vec::new());
        self.preds.push(Vec::new());
        self.n += 1;
        id
    }

    pub(crate) fn add_path(&mut self, nodes: &[T]) -> Vec<usize> {
        if nodes.is_empty() {
            return Vec::new();
        }
        let mut ids = Vec::with_capacity(nodes.len());
        let mut previous = self.add_node(nodes[0].clone());
        ids.push(previous);
        for node in &nodes[1..] {
            let id = self.add_node(node.clone());
            ids.push(id);
            insert_sorted(&mut self.succs[previous], id);
            insert_sorted(&mut self.preds[id], previous);
            previous = id;
        }
        ids
    }

    pub(crate) fn remove_node_id(&mut self, id: usize) {
        if id <= 1 || self.nodes.get(id).and_then(Option::as_ref).is_none() {
            return;
        }
        let Some(node) = self.nodes[id].take() else {
            return;
        };
        self.map.remove(&node);
        let succs = self.succs[id].clone();
        let preds = self.preds[id].clone();
        for succ in succs {
            self.preds[succ].retain(|&pred| pred != id);
        }
        for pred in preds {
            self.succs[pred].retain(|&succ| succ != id);
        }
        self.n -= 1;
    }

    pub(crate) fn unreachable(&self) -> Vec<usize> {
        let from_source = reachable(0, &self.succs);
        let to_sink = reachable(1, &self.preds);
        let mut lost: Vec<usize> = self
            .map
            .values()
            .copied()
            .filter(|&id| {
                let from = from_source.get(id).copied().unwrap_or(false);
                let to = to_sink.get(id).copied().unwrap_or(false);
                !from || !to
            })
            .collect();
        lost.sort_unstable();
        lost
    }

    pub(crate) fn is_valid(&self) -> bool {
        if self.nodes[0].is_none() || self.nodes[1].is_none() {
            return false;
        }
        let mut pred = vec![false; self.nodes.len()];
        let mut succ = vec![false; self.nodes.len()];
        pred[0] = true;
        succ[1] = true;
        let mut forward = vec![0usize];
        let mut reverse = vec![1usize];
        while !forward.is_empty() && !reverse.is_empty() {
            if forward.len() <= reverse.len() {
                let this = std::mem::take(&mut forward);
                for v in this {
                    for &w in &self.succs[v] {
                        if !pred[w] {
                            forward.push(w);
                            pred[w] = true;
                        }
                        if succ[w] {
                            return true;
                        }
                    }
                }
            } else {
                let this = std::mem::take(&mut reverse);
                for v in this {
                    for &w in &self.preds[v] {
                        if !succ[w] {
                            succ[w] = true;
                            reverse.push(w);
                        }
                        if pred[w] {
                            return true;
                        }
                    }
                }
            }
        }
        false
    }

    /// Depth-first topological order, highest successor id first.
    pub(crate) fn topo(&self) -> Result<Vec<usize>, String> {
        let n = self.nodes.len();
        let mut explored = vec![false; n];
        let mut order = Vec::new();
        let limit = n.saturating_mul(n).saturating_add(n).max(1);
        let mut steps = 0usize;
        for i in 0..n {
            if self.nodes[i].is_none() || explored[i] {
                continue;
            }
            let mut fringe = vec![i];
            while let Some(&j) = fringe.last() {
                steps += 1;
                if steps > limit {
                    return Err("transcript graph has a cycle".to_string());
                }
                if explored[j] {
                    fringe.pop();
                    continue;
                }
                let before = fringe.len();
                for &succ in &self.succs[j] {
                    if self.nodes[succ].is_some() && !explored[succ] {
                        fringe.push(succ);
                    }
                }
                if fringe.len() == before {
                    explored[j] = true;
                    order.push(j);
                    fringe.pop();
                }
            }
        }
        order.reverse();
        Ok(order)
    }

    pub(crate) fn components(&self) -> Vec<i32> {
        let mut component = vec![-1i32; self.nodes.len()];
        let mut seen = vec![false; self.nodes.len()];
        let mut next_id = 0i32;
        for start in self.node_ids() {
            if seen[start] {
                continue;
            }
            let mut stack = vec![start];
            while let Some(node) = stack.pop() {
                if seen[node] {
                    continue;
                }
                seen[node] = true;
                component[node] = next_id;
                for &next in self.succs[node].iter().rev() {
                    if !seen[next] {
                        stack.push(next);
                    }
                }
                for &next in self.preds[node].iter().rev() {
                    if !seen[next] {
                        stack.push(next);
                    }
                }
            }
            next_id += 1;
        }
        component
    }
}

fn reachable(start: usize, neighbors: &[Vec<usize>]) -> Vec<bool> {
    let mut seen = vec![false; neighbors.len()];
    let mut queue = vec![start];
    seen[start] = true;
    while let Some(node) = queue.pop() {
        for &next in &neighbors[node] {
            if !seen[next] {
                seen[next] = true;
                queue.push(next);
            }
        }
    }
    seen
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;

    #[test]
    fn topological_order_respects_edges() {
        let mut graph = Graph::new(-1, -2);
        graph.add_path(&[-1, 2, 4, -2]);
        graph.add_path(&[-1, 3, 4, -2]);
        let order = graph.topo().unwrap();
        let pos: HashMap<_, _> = order.iter().enumerate().map(|(i, &n)| (n, i)).collect();
        for (from, tos) in graph.succs.iter().enumerate() {
            for &to in tos {
                assert!(pos[&from] < pos[&to], "{from} -> {to} in {order:?}");
            }
        }
        assert!(graph.is_valid());
    }
}
