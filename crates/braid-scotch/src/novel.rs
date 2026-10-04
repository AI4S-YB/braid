//! Novel isoform discovery.
//!
//! Reads that match no annotated isoform become nodes in a similarity graph.
//! An edge weight is the number of sub-exons where both reads agree, and only
//! when read `i` does not claim an exon that read `j` skips. Communities come
//! from one pass of deterministic Louvain: nodes move in index order, and a
//! tie keeps the lowest community id. A community becomes an isoform when it
//! holds more than one read. The consensus exon is included when `1`s strictly
//! outnumber `-1`s. A discovery chunk is kept only when it assigns at least
//! `min_assigned` reads, which is 10 in SCOTCH.

use crate::gtf::{novel_name, subexon_indexes};

#[derive(Clone, Debug)]
pub struct NovelIsoform {
    pub name: String,
    pub exons: Vec<usize>,
    pub vector: Vec<i8>,
}

#[derive(Clone, Debug)]
pub enum NovelCall {
    Isoform(String),
    Uncategorized,
}

/// Assign every novel read to a discovered isoform, or leave it uncategorized.
pub fn discover(
    reads: &[(String, Vec<i8>)],
    min_assigned: usize,
    max_rounds: usize,
) -> (Vec<NovelIsoform>, Vec<NovelCall>) {
    if reads.len() <= 1 {
        return (
            Vec::new(),
            reads.iter().map(|_| NovelCall::Uncategorized).collect(),
        );
    }
    let mut pending: Vec<usize> = (0..reads.len()).collect();
    let mut isoforms: Vec<NovelIsoform> = Vec::new();
    let mut rounds = 0usize;
    while !pending.is_empty() && rounds < max_rounds {
        let chunk = if pending.len() > 1500 {
            pending[..1500].to_vec()
        } else {
            pending.clone()
        };
        let found = cluster(
            &chunk
                .iter()
                .map(|&index| reads[index].1.clone())
                .collect::<Vec<_>>(),
        );
        let assignments = assign_to_isoforms(&found, &pending, reads);
        let assigned: Vec<usize> = assignments
            .iter()
            .enumerate()
            .filter_map(|(slot, call)| call.is_some().then_some(pending[slot]))
            .collect();
        if assigned.len() >= min_assigned {
            for isoform in found {
                if isoforms.iter().all(|have| have.name != isoform.name) {
                    isoforms.push(isoform);
                }
            }
            pending.retain(|index| !assigned.contains(index));
            rounds += 1;
            if pending.is_empty() {
                break;
            }
        } else if chunk.len() == pending.len() {
            break;
        } else {
            pending.drain(..chunk.len());
        }
    }
    if isoforms.is_empty() {
        return (
            Vec::new(),
            reads.iter().map(|_| NovelCall::Uncategorized).collect(),
        );
    }
    isoforms.sort_by(|a, b| b.exons.len().cmp(&a.exons.len()).then(a.name.cmp(&b.name)));
    let calls = assign_to_isoforms(&isoforms, &(0..reads.len()).collect::<Vec<_>>(), reads)
        .into_iter()
        .map(|call| match call {
            Some(name) => NovelCall::Isoform(name),
            None => NovelCall::Uncategorized,
        })
        .collect();
    (isoforms, calls)
}

fn assign_to_isoforms(
    isoforms: &[NovelIsoform],
    indexes: &[usize],
    reads: &[(String, Vec<i8>)],
) -> Vec<Option<String>> {
    if isoforms.is_empty() {
        return vec![None; indexes.len()];
    }
    let mut hits = vec![vec![false; isoforms.len()]; indexes.len()];
    for (row, &index) in indexes.iter().enumerate() {
        for (col, isoform) in isoforms.iter().enumerate() {
            hits[row][col] = !vector_conflict(&isoform.vector, &reads[index].1);
        }
    }
    let mut weight = vec![0usize; isoforms.len()];
    for row in &hits {
        for (col, on) in row.iter().enumerate() {
            if *on {
                weight[col] += 1;
            }
        }
    }
    hits.into_iter()
        .map(|row| {
            let mut best: Option<(usize, usize)> = None;
            for (col, on) in row.into_iter().enumerate() {
                if on && best.is_none_or(|(count, _)| weight[col] > count) {
                    best = Some((weight[col], col));
                }
            }
            best.map(|(_, col)| isoforms[col].name.clone())
        })
        .collect()
}

fn vector_conflict(isoform: &[i8], read: &[i8]) -> bool {
    isoform
        .iter()
        .zip(read)
        .any(|(&left, &right)| (left == 1 && right == -1) || (left == -1 && right == 1))
}

fn cluster(vectors: &[Vec<i8>]) -> Vec<NovelIsoform> {
    let n = vectors.len();
    if n == 0 {
        return Vec::new();
    }
    let mut adj = vec![Vec::<(usize, f64)>::new(); n];
    for i in 0..n {
        for j in (i + 1)..n {
            let weight = similarity(&vectors[i], &vectors[j]);
            if weight > 0.0 {
                adj[i].push((j, weight));
                adj[j].push((i, weight));
            }
        }
    }
    let communities = louvain(&adj);
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for (node, community) in communities.iter().enumerate() {
        if groups.len() <= *community {
            groups.resize(community + 1, Vec::new());
        }
        groups[*community].push(node);
    }
    let mut isoforms = Vec::new();
    for group in groups {
        if group.len() <= 1 {
            continue;
        }
        let width = vectors[group[0]].len();
        let mut consensus = vec![-1i8; width];
        let mut bits = vec![false; width];
        for exon in 0..width {
            let ones = group
                .iter()
                .filter(|&&node| vectors[node][exon] == 1)
                .count();
            let minus = group
                .iter()
                .filter(|&&node| vectors[node][exon] == -1)
                .count();
            if ones > minus {
                consensus[exon] = 1;
                bits[exon] = true;
            }
        }
        if !bits.iter().any(|&on| on) {
            continue;
        }
        let name = novel_name(&bits);
        if isoforms.iter().any(|iso: &NovelIsoform| iso.name == name) {
            continue;
        }
        isoforms.push(NovelIsoform {
            name,
            exons: subexon_indexes(&bits),
            vector: consensus,
        });
    }
    isoforms
}

/// Similarity used by SCOTCH: agreement count, or zero when `left` has a `1`
/// where `right` has a `-1`. The check is not symmetric.
fn similarity(left: &[i8], right: &[i8]) -> f64 {
    if left.iter().zip(right).any(|(&a, &b)| a == 1 && b == -1) {
        return 0.0;
    }
    left.iter()
        .zip(right)
        .filter(|(&a, &b)| (a == 1 && b == 1) || (a == -1 && b == -1))
        .count() as f64
}

fn louvain(adj: &[Vec<(usize, f64)>]) -> Vec<usize> {
    let n = adj.len();
    let mut community: Vec<usize> = (0..n).collect();
    if n == 0 {
        return community;
    }
    let degree: Vec<f64> = adj
        .iter()
        .map(|edges| edges.iter().map(|edge| edge.1).sum())
        .collect();
    let m2: f64 = degree.iter().sum();
    if m2 == 0.0 {
        return community;
    }
    let m = m2 / 2.0;
    let mut sigma = degree.clone();
    let mut moved = true;
    let mut guard = 0;
    while moved && guard < 50 {
        moved = false;
        guard += 1;
        for node in 0..n {
            let current = community[node];
            sigma[current] -= degree[node];
            let mut kin: Vec<(usize, f64)> = Vec::new();
            for &(neighbor, weight) in &adj[node] {
                let comm = community[neighbor];
                if let Some(slot) = kin.iter_mut().find(|(id, _)| *id == comm) {
                    slot.1 += weight;
                } else {
                    kin.push((comm, weight));
                }
            }
            let stay = kin
                .iter()
                .find(|(id, _)| *id == current)
                .map(|(_, weight)| *weight)
                .unwrap_or(0.0);
            let stay_gain = gain(stay, sigma[current], degree[node], m);
            let mut best = current;
            let mut best_gain = stay_gain;
            for (comm, weight) in kin {
                if comm == current {
                    continue;
                }
                let next = gain(weight, sigma[comm], degree[node], m);
                if next > best_gain + 1e-12 || ((next - best_gain).abs() <= 1e-12 && comm < best) {
                    best = comm;
                    best_gain = next;
                }
            }
            if best != current && best_gain > 0.0 {
                community[node] = best;
                sigma[best] += degree[node];
                moved = true;
            } else {
                community[node] = current;
                sigma[current] += degree[node];
            }
        }
    }
    densify(&community)
}

fn gain(kin: f64, sigma_tot: f64, degree: f64, m: f64) -> f64 {
    kin / m - sigma_tot * degree / (2.0 * m * m)
}

fn densify(community: &[usize]) -> Vec<usize> {
    let mut map = Vec::new();
    community
        .iter()
        .map(|&id| {
            if let Some(pos) = map.iter().position(|&have| have == id) {
                pos
            } else {
                map.push(id);
                map.len() - 1
            }
        })
        .collect()
}

/// Group a shorter novel isoform into a longer one when the extra exons sit
/// entirely on the strand's outer end, which SCOTCH treats as truncation.
pub fn group_novel(isoforms: &[(String, Vec<usize>)], strand: char) -> Vec<(String, String)> {
    let mut pool: Vec<(String, Vec<usize>)> = isoforms.to_vec();
    pool.sort_by(|a, b| {
        a.1.len()
            .cmp(&b.1.len())
            .then_with(|| cmp_novel_id(&a.0, &b.0))
    });
    let mut groups: Vec<Vec<(String, Vec<usize>)>> = Vec::new();
    let mut extending = false;
    while !pool.is_empty() {
        if !extending {
            let child = pool.remove(0);
            if let Some(pos) = pool
                .iter()
                .position(|parent| truncation_parent(&child.1, &parent.1, strand))
            {
                let parent = pool.remove(pos);
                groups.push(vec![child, parent]);
                extending = true;
            } else {
                groups.push(vec![child]);
                extending = false;
            }
        } else {
            let child = groups.last().unwrap().last().unwrap().1.clone();
            if let Some(pos) = pool
                .iter()
                .position(|parent| truncation_parent(&child, &parent.1, strand))
            {
                groups.last_mut().unwrap().push(pool.remove(pos));
            } else {
                extending = false;
            }
        }
    }
    let mut map = Vec::new();
    for group in groups {
        let kept = group.last().unwrap().0.clone();
        for (name, _) in group {
            map.push((name, kept.clone()));
        }
    }
    map
}

fn novel_suffix(name: &str) -> &str {
    name.rsplit_once('_').map(|(_, id)| id).unwrap_or(name)
}

fn cmp_novel_id(left: &str, right: &str) -> std::cmp::Ordering {
    let (left_id, right_id) = (novel_suffix(left), novel_suffix(right));
    left_id
        .len()
        .cmp(&right_id.len())
        .then(left_id.cmp(right_id))
        .then(left.cmp(right))
}

fn truncation_parent(child: &[usize], parent: &[usize], strand: char) -> bool {
    if child.is_empty() || !child.iter().all(|exon| parent.contains(exon)) {
        return false;
    }
    let diff: Vec<usize> = parent
        .iter()
        .copied()
        .filter(|exon| !child.contains(exon))
        .collect();
    let min = *child.iter().min().unwrap();
    let max = *child.iter().max().unwrap();
    match strand {
        '+' => diff.iter().all(|exon| *exon < min),
        '-' => diff.iter().all(|exon| *exon > max),
        _ => true,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_vectors_form_one_isoform_and_truncations_group() {
        let vector = vec![1, -1, 1];
        let reads: Vec<_> = (0..10).map(|i| (format!("r{i}"), vector.clone())).collect();
        let (isoforms, calls) = discover(&reads, 10, 50);
        assert_eq!(isoforms.len(), 1);
        assert_eq!(isoforms[0].name, "novelIsoform_5");
        assert!(calls
            .iter()
            .all(|call| matches!(call, NovelCall::Isoform(_))));
        let map = group_novel(
            &[
                ("novelIsoform_2".to_string(), vec![1]),
                ("novelIsoform_3".to_string(), vec![0, 1]),
                ("novelIsoform_1".to_string(), vec![0]),
            ],
            '+',
        );
        let kept = |name: &str| map.iter().find(|(from, _)| from == name).unwrap().1.clone();
        assert_eq!(kept("novelIsoform_2"), "novelIsoform_3");
        assert_eq!(kept("novelIsoform_1"), "novelIsoform_1");
    }

    #[test]
    fn nine_reads_stay_uncategorized() {
        let reads: Vec<_> = (0..9).map(|i| (format!("r{i}"), vec![1, -1])).collect();
        let (isoforms, calls) = discover(&reads, 10, 50);
        assert!(isoforms.is_empty());
        assert_eq!(calls.len(), 9);
    }
}
