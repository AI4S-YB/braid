//! Native TAMA merge, following tama_merge.py at 2fa3c308.
//! BED coordinates stay zero-based here; shared Transcript coordinates are
//! converted to 1-based only at the public API boundary.
use crate::group::{gene_group, sort_position_keys, GroupedTranscript, TranscriptGeometry};
use crate::ids::{keys, reinsert, Ids};
use crate::Py27Dict;
use plena_io::suffix_path;
use plena_model::{Exon, MergeSource, Strand, Transcript};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs;
use std::path::Path;

fn ints(v: &[i64]) -> String {
    v.iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

#[derive(Clone, Debug)]
pub struct MergeSettings {
    pub ends: String,
    pub five_prime: i64,
    pub exon_diff: i64,
    pub three_prime: i64,
    pub duplicates: String,
    pub source_id: String,
    pub cds: String,
}
impl Default for MergeSettings {
    fn default() -> Self {
        Self {
            ends: "common_ends".into(),
            five_prime: 20,
            exon_diff: 10,
            three_prime: 20,
            duplicates: "no_merge".into(),
            source_id: "no_source_id".into(),
            cds: "no_cds".into(),
        }
    }
}
#[derive(Clone, Debug, Default)]
pub struct MergeTexts {
    pub bed: String,
    pub trans_report: String,
    pub gene_report: String,
    pub merge: String,
    pub transcripts: Vec<Transcript>,
}

#[derive(Clone)]
struct Model {
    id: String,
    gene: String,
    trans: String,
    source: String,
    chrom: String,
    strand: String,
    start: i64,
    end: i64,
    starts: Vec<i64>,
    ends: Vec<i64>,
    cds_start: i64,
    cds_end: i64,
    capped: bool,
    priority: [i64; 3],
}
impl Model {
    fn n(&self) -> usize {
        self.starts.len()
    }
}

impl GroupedTranscript for Model {
    fn geometry(&self) -> TranscriptGeometry<'_> {
        TranscriptGeometry {
            id: &self.id,
            starts: &self.starts,
            ends: &self.ends,
        }
    }
}

// Same, degraded with equal exon count, degraded with fewer exons, or no match.
#[derive(PartialEq)]
enum Match {
    Same,
    Short,
    Fewer,
    Different,
}
fn compare(a: &Model, b: &Model, p: &MergeSettings) -> Match {
    if a.capped && b.capped && a.n() != b.n() {
        return Match::Different;
    }
    let (long, short) = if a.capped {
        (a, b)
    } else if b.capped {
        (b, a)
    } else if a.n() > b.n()
        || (a.n() == b.n()
            && if a.strand == "+" {
                a.starts[0] <= b.starts[0]
            } else {
                a.ends.last() >= b.ends.last()
            })
    {
        (a, b)
    } else {
        (b, a)
    };
    if long.n() < short.n() {
        return Match::Different;
    }
    for i in 0..short.n() {
        let li = if a.strand == "+" { long.n() - 1 - i } else { i };
        let si = if a.strand == "+" {
            short.n() - 1 - i
        } else {
            i
        };
        let (ls, le, ss, se) = (
            long.starts[li],
            long.ends[li],
            short.starts[si],
            short.ends[si],
        );
        let st = if a.strand == "-" && i == 0 {
            p.three_prime
        } else if a.strand == "+" && i + 1 == long.n() {
            p.five_prime
        } else {
            p.exon_diff
        };
        let et = if a.strand == "+" && i == 0 {
            p.three_prime
        } else if a.strand == "-" && i + 1 == long.n() {
            p.five_prime
        } else {
            p.exon_diff
        };
        let sm = (ls - ss).abs() <= st;
        let em = (le - se).abs() <= et;
        if (a.capped && b.capped) || i + 1 < short.n() {
            if !sm || !em {
                return Match::Different;
            }
        } else {
            let five_match = if a.strand == "+" { sm } else { em };
            if a.strand == "+" {
                if !em || (!sm && ls > ss) {
                    return Match::Different;
                }
            } else if !sm || (!em && le < se) {
                return Match::Different;
            }
            return if long.n() != short.n() {
                Match::Fewer
            } else if five_match {
                Match::Same
            } else {
                Match::Short
            };
        }
    }
    Match::Same
}

#[derive(Default)]
struct Groups {
    count: i64,
    trans: Py27Dict<String, Py27Dict<i64, i32>>,
    members: Py27Dict<i64, Ids>,
    max: HashMap<i64, usize>,
}
impl Groups {
    fn new_group(&mut self, m: &Model) {
        self.count += 1;
        let n = self.count;
        self.trans
            .or_insert_with(m.id.clone(), Py27Dict::new)
            .insert(n, 1);
        self.members
            .or_insert_with(n, Ids::new)
            .insert(m.id.clone(), 1);
        self.max.insert(n, m.n());
    }
    fn groups(&self, id: &str) -> Vec<i64> {
        self.trans
            .get_str(id)
            .unwrap()
            .iter()
            .map(|(g, _)| *g)
            .collect()
    }
    fn same(&self, a: &str, b: &str) -> bool {
        self.groups(a)
            .iter()
            .any(|g| self.trans.get_str(b).unwrap().contains(g))
    }
    fn add(&mut self, a: &str, b: &str) {
        let gs = self.groups(a);
        if gs.len() == 1 && self.members.get(&gs[0]).unwrap().len() == 1 {
            self.members.remove(&gs[0]);
            self.trans.remove(&a.to_string());
        }
        for g in self.groups(b) {
            self.trans
                .or_insert_with(a.to_string(), Py27Dict::new)
                .insert(g, 1);
            self.members.get_mut(&g).unwrap().insert(a.to_string(), 1);
        }
    }
    fn merge(&mut self, a: &Model, b: &Model, nocap: bool) {
        self.count += 1;
        let dest = self.count;
        self.members.insert(dest, Ids::new());
        self.max.insert(dest, 0);
        for m in [a, b] {
            for g in self.groups(&m.id) {
                if g == dest || (nocap && self.max[&g] != m.n()) {
                    continue;
                }
                for id in keys(self.members.get(&g).unwrap()) {
                    self.members.get_mut(&dest).unwrap().insert(id.clone(), 1);
                    let gs = self.trans.get_str_mut(&id).unwrap();
                    gs.remove(&g);
                    gs.insert(dest, 1);
                }
                self.max.insert(dest, self.max[&dest].max(self.max[&g]));
                self.members.remove(&g);
            }
        }
        if self.members.get(&dest).unwrap().is_empty() {
            self.members.remove(&dest);
        }
    }
    // Upstream check_trans_group shadows its argument and checks ALL groups.
    fn any_singleton(&self) -> bool {
        self.trans.iter().any(|(_, gs)| {
            gs.len() == 1 && self.members.get(gs.iter().next().unwrap().0).unwrap().len() == 1
        })
    }
}

fn simplify(models: &[&Model], p: &MergeSettings) -> Py27Dict<i64, Ids> {
    let index: HashMap<&str, &Model> = models.iter().map(|m| (m.id.as_str(), *m)).collect();
    let mut g = Groups::default();
    for m in models {
        g.new_group(m);
    }
    let caps: Vec<&Model> = models.iter().copied().filter(|m| m.capped).collect();
    let nocaps: Vec<&Model> = models.iter().copied().filter(|m| !m.capped).collect();
    let mut ungrouped = Ids::new();
    for m in &caps {
        ungrouped.insert(m.id.clone(), 1);
    }
    let mut pending = Ids::new();
    while !ungrouped.is_empty() {
        let mut ids = keys(&ungrouped);
        ids.sort();
        let mut h = ids[0].clone();
        ungrouped.remove(&h);
        loop {
            for (id, _) in ungrouped.iter() {
                if id != &h
                    && !g.same(&h, id)
                    && compare(index[h.as_str()], index[id.as_str()], p) == Match::Same
                {
                    g.merge(index[h.as_str()], index[id.as_str()], false);
                    pending.insert(id.clone(), 1);
                }
            }
            for (id, _) in pending.iter() {
                ungrouped.remove(id);
            }
            if pending.is_empty() || ungrouped.is_empty() {
                break;
            }
            let mut ids = keys(&pending);
            ids.sort();
            h = ids[0].clone();
            pending.remove(&h);
        }
    }
    for h in &caps {
        for prey in &nocaps {
            if !g.same(&h.id, &prey.id) && compare(h, prey, p) != Match::Different {
                g.add(&prey.id, &h.id);
            }
        }
    }
    if !g.any_singleton() {
        return g.members;
    }
    let mut levels: BTreeMap<usize, Ids> = BTreeMap::new();
    for m in &nocaps {
        levels.entry(m.n()).or_default().insert(m.id.clone(), 1);
    }
    let mut used = HashSet::new();
    let mut degraded = HashSet::new();
    let mut pending = Ids::new();
    for (&level, ids) in levels.iter().rev() {
        let mut ungrouped = reinsert(ids);
        while !ungrouped.is_empty() {
            let mut by_five: BTreeMap<i64, Ids> = BTreeMap::new();
            for (id, _) in ungrouped.iter() {
                let m = index[id.as_str()];
                let pos = if m.strand == "+" { m.start } else { -m.end };
                by_five.entry(pos).or_default().insert(id.clone(), 1);
            }
            let mut h = by_five
                .first_key_value()
                .unwrap()
                .1
                .iter()
                .next()
                .unwrap()
                .0
                .clone();
            ungrouped.remove(&h);
            loop {
                used.insert(h.clone());
                let hunter = index[h.as_str()];
                if !degraded.contains(&h) {
                    for (id, _) in ungrouped.iter() {
                        if id == &h || g.same(&h, id) {
                            continue;
                        }
                        match compare(hunter, index[id.as_str()], p) {
                            Match::Same => g.merge(hunter, index[id.as_str()], true),
                            Match::Short => g.add(id, &h),
                            _ => continue,
                        }
                        if !used.contains(id) {
                            pending.insert(id.clone(), 1);
                        }
                    }
                }
                for (_, small) in levels.range(..level).rev() {
                    for (id, _) in reinsert(small).iter() {
                        if !g.same(&h, id) && compare(hunter, index[id.as_str()], p) == Match::Fewer
                        {
                            g.add(id, &h);
                            degraded.insert(id.clone());
                        }
                    }
                }
                if pending.is_empty() {
                    break;
                }
                let mut ids = keys(&pending);
                ids.sort();
                h = ids[0].clone();
                pending.remove(&h);
            }
        }
    }
    g.members
}

#[derive(Clone)]
struct Collapsed {
    members: Ids,
    order: Vec<String>,
    starts: Vec<i64>,
    ends: Vec<i64>,
    sw: Vec<i64>,
    ew: Vec<i64>,
    start_support: Vec<String>,
    end_support: Vec<String>,
}
fn collapse(
    ids: Vec<String>,
    all: &HashMap<String, Model>,
    longest: bool,
) -> Result<Collapsed, String> {
    let mut members = Ids::new();
    let mut order = Vec::new();
    for id in ids {
        if !members.contains(&id) {
            members.insert(id.clone(), 1);
            order.push(id);
        }
    }
    let models: Vec<&Model> = order.iter().map(|id| &all[id]).collect();
    let n = models
        .iter()
        .map(|m| m.n())
        .max()
        .ok_or("empty merge group")?;
    let plus = models[0].strand == "+";
    let mut starts = Vec::new();
    let mut ends = Vec::new();
    let mut sw = Vec::new();
    let mut ew = Vec::new();
    let mut start_ids: HashMap<i64, BTreeSet<String>> = HashMap::new();
    let mut end_ids: HashMap<i64, BTreeSet<String>> = HashMap::new();
    if models.len() == 1 {
        starts = models[0].starts.clone();
        ends = models[0].ends.clone();
        sw = vec![0; n];
        ew = vec![0; n];
        for &pos in &starts {
            start_ids.insert(pos, BTreeSet::from([models[0].id.clone()]));
        }
        for &pos in &ends {
            end_ids.insert(pos, BTreeSet::from([models[0].id.clone()]));
        }
    } else {
        for i in 0..n {
            let mut sc: BTreeMap<i64, BTreeMap<i64, usize>> = BTreeMap::new();
            let mut ec: BTreeMap<i64, BTreeMap<i64, usize>> = BTreeMap::new();
            for m in &models {
                if i >= m.n() {
                    continue;
                }
                let j = if plus { m.n() - 1 - i } else { i };
                let (a, b) = (m.starts[j], m.ends[j]);
                let [five, sj, three] = m.priority;
                let (sp, ep) = if i + 1 == n {
                    if plus {
                        (five, sj)
                    } else {
                        (sj, five)
                    }
                } else if i + 1 == m.n() {
                    if m.capped {
                        return Err("truncated capped model in merge group".into());
                    }
                    if plus {
                        (999, sj)
                    } else {
                        (sj, 999)
                    }
                } else if i > 0 {
                    (sj, sj)
                } else if plus {
                    (sj, three)
                } else {
                    (three, sj)
                };
                for (pos, pri, votes, support) in [
                    (a, sp, &mut sc, &mut start_ids),
                    (b, ep, &mut ec, &mut end_ids),
                ] {
                    let bucket = votes.entry(pri).or_default();
                    if !bucket.contains_key(&pos) {
                        support.insert(pos, BTreeSet::new());
                    }
                    *bucket.entry(pos).or_default() += 1;
                    support.get_mut(&pos).unwrap().insert(m.id.clone());
                }
            }
            let sc = sc.first_key_value().unwrap().1;
            let ec = ec.first_key_value().unwrap().1;
            let smin = *sc.first_key_value().unwrap().0;
            let smax = *sc.last_key_value().unwrap().0;
            let emin = *ec.first_key_value().unwrap().0;
            let emax = *ec.last_key_value().unwrap().0;
            let maxs = *sc.values().max().unwrap();
            let maxe = *ec.values().max().unwrap();
            let mut a = *sc.iter().find(|(_, v)| **v == maxs).unwrap().0;
            let mut b = *ec.iter().rev().find(|(_, v)| **v == maxe).unwrap().0;
            if longest {
                if (plus && i + 1 == n) || (!plus && i == 0) {
                    a = smin;
                }
                if (!plus && i + 1 == n) || (plus && i == 0) {
                    b = emax;
                }
            }
            starts.push(a);
            ends.push(b);
            sw.push(smax - smin);
            ew.push(emax - emin);
        }
    }
    let sorted = |coords: Vec<i64>, w: Vec<i64>| {
        let mut pairs: Vec<_> = coords.into_iter().zip(w).collect();
        pairs.sort();
        pairs.into_iter().unzip()
    };
    let (starts, sw): (Vec<_>, Vec<_>) = sorted(starts, sw);
    let (ends, ew): (Vec<_>, Vec<_>) = sorted(ends, ew);
    let support = |coords: &[i64], table: &HashMap<i64, BTreeSet<String>>| {
        coords
            .iter()
            .map(|c| table[c].iter().cloned().collect::<Vec<_>>().join(","))
            .collect()
    };
    Ok(Collapsed {
        start_support: support(&starts, &start_ids),
        end_support: support(&ends, &end_ids),
        members,
        order,
        starts,
        ends,
        sw,
        ew,
    })
}

fn position_key(c: &Collapsed) -> String {
    c.starts
        .iter()
        .zip(&c.ends)
        .flat_map(|(&a, &b)| [a, b])
        .map(|n| n.to_string())
        .collect::<Vec<_>>()
        .join(",")
}
fn collapse_groups(
    groups: Py27Dict<i64, Ids>,
    all: &HashMap<String, Model>,
    p: &MergeSettings,
) -> Result<Vec<Collapsed>, String> {
    let mut positions: Py27Dict<String, Collapsed> = Py27Dict::new();
    for (_, ids) in groups.iter() {
        // process_trans_group explicitly resets this to common_ends upstream.
        let c = collapse(keys(ids), all, false)?;
        let key = position_key(&c);
        if let Some(old) = positions.get_str_mut(&key) {
            if p.duplicates != "merge_dup" {
                return Err("duplicate merged models; use -d merge_dup".into());
            }
            for (id, _) in c.members.iter() {
                if !old.members.contains(id) {
                    old.members.insert(id.clone(), 1);
                    old.order.push(id.clone());
                }
            }
            let mut new = collapse(keys(&old.members), all, p.ends == "longest_ends")?;
            new.order = old.order.clone();
            new.members = old.members.clone();
            let mut ss = old.start_support.clone();
            ss.extend(new.start_support);
            new.start_support = ss;
            let mut es = old.end_support.clone();
            es.extend(new.end_support);
            new.end_support = es;
            *old = new;
        } else {
            positions.insert(key, c);
        }
    }
    let keys: Vec<_> = positions.iter().map(|(k, _)| k.clone()).collect();
    sort_position_keys(&keys)?
        .iter()
        .map(|k| {
            positions
                .remove(k)
                .ok_or_else(|| format!("missing sorted merge {k}"))
        })
        .collect()
}

struct BedStyle<'a> {
    cds: (i64, i64),
    colour: &'a str,
}

fn bed(
    chrom: &str,
    strand: &str,
    starts: &[i64],
    ends: &[i64],
    name: &str,
    style: BedStyle<'_>,
    bounds: Option<(i64, i64)>,
) -> String {
    let BedStyle { cds, colour } = style;
    let (start, end) = bounds.unwrap_or((starts[0], *ends.last().unwrap()));
    let sizes: Vec<_> = starts.iter().zip(ends).map(|(a, b)| b - a).collect();
    let rels: Vec<_> = starts.iter().map(|a| a - start).collect();
    format!(
        "{chrom}\t{start}\t{end}\t{name}\t40\t{strand}\t{}\t{}\t{colour}\t{}\t{}\t{}\n",
        cds.0,
        cds.1,
        starts.len(),
        ints(&sizes),
        ints(&rels)
    )
}
const COLOURS: [&str; 10] = [
    "255,0,0",
    "255,100,0",
    "255,200,0",
    "200,255,0",
    "0,255,200",
    "0,200,255",
    "0,100,255",
    "0,0,255",
    "100,0,255",
    "200,0,255",
];

pub fn read_merge_sources(path: &Path) -> Result<Vec<MergeSource>, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("read {}: {e}", path.display()))?;
    text.trim_end_matches('\n')
        .split('\n')
        .map(|line| {
            let cols: Vec<_> = line.split('\t').collect();
            if cols.len() != 4 || line.contains(['\r', ' ']) {
                return Err(
                    "merge filelist must have four tab-separated columns, without spaces or blank rows".into(),
                );
            }
            Ok(MergeSource {
                path: cols[0].into(),
                cap: cols[1].into(),
                priority: cols[2].into(),
                name: cols[3].into(),
            })
        })
        .collect()
}

pub fn merge_sources(sources: &[MergeSource], p: &MergeSettings) -> Result<MergeTexts, String> {
    if !["common_ends", "longest_ends"].contains(&p.ends.as_str()) {
        return Err("-e must be common_ends or longest_ends".into());
    }
    if !["no_merge", "merge_dup"].contains(&p.duplicates.as_str()) {
        return Err("-d must be no_merge or merge_dup".into());
    }
    let mut info: HashMap<String, (bool, [i64; 3])> = HashMap::new();
    let mut rows = Vec::new();
    let mut gene_ids: HashMap<(String, String), String> = HashMap::new();
    for source in sources {
        if !info.contains_key(&source.name) {
            if source.cap != "capped" && source.cap != "no_cap" {
                return Err("source cap must be capped or no_cap".into());
            }
            let priority = source
                .priority
                .split(',')
                .map(|s| {
                    s.parse::<i64>()
                        .map_err(|_| "bad source priority".to_string())
                })
                .collect::<Result<Vec<_>, _>>()?;
            let priority: [i64; 3] = priority
                .try_into()
                .map_err(|_| "source priority requires three integers")?;
            info.insert(source.name.clone(), (source.cap == "capped", priority));
        }
        let (capped, priority) = info[&source.name];
        let text = fs::read_to_string(&source.path)
            .map_err(|e| format!("read {}: {e}", source.path.display()))?;
        for line in text.trim_end_matches('\n').split('\n') {
            let c: Vec<_> = line.split('\t').collect();
            if c.len() != 12 {
                return Err(format!("{}: expected BED12", source.path.display()));
            }
            let num = |s: &str| {
                s.trim()
                    .parse::<i64>()
                    .map_err(|_| format!("bad BED integer {s}"))
            };
            let (start, end) = (num(c[1])?, num(c[2])?);
            let ids: Vec<_> = c[3].split(';').collect();
            if ids.len() < 2 {
                return Err("BED name must be gene_id;transcript_id".into());
            }
            if c[5] != "+" && c[5] != "-" {
                return Err("BED strand must be + or -".into());
            }
            let sizes = c[10]
                .trim_end_matches(',')
                .split(',')
                .map(num)
                .collect::<Result<Vec<_>, _>>()?;
            let rels = c[11]
                .trim_end_matches(',')
                .split(',')
                .map(num)
                .collect::<Result<Vec<_>, _>>()?;
            if sizes.is_empty() || sizes.len() != rels.len() || sizes.len() as i64 != num(c[9])? {
                return Err("inconsistent BED block counts".into());
            }
            let starts: Vec<_> = rels.iter().map(|x| start + x).collect();
            let ends = starts.iter().zip(sizes).map(|(a, n)| a + n).collect();
            let m = Model {
                id: format!("{}_{}", source.name, ids[1]),
                gene: ids[0].into(),
                trans: ids[1].into(),
                source: source.name.clone(),
                chrom: c[0].into(),
                strand: c[5].into(),
                start,
                end,
                starts,
                ends,
                cds_start: num(c[6])?,
                cds_end: num(c[7])?,
                capped,
                priority,
            };
            gene_ids.insert((m.source.clone(), m.trans.clone()), m.gene.clone());
            rows.push(m);
        }
    }
    if p.source_id != "no_source_id" {
        for name in p.source_id.split(',') {
            if !info.contains_key(name) {
                return Err(format!("unknown source for -s: {name}"));
            }
        }
    }
    let mut source_names: Vec<_> = info.keys().cloned().collect();
    source_names.sort();
    let colours: HashMap<_, _> = source_names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.clone(), COLOURS[i.min(9)]))
        .collect();
    rows.sort_by(|a, b| (&a.chrom, a.start, a.end).cmp(&(&b.chrom, b.start, b.end)));
    let mut out=MergeTexts{trans_report:"transcript_id\tnum_clusters\tsources\tstart_wobble_list\tend_wobble_list\texon_start_support\texon_end_support\tall_source_trans\n".into(),
        gene_report:"gene_id\tnum_clusters\tnum_final_trans\tsources\tchrom\tstart\tend\tsource_genes\tsource_summary\n".into(),..MergeTexts::default()};
    let mut gene_count = 0;
    let mut pos = 0;
    while pos < rows.len() {
        let mut end = rows[pos].end;
        let mut stop = pos + 1;
        while stop < rows.len() && rows[stop].chrom == rows[pos].chrom && rows[stop].start < end {
            end = end.max(rows[stop].end);
            stop += 1;
        }
        let locus = &rows[pos..stop];
        let mut all = HashMap::new();
        let mut forward = Vec::new();
        let mut reverse = Vec::new();
        for m in locus {
            if all.insert(m.id.clone(), m.clone()).is_some() {
                return Err(format!("duplicate source transcript {} in locus", m.id));
            }
            if m.strand == "+" {
                forward.push(m.id.clone());
            } else {
                reverse.push(m.id.clone());
            }
        }
        let (fg, fs) = gene_group(&forward, &all)?;
        let (rg, rs) = gene_group(&reverse, &all)?;
        let starts: BTreeSet<_> = fs.into_iter().chain(rs).collect();
        for start in starts {
            for dict in [&fg, &rg] {
                if let Some(ids) = dict.get(&start) {
                    gene_count += 1;
                    let gene = format!("G{gene_count}");
                    let ms: Vec<_> = ids.iter().map(|(id, _)| &all[id]).collect();
                    let merged = collapse_groups(simplify(&ms, p), &all, p)?;
                    let count = merged.len();
                    let mut tracks = BTreeSet::new();
                    for (i, c) in merged.into_iter().enumerate() {
                        let id = format!("{gene}.{}", i + 1);
                        let first = &all[&c.order[0]];
                        let sources: BTreeSet<_> = c
                            .members
                            .iter()
                            .map(|(id, _)| all[id].source.clone())
                            .collect();
                        let source_line = sources.iter().cloned().collect::<Vec<_>>().join(",");
                        out.trans_report.push_str(&format!(
                            "{id}\t{}\t{source_line}\t{}\t{}\t{}\t{}\t{}\n",
                            c.members.len(),
                            ints(&c.sw),
                            ints(&c.ew),
                            c.start_support.join(";"),
                            c.end_support.join(";"),
                            c.order.join(",")
                        ));
                        let mut name = format!("{gene};{id}");
                        let mut cds = (0, 0);
                        for (mid, _) in c.members.iter() {
                            let m = &all[mid];
                            let track = mid.rsplit_once('.').map_or(mid.as_str(), |(gene, _)| gene);
                            tracks.insert(track.to_string());
                            out.merge.push_str(&bed(
                                &m.chrom,
                                &m.strand,
                                &m.starts,
                                &m.ends,
                                &format!("{id};{mid}"),
                                BedStyle {
                                    cds: (m.start, m.end),
                                    colour: "255,0,0",
                                },
                                Some((m.start, m.end)),
                            ));
                            let (source, trans) = mid.split_once('_').unwrap();
                            if p.source_id != "no_source_id" {
                                for flag in p.source_id.split(',') {
                                    if source == flag {
                                        let g = gene_ids
                                            .get(&(source.into(), trans.into()))
                                            .ok_or("source transcript ID lookup failed")?;
                                        name.push_str(&format!(";{g};{trans}"));
                                    }
                                }
                            }
                            if p.cds != "no_cds" {
                                for flag in p.cds.split(',') {
                                    if source == flag {
                                        cds = (m.cds_start, m.cds_end);
                                    }
                                }
                            }
                        }
                        if cds.0 == 0 {
                            cds = (c.starts[0], *c.ends.last().unwrap());
                        }
                        let colour = if sources.len() > 1 {
                            COLOURS[9]
                        } else {
                            colours[sources.first().unwrap()]
                        };
                        out.bed.push_str(&bed(
                            &first.chrom,
                            &first.strand,
                            &c.starts,
                            &c.ends,
                            &name,
                            BedStyle { cds, colour },
                            None,
                        ));
                        out.transcripts.push(Transcript {
                            chrom: first.chrom.clone(),
                            strand: if first.strand == "+" {
                                Strand::Forward
                            } else {
                                Strand::Reverse
                            },
                            exons: c
                                .starts
                                .iter()
                                .zip(&c.ends)
                                .map(|(a, b)| Exon {
                                    start: a + 1,
                                    end: b + 1,
                                })
                                .collect(),
                            gene_id: gene.clone(),
                            transcript_id: id,
                            source: if sources.len() == 1 {
                                sources.first().cloned()
                            } else {
                                None
                            },
                            score: None,
                            cds_start: Some(cds.0 + 1),
                            cds_end: Some(cds.1 + 1),
                        });
                    }
                    let source_line = ms
                        .iter()
                        .map(|m| m.source.clone())
                        .collect::<BTreeSet<_>>()
                        .into_iter()
                        .collect::<Vec<_>>()
                        .join(",");
                    let mut summary: BTreeMap<String, usize> = BTreeMap::new();
                    for track in &tracks {
                        *summary
                            .entry(track.split('_').next().unwrap().into())
                            .or_default() += 1;
                    }
                    out.gene_report.push_str(&format!(
                        "{gene}\t{}\t{count}\t{source_line}\t{}\t{}\t{}\t{}\t{}\n",
                        ms.len(),
                        ms[0].chrom,
                        ms.iter().map(|m| m.start).min().unwrap(),
                        ms.iter().map(|m| m.end).max().unwrap(),
                        tracks.into_iter().collect::<Vec<_>>().join(","),
                        summary
                            .iter()
                            .map(|(k, v)| format!("{k}:{v}"))
                            .collect::<Vec<_>>()
                            .join(",")
                    ));
                }
            }
        }
        pos = stop;
    }
    Ok(out)
}

pub fn write_merge_texts(prefix: &Path, texts: &MergeTexts) -> Result<(), String> {
    for (suffix, body) in [
        (".bed", &texts.bed),
        ("_trans_report.txt", &texts.trans_report),
        ("_gene_report.txt", &texts.gene_report),
        ("_merge.txt", &texts.merge),
    ] {
        let path = suffix_path(prefix, suffix);
        fs::write(&path, body).map_err(|e| format!("write {}: {e}", path.display()))?;
    }
    Ok(())
}
