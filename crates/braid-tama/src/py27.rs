//! CPython 2.7 dict and string hash with `PYTHONHASHSEED=0`.
//!
//! Iteration follows table-slot order, not insertion order. The probe and
//! resize rules are the ones in `dictobject.c` at tag v2.7.18: open addressing,
//! `PERTURB_SHIFT` 5, growth at two-thirds fill, and a new table of the
//! smallest power of two above `4 * used` (or `2 * used` past 50_000).
//! Deletion leaves a dummy slot and does not resize.
//!
//! This stays inside `braid-tama`. Other algorithms must not depend on it.

const MINSIZE: usize = 8;
const PERTURB_SHIFT: u32 = 5;

/// `string_hash` from CPython 2.7 `stringobject.c`, secret prefix and suffix 0.
/// Bytes, not Unicode scalar values. Empty string hashes to 0. A result of -1
/// is rewritten to -2, matching the interpreter.
pub fn py27_str_hash(s: &str) -> i64 {
    let bytes = s.as_bytes();
    let len = bytes.len();
    if len == 0 {
        return 0;
    }
    let mut x: i64 = (bytes[0] as i64) << 7;
    let mut i = 0;
    let mut remaining = len as isize;
    while {
        remaining -= 1;
        remaining >= 0
    } {
        x = x.wrapping_mul(1_000_003) ^ (bytes[i] as i64);
        i += 1;
    }
    x ^= len as i64;
    if x == -1 {
        -2
    } else {
        x
    }
}

/// Python 2 `int` hash on a 64-bit build: the value itself, except -1 becomes -2.
pub fn py27_int_hash(value: i64) -> i64 {
    if value == -1 {
        -2
    } else {
        value
    }
}

pub trait Py27Key: Eq {
    fn py_hash(&self) -> i64;
}

impl Py27Key for String {
    fn py_hash(&self) -> i64 {
        py27_str_hash(self)
    }
}

impl Py27Key for i64 {
    fn py_hash(&self) -> i64 {
        py27_int_hash(*self)
    }
}

#[derive(Clone, Debug)]
enum Slot<K, V> {
    Empty,
    Dummy,
    Live { hash: i64, key: K, value: V },
}

#[derive(Clone, Debug)]
pub struct Py27Dict<K, V> {
    table: Vec<Slot<K, V>>,
    mask: u64,
    used: usize,
    fill: usize,
}

impl<K, V> Default for Py27Dict<K, V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K, V> Py27Dict<K, V> {
    pub fn new() -> Self {
        Self {
            table: empty_table(MINSIZE),
            mask: (MINSIZE as u64) - 1,
            used: 0,
            fill: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.used
    }

    pub fn is_empty(&self) -> bool {
        self.used == 0
    }

    pub fn capacity(&self) -> usize {
        self.table.len()
    }

    pub fn iter(&self) -> Iter<'_, K, V> {
        Iter {
            slots: self.table.iter(),
        }
    }

    fn live_at(&self, index: usize) -> Option<(&K, &V)> {
        match &self.table[index] {
            Slot::Live { key, value, .. } => Some((key, value)),
            _ => None,
        }
    }
}

impl<K: Py27Key, V> Py27Dict<K, V> {
    pub fn get(&self, key: &K) -> Option<&V> {
        let index = self.probe(key.py_hash(), |found| found == key);
        self.live_at(index).map(|(_, value)| value)
    }

    pub fn get_mut(&mut self, key: &K) -> Option<&mut V> {
        let index = self.probe(key.py_hash(), |found| found == key);
        match &mut self.table[index] {
            Slot::Live { value, .. } => Some(value),
            _ => None,
        }
    }

    pub fn contains(&self, key: &K) -> bool {
        self.get(key).is_some()
    }

    pub fn insert(&mut self, key: K, value: V) -> Option<V> {
        let hash = key.py_hash();
        let index = self.probe(hash, |found| found == &key);
        match &mut self.table[index] {
            Slot::Live { value: slot, .. } => {
                return Some(std::mem::replace(slot, value));
            }
            Slot::Empty => {
                self.fill += 1;
                self.used += 1;
            }
            Slot::Dummy => {
                self.used += 1;
            }
        }
        self.table[index] = Slot::Live { hash, key, value };
        let size = self.mask as usize + 1;
        if self.fill * 3 >= size * 2 {
            let factor = if self.used > 50_000 { 2 } else { 4 };
            self.resize(factor * self.used);
        }
        None
    }

    /// Insert `key` when it is absent. `K: Clone` lets the lookup survive a resize.
    pub fn or_insert_with<F>(&mut self, key: K, make: F) -> &mut V
    where
        F: FnOnce() -> V,
        K: Clone,
    {
        if self.get(&key).is_none() {
            self.insert(key.clone(), make());
        }
        self.get_mut(&key).expect("key present after insert")
    }

    pub fn remove(&mut self, key: &K) -> Option<V> {
        let index = self.probe(key.py_hash(), |found| found == key);
        if !matches!(self.table[index], Slot::Live { .. }) {
            return None;
        }
        match std::mem::replace(&mut self.table[index], Slot::Dummy) {
            Slot::Live { value, .. } => {
                self.used -= 1;
                Some(value)
            }
            _ => None,
        }
    }

    fn probe(&self, hash: i64, eq: impl Fn(&K) -> bool) -> usize {
        let mask = self.mask;
        let mut i = (hash as u64) & mask;
        let mut freeslot: Option<usize> = None;
        match &self.table[i as usize] {
            Slot::Empty => return i as usize,
            Slot::Dummy => freeslot = Some(i as usize),
            Slot::Live {
                hash: found, key, ..
            } => {
                if *found == hash && eq(key) {
                    return i as usize;
                }
            }
        }
        let mut perturb = hash as u64;
        loop {
            i = i.wrapping_mul(5).wrapping_add(perturb).wrapping_add(1);
            perturb >>= PERTURB_SHIFT;
            let index = (i & mask) as usize;
            match &self.table[index] {
                Slot::Empty => return freeslot.unwrap_or(index),
                Slot::Dummy => {
                    if freeslot.is_none() {
                        freeslot = Some(index);
                    }
                }
                Slot::Live {
                    hash: found, key, ..
                } => {
                    if *found == hash && eq(key) {
                        return index;
                    }
                }
            }
        }
    }

    fn resize(&mut self, minused: usize) {
        let mut newsize = MINSIZE;
        while newsize <= minused {
            newsize <<= 1;
            assert!(newsize > 0, "Py27Dict size overflow");
        }
        if newsize == MINSIZE && self.table.len() == MINSIZE && self.fill == self.used {
            return;
        }
        let old = std::mem::replace(&mut self.table, empty_table(newsize));
        self.mask = (newsize as u64) - 1;
        self.used = 0;
        self.fill = 0;
        for slot in old {
            if let Slot::Live { hash, key, value } = slot {
                self.insert_clean(hash, key, value);
            }
        }
    }

    fn insert_clean(&mut self, hash: i64, key: K, value: V) {
        let mask = self.mask;
        let mut i = (hash as u64) & mask;
        if matches!(self.table[i as usize], Slot::Empty) {
            self.place(i as usize, hash, key, value);
            return;
        }
        let mut perturb = hash as u64;
        loop {
            i = i.wrapping_mul(5).wrapping_add(perturb).wrapping_add(1);
            perturb >>= PERTURB_SHIFT;
            let index = (i & mask) as usize;
            if matches!(self.table[index], Slot::Empty) {
                self.place(index, hash, key, value);
                return;
            }
        }
    }

    fn place(&mut self, index: usize, hash: i64, key: K, value: V) {
        self.table[index] = Slot::Live { hash, key, value };
        self.fill += 1;
        self.used += 1;
    }
}

impl<V> Py27Dict<String, V> {
    pub fn get_str(&self, key: &str) -> Option<&V> {
        let index = self.probe(py27_str_hash(key), |found| found == key);
        self.live_at(index).map(|(_, value)| value)
    }

    pub fn get_str_mut(&mut self, key: &str) -> Option<&mut V> {
        let index = self.probe(py27_str_hash(key), |found| found == key);
        match &mut self.table[index] {
            Slot::Live { value, .. } => Some(value),
            _ => None,
        }
    }
}

fn empty_table<K, V>(size: usize) -> Vec<Slot<K, V>> {
    let mut table = Vec::with_capacity(size);
    for _ in 0..size {
        table.push(Slot::Empty);
    }
    table
}

pub struct Iter<'a, K, V> {
    slots: std::slice::Iter<'a, Slot<K, V>>,
}

impl<'a, K, V> Iterator for Iter<'a, K, V> {
    type Item = (&'a K, &'a V);

    fn next(&mut self) -> Option<Self::Item> {
        for slot in self.slots.by_ref() {
            if let Slot::Live { key, value, .. } = slot {
                return Some((key, value));
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use std::path::Path;

    #[derive(Clone, Debug)]
    enum Scalar {
        Int(i64),
        Str(String),
    }

    fn parity(name: &str) -> String {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/parity")
            .join(name);
        std::fs::read_to_string(&path)
            .unwrap_or_else(|err| panic!("read {}: {err}", path.display()))
    }

    fn parse_seq(text: &str) -> Vec<Scalar> {
        let bytes = text.as_bytes();
        let mut i = 0;
        let mut out = Vec::new();
        while i < bytes.len() {
            if bytes[i] == b'\'' {
                i += 1;
                let start = i;
                while bytes[i] != b'\'' {
                    i += 1;
                }
                out.push(Scalar::Str(text[start..i].to_string()));
                i += 1;
            } else {
                let start = i;
                while i < bytes.len() && bytes[i] != b',' {
                    i += 1;
                }
                let token = text[start..i].trim();
                if !token.is_empty() {
                    out.push(Scalar::Int(
                        token.parse().unwrap_or_else(|_| panic!("int {token}")),
                    ));
                }
            }
            if i < bytes.len() && bytes[i] == b',' {
                i += 1;
            }
        }
        out
    }

    fn repr_scalar(value: &Scalar) -> String {
        match value {
            Scalar::Int(n) => n.to_string(),
            Scalar::Str(s) => format!("'{s}'"),
        }
    }

    struct Case {
        insert: Vec<Scalar>,
        order: String,
        hashes: Vec<i64>,
    }

    fn order_of(keys: &[Scalar]) -> String {
        match keys.first() {
            Some(Scalar::Str(_)) => {
                let mut dict = Py27Dict::new();
                for key in keys {
                    let Scalar::Str(text) = key else {
                        unreachable!()
                    };
                    dict.insert(text.clone(), 1);
                }
                dict.iter()
                    .map(|(key, _)| format!("'{key}'"))
                    .collect::<Vec<_>>()
                    .join(",")
            }
            Some(Scalar::Int(_)) => {
                let mut dict = Py27Dict::new();
                for key in keys {
                    let Scalar::Int(n) = key else { unreachable!() };
                    dict.insert(*n, 1);
                }
                dict.iter()
                    .map(|(key, _)| key.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            }
            None => String::new(),
        }
    }

    #[test]
    fn sixth_insert_grows_from_8_to_32() {
        let mut dict = Py27Dict::new();
        for n in 0..5 {
            dict.insert(n, ());
        }
        assert_eq!(dict.capacity(), 8);
        dict.insert(5, ());
        assert_eq!(dict.capacity(), 32);
        assert_eq!(dict.len(), 6);
    }

    #[test]
    fn string_and_int_hashes_match_seed0() {
        assert_eq!(py27_str_hash(""), 0);
        assert_eq!(py27_int_hash(-1), -2);
        assert_eq!(py27_int_hash(0), 0);
        assert_eq!(py27_int_hash(2147483647), 2147483647);
    }

    #[test]
    fn matches_seed0_probe() {
        let text = parity("dict_seed0.txt");
        let mut cases: HashMap<String, Case> = HashMap::new();
        let mut update = String::new();
        let mut popped = String::new();
        let mut nest_pri = String::new();
        let mut nest_coord = String::new();
        let mut nest_err = String::new();
        for line in text.lines() {
            let mut cols = line.split('\t');
            let kind = cols.next().unwrap();
            match kind {
                "VERSION" | "HASHSEED" => {}
                "CASE" => {
                    let label = cols.next().unwrap().to_string();
                    let _same = cols.next().unwrap();
                    let insert = parse_seq(cols.next().unwrap());
                    let order = cols.next().unwrap().to_string();
                    cases.insert(
                        label,
                        Case {
                            insert,
                            order,
                            hashes: Vec::new(),
                        },
                    );
                }
                "HASH" => {
                    let label = cols.next().unwrap();
                    let hashes: Vec<i64> = cols
                        .next()
                        .unwrap()
                        .split(',')
                        .map(|n| n.parse().unwrap())
                        .collect();
                    cases.get_mut(label).unwrap().hashes = hashes;
                }
                "UPDATE" => update = cols.next().unwrap().to_string(),
                "POP" => popped = cols.next().unwrap().to_string(),
                "NEST_PRI" => nest_pri = cols.next().unwrap().to_string(),
                "NEST_COORD_1" => nest_coord = cols.next().unwrap().to_string(),
                "NEST_ERR" => nest_err = cols.next().unwrap().to_string(),
                other => panic!("unexpected probe line {other}"),
            }
        }

        for (label, case) in &cases {
            assert!(!case.hashes.is_empty(), "{label} missing hashes");
            for (key, hash) in case.insert.iter().zip(&case.hashes) {
                let got = match key {
                    Scalar::Str(text) => py27_str_hash(text),
                    Scalar::Int(n) => py27_int_hash(*n),
                };
                assert_eq!(got, *hash, "{label} hash of {}", repr_scalar(key));
            }
            assert_eq!(order_of(&case.insert), case.order, "{label} iter order");
        }

        let read_ids = &cases["read_ids"].insert;
        let mut updated = Py27Dict::new();
        for key in read_ids {
            let Scalar::Str(text) = key else {
                unreachable!()
            };
            updated.insert(text.clone(), 1);
        }
        let Scalar::Str(first) = &read_ids[0] else {
            unreachable!()
        };
        updated.insert(first.clone(), 2);
        updated.insert("new_id".to_string(), 1);
        let got = updated
            .iter()
            .map(|(key, _)| format!("'{key}'"))
            .collect::<Vec<_>>()
            .join(",");
        assert_eq!(got, update);

        let mut removed = Py27Dict::new();
        for key in read_ids {
            let Scalar::Str(text) = key else {
                unreachable!()
            };
            removed.insert(text.clone(), 1);
        }
        let Scalar::Str(second) = &read_ids[1] else {
            unreachable!()
        };
        assert!(removed.remove(&second.clone()).is_some());
        let got = removed
            .iter()
            .map(|(key, _)| format!("'{key}'"))
            .collect::<Vec<_>>()
            .join(",");
        assert_eq!(got, popped);

        let mut nested: Py27Dict<i64, Py27Dict<i64, Py27Dict<String, i32>>> = Py27Dict::new();
        for (priority, coord, err) in [
            (1, 100, "na"),
            (1, 100, "3.A"),
            (0, 90, "0"),
            (1, 110, "na"),
            (2, 80, "5I"),
        ] {
            nested.or_insert_with(priority, Py27Dict::new);
            nested
                .get_mut(&priority)
                .unwrap()
                .or_insert_with(coord, Py27Dict::new);
            nested
                .get_mut(&priority)
                .unwrap()
                .get_mut(&coord)
                .unwrap()
                .insert(err.to_string(), 1);
        }
        let pri = nested
            .iter()
            .map(|(key, _)| key.to_string())
            .collect::<Vec<_>>()
            .join(",");
        assert_eq!(pri, nest_pri);
        let coord = nested
            .get(&1)
            .unwrap()
            .iter()
            .map(|(key, _)| key.to_string())
            .collect::<Vec<_>>()
            .join(",");
        assert_eq!(coord, nest_coord);
        let err = nested
            .get(&1)
            .unwrap()
            .get(&100)
            .unwrap()
            .iter()
            .map(|(key, _)| format!("'{key}'"))
            .collect::<Vec<_>>()
            .join(",");
        assert_eq!(err, nest_err);
    }
}
