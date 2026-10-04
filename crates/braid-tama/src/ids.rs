//! Transcript ID collections whose traversal follows CPython 2.7.
use crate::Py27Dict;

pub(crate) type Ids = Py27Dict<String, i32>;
pub(crate) type GeneGroups = Py27Dict<i64, Ids>;

pub(crate) fn keys(ids: &Ids) -> Vec<String> {
    ids.iter().map(|(id, _)| id.clone()).collect()
}

/// Reinsertion rebuilds the hash table and can change its traversal order.
/// A clone would preserve the original table and is not equivalent.
pub(crate) fn reinsert(ids: &Ids) -> Ids {
    let mut result = Ids::new();
    for (id, _) in ids.iter() {
        result.insert(id.clone(), 1);
    }
    result
}
