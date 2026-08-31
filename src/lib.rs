pub use intervaltree;
use regex::Regex;
use std::{
    borrow::Borrow,
    cmp::Reverse,
    collections::{BTreeMap, HashMap},
    fmt::Debug,
    hash::Hash,
    ops::{Deref, DerefMut},
};

static BASES: [u8; 4] = ['A' as u8, 'C' as u8, 'G' as u8, 'T' as u8];

pub struct UnitAndRepeats {
    unit_size: u8,
    min_repeats: u8,
}

impl UnitAndRepeats {
    pub fn new(unit_size: u8, min_repeats: u8) -> Self {
        Self {
            unit_size,
            min_repeats,
        }
    }

    pub fn build_finder_regrex(&self) -> HashMap<String, Regex> {
        let motifs = generate_motifs(self.unit_size);
        motifs
            .into_iter()
            .map(|motif| {
                let regex_str = format!("({}){{{},}}", motif, self.min_repeats);
                let reg = Regex::new(&regex_str).unwrap();
                (motif, reg)
            })
            .collect()
    }
}

fn generate_motifs(unit_size: u8) -> Vec<String> {
    let mut tracer = vec![];
    let mut motifs = vec![];
    generate_motif_core(
        unit_size as usize,
        &mut tracer,
        &mut motifs,
        0,
        unit_size as usize,
    );
    motifs
}

fn generate_motif_core(
    remain_base_num: usize,
    tracer: &mut Vec<u8>,
    result: &mut Vec<String>,
    eq_first_cnt: usize,
    tot_base_num: usize,
) {
    if remain_base_num == 0 {
        if eq_first_cnt == tot_base_num && tot_base_num > 1 {
            return;
        }
        result.push(String::from_utf8(tracer.clone()).unwrap());
        return;
    }

    for cur_base in BASES {
        tracer.push(cur_base);
        generate_motif_core(
            remain_base_num - 1,
            tracer,
            result,
            eq_first_cnt + if cur_base == tracer[0] { 1 } else { 0 },
            tot_base_num,
        );
        tracer.pop();
    }
}

#[derive(Debug)]
pub struct Region2Motif<T> {
    // (usize, usize) -> (start, end)
    value: HashMap<(usize, usize), T>,
}

impl<T> Region2Motif<T>
where
    T: Clone,
{
    pub fn to_interval_search_tree(&self) -> intervaltree::IntervalTree<usize, T> {
        intervaltree::IntervalTree::from_iter(
            self.value
                .iter()
                .map(|(key, value)| (key.0..key.1, value.clone())),
        )
    }

    pub fn flatten(&self) -> HashMap<usize, Vec<((usize, usize), T)>> {
        let mut result = HashMap::new();

        self.value.iter().for_each(|(key, value)| {
            (key.0..key.1).into_iter().for_each(|pos| {
                result
                    .entry(pos)
                    .or_insert(vec![])
                    .push(((key.0, key.1), value.clone()));
            });
        });

        result
    }
}

impl<T> Default for Region2Motif<T> {
    fn default() -> Self {
        Self {
            value: HashMap::new(),
        }
    }
}

impl<T> Deref for Region2Motif<T>
where
    T: Sized,
{
    type Target = HashMap<(usize, usize), T>;
    fn deref(&self) -> &Self::Target {
        &self.value
    }
}
impl<T> DerefMut for Region2Motif<T>
where
    T: Sized,
{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.value
    }
}

pub fn all_seq_hp_tr_finder<RegK, SeqK, SeqV, MotifT>(
    all_regs: &Vec<HashMap<RegK, Regex>>,
    seqs: &HashMap<SeqK, SeqV>,
) -> HashMap<SeqK, Region2Motif<MotifT>>
where
    RegK: std::borrow::Borrow<String>,
    SeqK: Clone + Eq + Hash,
    SeqV: std::borrow::Borrow<String>,
    MotifT: From<String> + Clone + Borrow<String> + Debug,
{
    let mut match_patterns = HashMap::new();

    seqs.iter()
        .map(|(seq_name, seq)| {
            let region2motif = single_seq_hp_tr_finder(all_regs, &mut match_patterns, seq.borrow());

            (seq_name.clone(), region2motif)
        })
        .collect()
}

pub fn single_seq_hp_tr_finder<RegK, MotifT>(
    all_regs: &Vec<HashMap<RegK, Regex>>,
    match_patterns: &mut HashMap<String, MotifT>,
    seq: &str,
) -> Region2Motif<MotifT>
where
    RegK: std::borrow::Borrow<String>,
    MotifT: From<String> + Clone + Borrow<String> + Debug,
{
    let mut region2motif = Region2Motif::default();
    all_regs.iter().for_each(|regs| {
        hp_tr_finder(regs, seq, &mut region2motif, match_patterns);
    });

    dedup_overlaps(region2motif)
}

/// Keep only non-overlapping annotations, longest-span first (ties: smaller
/// unit, then leftmost start). A TR region admits several period
/// representations (e.g. a 46 bp pure-AT run is both `(AT)23` and `(ATAT)11`);
/// we keep the single canonical one, discarding the rest entirely.
fn dedup_overlaps<Pat>(region2motif: Region2Motif<Pat>) -> Region2Motif<Pat>
where
    Pat: From<String> + Clone + Borrow<String>,
{
    let mut matches: Vec<_> = region2motif
        .value
        .into_iter()
        .map(|((start, end), pat)| {
            // pat is `(motif)N`; DNA motifs never contain parentheses.
            let motif: &String = pat.borrow();
            let motif = &motif[motif.find('(').unwrap() + 1..motif.rfind(')').unwrap()];
            (start, end, motif.len(), pat)
        })
        .collect();

    // span desc, unit asc, start asc — a strict total order, so the output is
    // fully deterministic regardless of the HashMap's iteration order.
    matches.sort_by_key(|(start, end, unit, _)| (Reverse(end - start), *unit, *start));

    // `start -> end` of the intervals accepted so far. Mutually non-overlapping,
    // so the starts are unique and the ends come out in the same order.
    let mut accepted: BTreeMap<usize, usize> = BTreeMap::new();
    let mut result = Region2Motif::default();
    for (start, end, _, pat) in matches {
        // Intervals are half-open, so one starting exactly at `end` merely
        // touches this candidate and cannot overlap it: bound the lookup by
        // start alone, exclusively. Accepted intervals are non-overlapping, so
        // of those the greatest start also has the greatest end — checking only
        // it is complete.
        let overlaps = accepted
            .range(..end)
            .next_back()
            .is_some_and(|(_, &other_end)| other_end > start);
        if !overlaps {
            accepted.insert(start, end);
            result.insert((start, end), pat);
        }
    }

    result
}

pub fn hp_tr_finder<RegK, Pat>(
    regs: &HashMap<RegK, Regex>,
    seq: &str,
    region2motif: &mut Region2Motif<Pat>,
    match_patterns: &mut HashMap<String, Pat>,
) where
    RegK: std::borrow::Borrow<String>,
    Pat: From<String> + Clone + Debug,
{
    for (motif, reg) in regs {
        for m in reg.find_iter(&seq) {
            let (s, e) = (m.start(), m.end());
            let match_pat = format!("({}){}", motif.borrow(), (e - s) / motif.borrow().len());
            if !match_patterns.contains_key(&match_pat) {
                let m_pat_ = match_pat.clone();
                match_patterns.insert(m_pat_, match_pat.clone().into());
            }
            let start_end = (m.start(), m.end());

            //
            if region2motif.contains_key(&start_end) {
                continue;
            }
            // assert!(
            //     !region2motif.contains_key(&start_end),
            //     "duplicated start,end. {:?}, exists_pat: {:?}, new_pat:{}",
            //     start_end, region2motif.get(&start_end).unwrap(), match_pat
            // );

            region2motif.insert(
                (m.start(), m.end()),
                match_patterns.get(&match_pat).unwrap().clone(),
            );
        }
    }

    // regions
}

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, sync::Arc};

    use crate::{
        Region2Motif, UnitAndRepeats, all_seq_hp_tr_finder, dedup_overlaps, generate_motifs,
    };

    #[test]
    fn test_generate_motifs() {
        let motifs = generate_motifs(2);
        assert_eq!(
            motifs,
            vec![
                "AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"
            ]
        );

        let motifs = generate_motifs(3);
        assert_eq!(
            motifs,
            vec![
                "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA",
                "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCG", "CCT", "CGA", "CGC",
                "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC",
                "GCG", "GCT", "GGA", "GGC", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG",
                "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG"
            ]
        );
    }

    #[test]
    fn test_unit_and_repeats() {
        let unit_and_repeats = UnitAndRepeats::new(2, 2);
        println!("{:?}", unit_and_repeats.build_finder_regrex());
    }

    #[test]
    fn test_tr_finder() {
        let all_regs = vec![
            UnitAndRepeats::new(1, 3).build_finder_regrex(),
            UnitAndRepeats::new(2, 3).build_finder_regrex(),
            UnitAndRepeats::new(3, 3).build_finder_regrex(),
            UnitAndRepeats::new(4, 3).build_finder_regrex(),
        ];

        let mut seqs = HashMap::new();
        seqs.insert("seq1".to_string(), "ACGTACGTAAACGT".to_string());
        seqs.insert("seq2".to_string(), "ACACACCCGCGCG".to_string());
        let res: HashMap<String, crate::Region2Motif<Arc<String>>> =
            all_seq_hp_tr_finder(&all_regs, &seqs);
        println!("{res:?}");

        res.iter().for_each(|(_key, value)| {
            let mut result = value.flatten().into_iter().collect::<Vec<_>>();
            result.sort_by_key(|v| v.0);
            println!("{:?}", result);
        });
    }

    #[test]
    fn test_dedup_overlaps_canonical_annotation() {
        let all_regs = vec![
            UnitAndRepeats::new(1, 3).build_finder_regrex(),
            UnitAndRepeats::new(2, 3).build_finder_regrex(),
            UnitAndRepeats::new(3, 3).build_finder_regrex(),
            UnitAndRepeats::new(4, 3).build_finder_regrex(),
        ];

        // A 46 bp pure-AT run is also (ATAT)11 [0,44), (TATA)11 [2,46),
        // (TA)22 [1,45), (ATA)15 [0,45); the longest annotation (AT)23 [0,46)
        // overlaps and discards all of them.
        let mut seqs = HashMap::new();
        seqs.insert("pure_at".to_string(), "AT".repeat(23));
        let res: HashMap<String, crate::Region2Motif<Arc<String>>> =
            all_seq_hp_tr_finder(&all_regs, &seqs);
        let regions = res.get("pure_at").unwrap();
        assert_eq!(
            regions.value.clone(),
            HashMap::from([((0usize, 46usize), Arc::new("(AT)23".to_string()))])
        );

        // A 13 bp alternating run is both (TA)6 [0,12) and (AT)6 [1,13);
        // same span and unit, so the leftmost start wins.
        let mut seqs = HashMap::new();
        seqs.insert("alt".to_string(), "TATATATATATAT".to_string());
        let res: HashMap<String, crate::Region2Motif<Arc<String>>> =
            all_seq_hp_tr_finder(&all_regs, &seqs);
        let regions = res.get("alt").unwrap();
        assert_eq!(
            regions.value.clone(),
            HashMap::from([((0usize, 12usize), Arc::new("(TA)6".to_string()))])
        );
    }

    /// Feeds `dedup_overlaps` a set of hand-built annotations and returns what
    /// survived, sorted by start, so a test can spell out exactly which regions
    /// are kept. Patterns are `(motif)N` as the real finder emits them, with the
    /// copy count matching the span.
    fn dedup_keep(entries: &[((usize, usize), &str)]) -> Vec<(usize, usize, String)> {
        let mut input: Region2Motif<Arc<String>> = Region2Motif::default();
        for ((start, end), pattern) in entries {
            input.insert((*start, *end), Arc::new(pattern.to_string()));
        }

        let mut kept: Vec<_> = dedup_overlaps(input)
            .value
            .into_iter()
            .map(|((start, end), pattern)| (start, end, pattern.to_string()))
            .collect();
        kept.sort();
        kept
    }

    #[test]
    fn test_touching_annotations_are_not_overlaps() {
        // Intervals are half-open: [0,6) and [6,16) share no base, so both must
        // survive. The long run is accepted first in both input orders (sorted
        // by span), and used to reject its left neighbour for starting exactly
        // where that neighbour ended.
        let both = vec![(0, 6, "(A)6".to_string()), (6, 16, "(C)10".to_string())];
        assert_eq!(dedup_keep(&[((0, 6), "(A)6"), ((6, 16), "(C)10")]), both);
        assert_eq!(dedup_keep(&[((6, 16), "(C)10"), ((0, 6), "(A)6")]), both);

        // A chain of mutually touching runs — each one abuts an already-accepted
        // neighbour with a larger span.
        assert_eq!(
            dedup_keep(&[((0, 10), "(A)10"), ((10, 25), "(C)15"), ((25, 35), "(G)10"),]),
            vec![
                (0, 10, "(A)10".to_string()),
                (10, 25, "(C)15".to_string()),
                (25, 35, "(G)10".to_string()),
            ]
        );
    }

    #[test]
    fn test_real_overlaps_are_still_dropped() {
        // The touching-neighbour exemption must not blind the check to a
        // genuine overlap further left: [10,20) overlaps accepted [0,12) even
        // though [20,31) — a larger key — only touches it at 20.
        assert_eq!(
            dedup_keep(&[((0, 12), "(AT)6"), ((20, 31), "(C)11"), ((10, 20), "(AC)5")]),
            vec![(0, 12, "(AT)6".to_string()), (20, 31, "(C)11".to_string())]
        );

        // Nested and partial overlaps still lose to the longest span.
        assert_eq!(
            dedup_keep(&[
                ((0, 16), "(AT)8"),
                ((0, 6), "(AT)3"),
                ((5, 15), "(TA)5"),
                ((4, 12), "(ATAT)2"),
            ]),
            vec![(0, 16, "(AT)8".to_string())]
        );
    }

    #[test]
    fn test_adjacent_trs_survive_the_full_pipeline() {
        let all_regs = vec![
            UnitAndRepeats::new(1, 3).build_finder_regrex(),
            UnitAndRepeats::new(2, 3).build_finder_regrex(),
            UnitAndRepeats::new(3, 3).build_finder_regrex(),
            UnitAndRepeats::new(4, 3).build_finder_regrex(),
        ];

        // Three homopolymer runs with no filler between them: (A)6|(C)10|(T)5.
        let mut seqs = HashMap::new();
        seqs.insert(
            "abutting_runs".to_string(),
            format!("{}{}{}", "A".repeat(6), "C".repeat(10), "T".repeat(5)),
        );
        let res: HashMap<String, crate::Region2Motif<Arc<String>>> =
            all_seq_hp_tr_finder(&all_regs, &seqs);

        let mut kept: Vec<_> = res["abutting_runs"]
            .value
            .iter()
            .map(|((start, end), pattern)| (*start, *end, pattern.to_string()))
            .collect();
        kept.sort();
        assert_eq!(
            kept,
            vec![
                (0, 6, "(A)6".to_string()),
                (6, 16, "(C)10".to_string()),
                (16, 21, "(T)5".to_string()),
            ]
        );
    }
}
