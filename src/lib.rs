use intervaltree;
use regex::Regex;
use std::{
    collections::HashMap,
    hash::Hash,
    ops::{Deref, DerefMut},
    rc::Rc,
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

#[derive(Debug, Default)]
pub struct Region2Motif {
    // (usize, usize) -> (start, end)
    value: HashMap<(usize, usize), Rc<String>>,
}

impl Region2Motif {
    pub fn to_interval_search_tree(&self) -> intervaltree::IntervalTree<usize, Rc<String>> {
        intervaltree::IntervalTree::from_iter(
            self.value
                .iter()
                .map(|(key, value)| (key.0..key.1, value.clone())),
        )
    }
}

impl Deref for Region2Motif {
    type Target = HashMap<(usize, usize), Rc<String>>;
    fn deref(&self) -> &Self::Target {
        &self.value
    }
}
impl DerefMut for Region2Motif {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.value
    }
}

pub fn all_seq_hp_tr_finder<SeqK, SeqV>(
    all_regs: &Vec<HashMap<String, Regex>>,
    seqs: &HashMap<SeqK, SeqV>,
) -> HashMap<SeqK, Region2Motif>
where
    SeqK: Clone + Eq + Hash,
    SeqV: std::borrow::Borrow<str>,
{
    let mut match_patterns = HashMap::new();

    seqs.iter()
        .map(|(seq_name, seq)| {
            let region2motif = single_seq_hp_tr_finder(all_regs, &mut match_patterns, seq.borrow());

            (seq_name.clone(), region2motif)
        })
        .collect()
}

pub fn single_seq_hp_tr_finder(
    all_regs: &Vec<HashMap<String, Regex>>,
    match_patterns: &mut HashMap<String, Rc<String>>,
    seq: &str,
) -> Region2Motif {
    let mut region2motif = Region2Motif::default();
    all_regs.iter().for_each(|regs| {
        hp_tr_finder(regs, seq, &mut region2motif, match_patterns);
    });

    region2motif
}

pub fn hp_tr_finder(
    regs: &HashMap<String, Regex>,
    seq: &str,
    region2motif: &mut Region2Motif,
    match_patterns: &mut HashMap<String, Rc<String>>,
) {
    for (motif, reg) in regs {
        for m in reg.find_iter(&seq) {
            let (s, e) = (m.start(), m.end());
            let match_pat = format!("({}){}", motif, (e - s) / motif.len());
            if !match_patterns.contains_key(&match_pat) {
                let m_pat_ = match_pat.clone();
                match_patterns.insert(m_pat_, Rc::new(match_pat.clone()));
            }
            let start_end = (m.start(), m.end());
            assert!(
                !region2motif.contains_key(&start_end),
                "duplicated start,end. {:?}",
                start_end
            );

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
    use std::collections::HashMap;

    use crate::{UnitAndRepeats, all_seq_hp_tr_finder, generate_motifs};

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
        let res = all_seq_hp_tr_finder(&all_regs, &seqs);
        println!("{res:?}");
    }
}
