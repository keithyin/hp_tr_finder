pub use intervaltree;
use regex::Regex;
use std::{
    borrow::Borrow,
    collections::HashMap,
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
    generate_motif_core(unit_size as usize, &mut tracer, &mut motifs);
    motifs
}

fn generate_motif_core(remain_base_num: usize, tracer: &mut Vec<u8>, result: &mut Vec<String>) {
    if remain_base_num == 0 {
        let motif = String::from_utf8(tracer.clone()).unwrap();
        // 过滤非最简（non-primitive）motif：本身是更短单位串联重复的 motif
        // （含同聚物，如 AAAA、GTGT），避免同一 TR 被长短 unit 重复报告
        if tracer.len() > 1 && !is_primitive_motif(&motif) {
            return;
        }
        result.push(motif);
        return;
    }

    for cur_base in BASES {
        tracer.push(cur_base);
        generate_motif_core(remain_base_num - 1, tracer, result);
        tracer.pop();
    }
}

/// motif 是最简（primitive）当且仅当它不是任何更短前缀的串联重复。
/// 例如 AAAA（周期 1）、GTGT（周期 2）都是非最简，ACGT、ACGAC 是最简
/// （后者周期 3 不整除长度 5，不得误过滤）。
pub fn is_primitive_motif(motif: &str) -> bool {
    let bytes = motif.as_bytes();
    let n = bytes.len();
    (1..n).all(|d| {
        if !n.is_multiple_of(d) {
            return true;
        }
        (d..n).any(|i| bytes[i] != bytes[i % d])
    })
}

/// 返回 motif 的最小重复周期及拷贝数。GTGT -> (GT, 2)，GT -> (GT, 1)，
/// ACGAC 无整除周期 -> (ACGAC, 1)。空 motif 返回 ("", 1)（调用方应自行拒绝空
/// motif；此处的兜底只保证不会除零 panic）。
pub fn simplify_motif(motif: &str) -> (String, usize) {
    let bytes = motif.as_bytes();
    let n = bytes.len();
    if n == 0 {
        return (String::new(), 1);
    }
    let period = (1..=n)
        .find(|&d| n.is_multiple_of(d) && (d..n).all(|i| bytes[i] == bytes[i % d]))
        .unwrap_or(n);
    (motif[..period].to_string(), n / period)
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

    dedup_rotation_hits(region2motif)
}

/// tr-finder v2 旋转冗余去重：同周期 p 内，起点错位量 shift ∈ [0, p-1] 的命中
/// 视为同一 TR 的循环旋转冗余（如 (AT) 6-15 与 (TA) 7-16），保留起点最靠左的
/// 一条（同起点取 end 更大者）；shift ≥ p 一律保留，即使区间重叠（如边界共享型
/// 相邻 TR：AAT×2 0-6 与 TGG×2 5-11 是两个真实 TR，shift=5 ≥ p=3）。跨周期
/// 重叠（如 (AT) 与 (ATG)）是真实的不同注释，天然不受影响。shift=0 在真实数据
/// 中不可能出现（同起点同周期意味着同一 motif，而贪婪匹配下同 motif 命中互不
/// 重叠），仅作防御性 tie-break 处理。
///
/// 正确性依据（详见 tr-finder docs/rotation-dedup.md）：
/// 1. 同一 TR 的旋转命中起点必落在第一个命中的第一个重复单元内，故 shift ∈ [1, p-1]；
///    shift=0 意味着同一 motif，而贪婪匹配下同一 motif 命中互不重叠，矛盾。
/// 2. min_repeats ≥ 2 时，shift < p 的两个不同 motif 命中，后者的第一个完整单元
///    整体落在前者的周期区间内部，必为前者的旋转，即同一 TR。
///
/// 前置条件：motif 均 primitive（motif 长度 = 最小周期，由 generate_motifs /
/// simplify_motif 保证），min_repeats ≥ 2。
fn dedup_rotation_hits<Pat>(region2motif: Region2Motif<Pat>) -> Region2Motif<Pat>
where
    Pat: From<String> + Clone + Borrow<String> + Debug,
{
    let mut hits: Vec<_> = region2motif
        .value
        .into_iter()
        .map(|((start, end), pat)| {
            // pat 是 `(motif)N`；DNA motif 不含括号，period = motif 长度。
            let pat_str: &String = pat.borrow();
            let motif = &pat_str[pat_str.find('(').unwrap() + 1..pat_str.rfind(')').unwrap()];
            (start, end, motif.len(), pat)
        })
        .collect();

    // start 升序（同起点取 end 更大者）——保留起点最靠左的一条
    hits.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)));

    let mut kept_starts: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut result = Region2Motif::default();
    for (start, end, period, pat) in hits {
        // 每个 period 已保留命中的起点按 start 升序（恒 <= 当前命中起点）。
        // 同周期已保留命中的相邻起点间距恒 >= period，故 shift < period 的窗口内
        // 最多只有最后一个起点，只需与它比较（O(1)），无需线性扫描整个列表；
        // 保留起点列表也无需按周期清理：shift 比较天然跨不同 TR 位置安全。
        let starts = kept_starts.entry(period).or_default();
        if starts.last().is_some_and(|&s| start - s < period) {
            continue;
        }
        starts.push(start);
        result.insert((start, end), pat);
    }

    result
}

/// 将命中区间左右各扩区 left/right bp，边界 clamp 到 [0, len]。
pub fn expand_regions(regions: &mut [(usize, usize)], left: usize, right: usize, len: usize) {
    for region in regions {
        region.0 = region.0.saturating_sub(left).min(len);
        region.1 = (region.1 + right).min(len);
    }
}

/// 合并重叠或相邻的区间，返回按起点升序的合并结果。
pub fn merge_region(regions: &mut [(usize, usize)]) -> Vec<(usize, usize)> {
    let mut result = vec![];
    if regions.is_empty() {
        return result;
    }

    regions.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    result.push(regions[0]);

    for region in &regions[1..] {
        let last = result.last_mut().unwrap();
        if last.1 >= region.0 {
            last.1 = last.1.max(region.1);
        } else {
            result.push(*region);
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
        Region2Motif, UnitAndRepeats, all_seq_hp_tr_finder, dedup_rotation_hits, expand_regions,
        generate_motifs, is_primitive_motif, merge_region, simplify_motif,
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
    fn test_generate_motifs_filters_non_primitive() {
        // unit 4：非最简 motif（含同聚物）不应出现
        let motifs = generate_motifs(4);
        for m in ["AAAA", "CCCC", "GTGT", "ACAC", "ATAT"] {
            assert!(!motifs.contains(&m.to_string()), "{m} should be filtered");
        }
        // 最简 motif 应保留（ACGTAC 是最简 6-mer，只应在 unit 6 中出现）
        for m in ["GTAT", "ACGT", "AAGT"] {
            assert!(motifs.contains(&m.to_string()), "{m} should be kept");
        }
        assert!(!motifs.contains(&"ACGTAC".to_string()));
        assert!(generate_motifs(6).contains(&"ACGTAC".to_string()));

        // 周期不整除长度的 motif 不应被误过滤（ACGAC 周期 3，n=5）
        let motifs = generate_motifs(5);
        assert!(motifs.contains(&"ACGAC".to_string()));
        assert!(motifs.contains(&"AAAAT".to_string())); // 周期 3 不整除 5，仍是最简
        assert!(!motifs.contains(&"AAAAA".to_string())); // 同聚物被过滤

        // 过滤后 unit 模式不再生成重复正则：GTGTGTGTGTGT 只被 (GT) 命中一次
        let regs = UnitAndRepeats::new(2, 2).build_finder_regrex();
        assert!(regs.contains_key("GT") && !regs.contains_key("GTGT"));
        let hits: Vec<_> = regs["GT"].find_iter("GTGTGTGTGTGT").collect();
        assert_eq!(hits.len(), 1);
        assert_eq!((hits[0].start(), hits[0].end()), (0, 12));
    }

    #[test]
    fn test_is_primitive_motif() {
        assert!(!is_primitive_motif("AAAA"));
        assert!(!is_primitive_motif("GTGT"));
        assert!(!is_primitive_motif("ACAC"));
        assert!(is_primitive_motif("ACGT"));
        assert!(is_primitive_motif("ACGAC")); // 周期 3 不整除 5
        assert!(is_primitive_motif("A"));
    }

    #[test]
    fn test_simplify_motif() {
        assert_eq!(simplify_motif("GTGT"), ("GT".to_string(), 2));
        assert_eq!(simplify_motif("A"), ("A".to_string(), 1));
        assert_eq!(simplify_motif("AAAA"), ("A".to_string(), 4));
        assert_eq!(simplify_motif("ACGT"), ("ACGT".to_string(), 1));
        assert_eq!(simplify_motif("ACGAC"), ("ACGAC".to_string(), 1)); // 周期不整除，保持原样
        // 空 motif 不得除零 panic（CLI 层会先拒绝）
        assert_eq!(simplify_motif(""), (String::new(), 1));
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
    fn test_rotation_dedup_pipeline() {
        let all_regs = vec![
            UnitAndRepeats::new(1, 3).build_finder_regrex(),
            UnitAndRepeats::new(2, 3).build_finder_regrex(),
            UnitAndRepeats::new(3, 3).build_finder_regrex(),
            UnitAndRepeats::new(4, 3).build_finder_regrex(),
        ];

        // A 46 bp pure-AT run only yields (AT)[0,46) and its rotation
        // (TA)[1,45) — (ATAT)/(TATA) regexes no longer exist after the
        // primitive filter, and (ATA)/... never match a pure AT repeat
        // (their 2nd repeat unit is TAT). v2: (TA) is dropped (shift < p).
        let mut seqs = HashMap::new();
        seqs.insert("pure_at".to_string(), "AT".repeat(23));
        let res: HashMap<String, crate::Region2Motif<Arc<String>>> =
            all_seq_hp_tr_finder(&all_regs, &seqs);
        let mut kept: Vec<_> = res["pure_at"]
            .value
            .iter()
            .map(|((start, end), pat)| (*start, *end, pat.to_string()))
            .collect();
        kept.sort();
        assert_eq!(kept, vec![(0, 46, "(AT)23".to_string())]);

        // A 13 bp alternating run is (TA)[0,12) + (AT)[1,13): rotation dups,
        // only the leftmost stays. Cross-period overlap behaviour is covered
        // by test_dedup_cross_period_overlap_kept (motif mode).
        let mut seqs = HashMap::new();
        seqs.insert("alt".to_string(), "TATATATATATAT".to_string());
        let res: HashMap<String, crate::Region2Motif<Arc<String>>> =
            all_seq_hp_tr_finder(&all_regs, &seqs);
        let mut kept: Vec<_> = res["alt"]
            .value
            .iter()
            .map(|((start, end), _)| (*start, *end))
            .collect();
        kept.sort();
        assert_eq!(kept, vec![(0, 12)]);
    }

    /// Feeds `dedup_rotation_hits` a set of hand-built annotations and returns
    /// what survived, sorted by start, so a test can spell out exactly which
    /// regions are kept. Patterns are `(motif)N` as the real finder emits them,
    /// with the copy count matching the span; the period is derived from the
    /// motif length, mirroring the pipeline.
    fn dedup_keep(entries: &[((usize, usize), &str)]) -> Vec<(usize, usize, String)> {
        let mut input: Region2Motif<Arc<String>> = Region2Motif::default();
        for ((start, end), pattern) in entries {
            input.insert((*start, *end), Arc::new(pattern.to_string()));
        }

        let mut kept: Vec<_> = dedup_rotation_hits(input)
            .value
            .into_iter()
            .map(|((start, end), pattern)| (start, end, pattern.to_string()))
            .collect();
        kept.sort();
        kept
    }

    #[test]
    fn test_dedup_2mer_rotation() {
        // 复现 STR27-1 场景：(AT)×5 与 (TA)×5 是同一 TR 的旋转，保留最左的 (AT)
        let kept = dedup_keep(&[((6, 16), "(AT)5"), ((7, 17), "(TA)5")]);
        assert_eq!(kept, vec![(6, 16, "(AT)5".to_string())]);
    }

    #[test]
    fn test_dedup_3mer_rotation() {
        // 周期 3 的三种旋转互为冗余，只保留起点最靠左的一条
        let kept = dedup_keep(&[
            ((0, 12), "(ATG)4"),
            ((1, 13), "(TGA)4"),
            ((2, 14), "(GAT)4"),
        ]);
        assert_eq!(kept, vec![(0, 12, "(ATG)4".to_string())]);
    }

    #[test]
    fn test_dedup_4mer_rotation() {
        // 周期 4 的四种旋转只保留一条
        let kept = dedup_keep(&[
            ((0, 12), "(ATGC)3"),
            ((1, 13), "(TGCA)3"),
            ((2, 14), "(GCAT)3"),
            ((3, 15), "(CATG)3"),
        ]);
        assert_eq!(kept, vec![(0, 12, "(ATGC)3".to_string())]);
    }

    #[test]
    fn test_dedup_same_period_non_overlap_kept() {
        // 同周期不重叠的两个不同 TR 都保留
        let kept = dedup_keep(&[((0, 6), "(AT)3"), ((20, 26), "(GT)3")]);
        assert_eq!(
            kept,
            vec![
                (0, 6, "(AT)3".to_string()),
                (20, 26, "(GT)3".to_string()),
            ]
        );
    }

    #[test]
    fn test_dedup_boundary_sharing_distinct_trs_kept() {
        // 边界共享型相邻同周期 TR：shift=5 >= p=3，即使重叠也必须都保留
        // AAT×2 0-6 与 TGG×2 5-11（TGG 不是 AAT 的旋转）
        let kept = dedup_keep(&[((0, 6), "(AAT)2"), ((5, 11), "(TGG)2")]);
        assert_eq!(
            kept,
            vec![
                (0, 6, "(AAT)2".to_string()),
                (5, 11, "(TGG)2".to_string()),
            ]
        );
    }

    #[test]
    fn test_dedup_cross_period_overlap_kept() {
        // 不同周期的重叠命中是真实的不同注释，保留
        let kept = dedup_keep(&[((0, 8), "(AT)4"), ((2, 8), "(ATG)2")]);
        assert_eq!(
            kept,
            vec![
                (0, 8, "(AT)4".to_string()),
                (2, 8, "(ATG)2".to_string()),
            ]
        );
    }

    #[test]
    fn test_dedup_same_start_keep_longer() {
        // 同起点时保留区间更长（信息量更大）的一条；输入乱序也应正确
        let kept = dedup_keep(&[((1, 11), "(AT)5"), ((1, 13), "(TA)6")]);
        assert_eq!(kept, vec![(1, 13, "(TA)6".to_string())]);
    }

    #[test]
    fn test_dedup_unordered_input_and_mixed_periods() {
        // 混合周期 + 乱序输入：两个独立 TR 各自的旋转冗余均被消除
        let kept = dedup_keep(&[
            ((7, 13), "(TA)3"),
            ((40, 49), "(TGG)3"),
            ((6, 12), "(AT)3"),
            ((41, 50), "(GGT)3"),
            ((60, 65), "(A)5"),
        ]);
        assert_eq!(
            kept,
            vec![
                (6, 12, "(AT)3".to_string()),
                (40, 49, "(TGG)3".to_string()),
                (60, 65, "(A)5".to_string()),
            ]
        );
    }

    #[test]
    fn test_shift_ge_period_kept_even_if_overlapping() {
        // shift >= p 的同周期命中是真实的不同 TR，即使重叠也保留；
        // 同起点 shift=0 触发防御性 tie-break，只留更长的一条。
        let kept = dedup_keep(&[
            ((0, 16), "(AT)8"),
            ((0, 6), "(AT)3"),
            ((5, 15), "(TA)5"),
            ((4, 12), "(ATAT)2"),
        ]);
        assert_eq!(
            kept,
            vec![
                (0, 16, "(AT)8".to_string()),
                (4, 12, "(ATAT)2".to_string()),
                (5, 15, "(TA)5".to_string()),
            ]
        );
    }

    #[test]
    fn test_touching_annotations_are_kept() {
        // 相邻（不重叠）的同周期命中 shift >= p，全部保留
        let kept = dedup_keep(&[((0, 6), "(A)6"), ((6, 16), "(C)10")]);
        assert_eq!(
            kept,
            vec![
                (0, 6, "(A)6".to_string()),
                (6, 16, "(C)10".to_string()),
            ]
        );

        // 一条首尾相接的同聚物链
        let kept = dedup_keep(&[((0, 10), "(A)10"), ((10, 25), "(C)15"), ((25, 35), "(G)10")]);
        assert_eq!(
            kept,
            vec![
                (0, 10, "(A)10".to_string()),
                (10, 25, "(C)15".to_string()),
                (25, 35, "(G)10".to_string()),
            ]
        );
    }

    #[test]
    fn test_expand_and_merge_regions() {
        // 扩区保持输入顺序，不排序
        let mut regions = vec![(10, 20), (0, 5)];
        expand_regions(&mut regions, 3, 3, 30);
        assert_eq!(regions, vec![(7, 23), (0, 8)]);

        // 重叠与相邻合并
        let mut regions = vec![(4, 9), (2, 3), (1, 4)];
        assert_eq!(merge_region(&mut regions), vec![(1, 9)]);
        let mut regions = vec![(0, 5), (5, 10), (20, 25)];
        assert_eq!(merge_region(&mut regions), vec![(0, 10), (20, 25)]);
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
