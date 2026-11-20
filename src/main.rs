use clap::{Parser, ValueEnum};
use hp_tr_finder::{
    all_seq_hp_tr_finder, expand_regions, merge_region, simplify_motif, Region2Motif,
    UnitAndRepeats,
};
use noodles_fasta as fasta;
use regex::Regex;
use std::{
    collections::HashMap,
    error::Error,
    fs::File,
    io::{BufReader, BufWriter, Write},
    sync::Arc,
};

#[derive(Copy, Clone, PartialEq, Eq, ValueEnum)]
enum Scenario {
    Ref,
    Called,
}

#[derive(Parser)]
#[command(
    name = "hp_tr_finder",
    version,
    about = "Find homopolymers and tandem repeats in a FASTA file.",
    long_about = "\
Given the fasta file and units_and_min_repeats setting, the program will output the tandem-repeats area in the fasta file.

Algorithm overview: for each seq in the fasta file, do the following steps
    1) generate regex for interested region. for example 1-2,2-2 will generate [ACGT]{2,} (AC){2,}, (AG){2,}, (AT){2,}, ... regexes.
    2) for each seq in the fasta file, extract the regions that contains the regex pattern, then output the result.
"
)]
struct Cli {
    /// which scenario the program runs on
    #[arg(help = "which scenario the program runs on")]
    scenario: Scenario,

    /// input fasta file
    #[arg(help = "input fasta file")]
    inp: String,

    /// unit and min repeats. 1-2,2-2,3-2,4-2
    #[arg(long = "unitAndRepeats", help = "unit and min repeats. 1-2,2-2,3-2,4-2", group = "spec")]
    unit_and_min_repeats: Option<String>,

    /// motif and min repeats. A-2,AT-2
    #[arg(long = "motifAndRepeats", help = "motif and min repeats. A-2,AT-2", group = "spec")]
    motif_and_min_repeats: Option<String>,

    /// output filepath; if not given, ${inp}.gff will be generated
    #[arg(short = 'o', long = "oupPath", help = "output filepath; if not given, ${inp}.gff will be generated")]
    oup_filepath: Option<String>,

    /// extending the tr-region. Only effective in 'called' scenario
    #[arg(long = "extLeft", default_value_t = 0, help = "extending the tr-region. Only effective in 'called' scenario")]
    extend_left: usize,

    /// extending the tr-region. Only effective in 'called' scenario
    #[arg(long = "extRight", default_value_t = 0, help = "extending the tr-region. Only effective in 'called' scenario")]
    extend_right: usize,

    /// extend left/right according to origin_tr_region_len * extRatio; if
    /// provided, it will invalidate extLeft/extRight. Only effective in 'called' scenario
    #[arg(long = "extRatio", help = "extend left right according to origin_tr_region_len * extRatio. if provided, it will invalidate extLeft/extRight. Only effective in 'called' scenario")]
    extend_ratio: Option<f32>,

    /// extMax and extRatio are used together. Only effective in 'called' scenario
    #[arg(long = "extMax", help = "extMax and extRatio are used together. Only effective in 'called' scenario")]
    extend_max: Option<usize>,
}

impl Cli {
    fn get_oup_filepath(&self) -> String {
        match self.oup_filepath.as_ref() {
            Some(path) => path.clone(),
            None => format!(
                "{}.gff",
                self.inp.rsplit_once('.').map_or(self.inp.as_str(), |(s, _)| s)
            ),
        }
    }
}

/// Build the regex maps used for detection.
///
/// unit mode: one `HashMap<motif, Regex>` per unit size, via the lib's
/// `UnitAndRepeats` (which enumerates primitive/non-degenerate motifs of that
/// size, so `AAAA`/`GTGT` no longer yield duplicate regexes).
/// motif mode: a single map built from the explicit `motif-N` pairs, with each
/// motif simplified to its primitive form and `min_repeats` scaled by the copy
/// count (`GTGT-3` == `GT-6`), mirroring tr-finder's
/// `build_tr_finder_regex_through_modif_min_repeats`.
fn build_finder_regs(cli: &Cli) -> Result<Vec<HashMap<String, Regex>>, String> {
    if let Some(spec) = &cli.unit_and_min_repeats {
        let mut regs = Vec::new();
        for item in spec.trim().split(',') {
            let (unit, min_repeats) = item
                .split_once('-')
                .ok_or_else(|| format!("invalid unit-repeats spec: '{item}', expected N-N"))?;
            let unit: u8 = unit
                .parse()
                .map_err(|_| format!("invalid unit size: '{unit}'"))?;
            let min_repeats: u8 = min_repeats
                .parse()
                .map_err(|_| format!("invalid min repeats: '{min_repeats}'"))?;
            if unit == 0 {
                return Err(format!("unit size must be >= 1, got '{unit}'"));
            }
            if min_repeats == 0 {
                // `(motif){0,}` 会匹配空串，产出无意义的零长度注释
                return Err(format!("min repeats must be >= 1, got '{min_repeats}'"));
            }
            regs.push(UnitAndRepeats::new(unit, min_repeats).build_finder_regrex());
        }
        Ok(regs)
    } else if let Some(spec) = &cli.motif_and_min_repeats {
        // 将每个 motif 化简为最简（primitive）形式，min_repeats 按拷贝数等比放大：
        // 例如 GTGT-3 等价于 GT-6（同一段 ≥12bp 的 GT 串联重复）。
        // 化简后出现重复 motif（如 GT-2 与 GTGT-3 都化简为 GT）时，保留更严格的 min_repeats。
        let mut simplified: HashMap<String, usize> = HashMap::new();
        for item in spec.trim().split(',') {
            let item = item.trim();
            let (motif, min_repeats) = item
                .split_once('-')
                .ok_or_else(|| format!("invalid motif-repeats spec: '{item}', expected N-N"))?;
            let min_repeats: usize = min_repeats
                .parse()
                .map_err(|_| format!("invalid min repeats: '{min_repeats}'"))?;
            if motif.is_empty() {
                // 例如 `--motifAndRepeats=-2`：split_once 把 '-' 当分隔符，motif 为空
                return Err(format!("motif must not be empty, got spec '{item}'"));
            }
            if min_repeats == 0 {
                // `(motif){0,}` 会匹配空串，产出无意义的零长度注释
                return Err(format!("min repeats must be >= 1, got '{min_repeats}'"));
            }
            let (root, copy_num) = simplify_motif(motif);
            let entry = simplified.entry(root).or_insert(0);
            *entry = (*entry).max(min_repeats * copy_num);
        }

        let mut regs = HashMap::new();
        for (motif, min_repeats) in simplified {
            let regex_str = format!("({motif}){{{min_repeats},}}");
            let reg = Regex::new(&regex_str)
                .map_err(|e| format!("invalid regex from motif '{motif}': {e}"))?;
            regs.insert(motif, reg);
        }
        Ok(vec![regs])
    } else {
        Err("must specify --unitAndRepeats or --motifAndRepeats".to_string())
    }
}

/// Parse a FASTA file into (id, seq) records in file order via noodles-fasta.
/// Duplicate ids are merged by concatenation at parse time (scaffold-style
/// FASTA can repeat a contig name), first-appearance order is preserved; the
/// FASTA header's first whitespace-separated token is the id.
fn parse_fasta(path: &str) -> Result<Vec<(String, String)>, Box<dyn Error>> {
    let file = File::open(path)?;
    let mut reader = fasta::io::Reader::new(BufReader::new(file));

    let mut records: Vec<(String, String)> = Vec::new();
    let mut index: HashMap<String, usize> = HashMap::new();

    for result in reader.records() {
        let record = result?;
        let id = str::from_utf8(record.name())
            .map_err(|e| format!("FASTA record name in {path} is not valid UTF-8: {e}"))?
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        let seq = str::from_utf8(record.sequence().as_ref())
            .map_err(|e| format!("FASTA sequence for {id} is not valid UTF-8: {e}"))?
            .to_string();

        match index.get(&id) {
            Some(&i) => records[i].1.push_str(&seq),
            None => {
                index.insert(id.clone(), records.len());
                records.push((id, seq));
            }
        }
    }

    Ok(records)
}

/// 计算一个命中区间的左右扩区量：提供 extRatio 时按区域长度比例扩区
/// （左右各 `region_len * ratio`，受 extMax 上限约束），且优先于 extLeft/extRight；
/// 否则使用 extLeft/extRight。
fn extend_amounts(cli: &Cli, region_len: usize) -> (usize, usize) {
    match cli.extend_ratio {
        Some(ratio) => {
            let ext = (region_len as f32 * ratio) as usize;
            let ext = cli.extend_max.map_or(ext, |max| ext.min(max));
            (ext, ext)
        }
        None => (cli.extend_left, cli.extend_right),
    }
}

fn tr_finder(cli: &Cli) -> Result<(), Box<dyn Error>> {
    // 输出按 FASTA 文件顺序（重复 id 已在解析时合并）。
    let records = parse_fasta(&cli.inp)?;
    let seqs: HashMap<String, String> = records.iter().cloned().collect();
    let regs = build_finder_regs(cli).map_err(|e| -> Box<dyn Error> { e.into() })?;

    let res: HashMap<String, Region2Motif<Arc<String>>> = all_seq_hp_tr_finder(&regs, &seqs);

    let cmd_line = std::env::args().collect::<Vec<_>>().join(" ");
    let mut writer = BufWriter::new(File::create(cli.get_oup_filepath())?);
    writeln!(&mut writer, "##gff-version 3")?;
    writeln!(&mut writer, "# {cmd_line}")?;

    // Iterate in FASTA order, not HashMap order.
    for (name, seq) in &records {
        match cli.scenario {
            Scenario::Ref => {
                let mut annotations: Vec<_> =
                    res.get(name).into_iter().flat_map(|r| r.iter()).collect();
                // (start, end) 全序，避免同起点时输出顺序依赖 HashMap 迭代顺序
                annotations.sort_unstable_by_key(|&(&(start, end), _)| (start, end));

                for (&(start, end), pattern) in annotations {
                    // pattern is `(motif)N`; DNA motifs never contain parentheses.
                    let motif =
                        &pattern[pattern.find('(').unwrap() + 1..pattern.rfind(')').unwrap()];
                    let copies = (end - start) / motif.len();
                    writeln!(
                        &mut writer,
                        "{name}\t.\t.\t{}\t{end}\t.\t.\t.\tNote=({motif}){copies},{}",
                        start + 1, // GFF3 start is 1-based; end stays 0-based like tr-finder
                        &seq[start..end]
                    )?;
                }
            }
            Scenario::Called => {
                // 命中区间按参数扩区（边界 clamp 到 [0, seq.len()]），重叠合并后输出
                let mut regions: Vec<(usize, usize)> = res
                    .get(name)
                    .into_iter()
                    .flat_map(|r| r.keys().copied())
                    .collect();
                if cli.extend_ratio.is_none() {
                    expand_regions(&mut regions, cli.extend_left, cli.extend_right, seq.len());
                } else {
                    for region in &mut regions {
                        let (left, right) = extend_amounts(cli, region.1 - region.0);
                        region.0 = region.0.saturating_sub(left).min(seq.len());
                        region.1 = (region.1 + right).min(seq.len());
                    }
                }
                let regions = merge_region(&mut regions);

                for (start, end) in regions {
                    writeln!(
                        &mut writer,
                        "{name}\t.\t.\t{}\t{end}\t.\t.\t.\t{}",
                        start + 1, // GFF3 start is 1-based; end stays 0-based like tr-finder
                        &seq[start..end]
                    )?;
                }
            }
        }
    }

    writer.flush()?;
    Ok(())
}

fn main() {
    let cli = Cli::parse();
    if let Err(e) = tr_finder(&cli) {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::{build_finder_regs, Cli, Scenario};
    use clap::Parser;

    /// Construct a Cli from an argv slice (program name + args).
    fn cli_from(args: &[&str]) -> Cli {
        Cli::parse_from(std::iter::once("hp_tr_finder").chain(args.iter().copied()))
    }

    #[test]
    fn unit_mode_builds_one_regex_map_per_unit() {
        let cli = cli_from(&[
            "ref", "t.fa", "--unitAndRepeats", "1-3,2-2",
        ]);
        let regs = build_finder_regs(&cli).unwrap();
        assert_eq!(regs.len(), 2);

        // unit 1: homopolymer motifs only, min repeats 3.
        let unit1 = &regs[0];
        assert_eq!(
            unit1.keys().cloned().collect::<std::collections::BTreeSet<_>>(),
            ["A", "C", "G", "T"].into_iter().map(str::to_string).collect()
        );
        let a3 = &unit1["A"];
        assert_eq!(a3.as_str(), "(A){3,}");
        assert!(a3.is_match("AAAA"));
        assert!(!a3.is_match("AA"));

        // unit 2: all non-degenerate dinucleotide motifs, min repeats 2.
        let unit2 = &regs[1];
        assert_eq!(unit2.len(), 12);
        let ac2 = &unit2["AC"];
        assert_eq!(ac2.as_str(), "(AC){2,}");
        assert!(ac2.is_match("ACAC"));
        assert!(!ac2.is_match("AC"));
    }

    #[test]
    fn motif_mode_builds_a_single_map_from_explicit_motifs() {
        let cli = cli_from(&[
            "ref", "t.fa", "--motifAndRepeats", "A-2,AT-2",
        ]);
        let regs = build_finder_regs(&cli).unwrap();
        assert_eq!(regs.len(), 1);

        let map = &regs[0];
        assert_eq!(map.len(), 2);
        let at2 = &map["AT"];
        assert_eq!(at2.as_str(), "(AT){2,}");
        assert!(at2.is_match("ATATAT"));
        assert!(!at2.is_match("AT"));
    }

    #[test]
    fn motif_mode_simplifies_non_primitive_motifs() {
        // GTGT-3 化简为 GT-6
        let cli = cli_from(&["ref", "t.fa", "--motifAndRepeats", "GTGT-3"]);
        let regs = build_finder_regs(&cli).unwrap();
        assert!(regs[0].contains_key("GT") && !regs[0].contains_key("GTGT"));
        let gt = &regs[0]["GT"];
        assert!(gt.is_match("GTGTGTGTGTGT")); // 6 拷贝，命中
        assert!(!gt.is_match("GTGTGTGTGT")); // 5 拷贝，不命中

        // ATATAT-2 化简为 AT-6
        let cli = cli_from(&["ref", "t.fa", "--motifAndRepeats", "ATATAT-2"]);
        let regs = build_finder_regs(&cli).unwrap();
        assert!(regs[0].contains_key("AT") && !regs[0].contains_key("ATATAT"));

        // 重复化简到同一 motif 时保留更严格的 min_repeats：GT-2 + GTGT-3 -> GT-6
        let cli = cli_from(&["ref", "t.fa", "--motifAndRepeats", "GT-2,GTGT-3"]);
        let regs = build_finder_regs(&cli).unwrap();
        assert_eq!(regs[0].len(), 1);
        assert!(!regs[0]["GT"].is_match("GTGTGTGTGT")); // 5 拷贝，不满足 GT-6
    }

    #[test]
    fn extend_amounts_semantics() {
        // 默认使用 extLeft/extRight
        let cli = cli_from(&[
            "called",
            "t.fa",
            "--unitAndRepeats",
            "2-2",
            "--extLeft",
            "3",
            "--extRight",
            "5",
        ]);
        assert_eq!(super::extend_amounts(&cli, 10), (3, 5));

        // extRatio 优先于 extLeft/extRight：10bp 区域 * 0.5 → 左右各 5
        let cli = cli_from(&[
            "called",
            "t.fa",
            "--unitAndRepeats",
            "2-2",
            "--extLeft",
            "3",
            "--extRatio",
            "0.5",
        ]);
        assert_eq!(super::extend_amounts(&cli, 10), (5, 5));

        // extMax 封顶
        let cli = cli_from(&[
            "called",
            "t.fa",
            "--unitAndRepeats",
            "2-2",
            "--extRatio",
            "0.5",
            "--extMax",
            "3",
        ]);
        assert_eq!(super::extend_amounts(&cli, 10), (3, 3));
        assert_eq!(super::extend_amounts(&cli, 4), (2, 2));
    }

    #[test]
    fn missing_spec_is_an_error() {
        let cli = cli_from(&["ref", "t.fa"]);
        let err = build_finder_regs(&cli).unwrap_err();
        assert!(err.contains("must specify"), "unexpected error: {err}");
    }

    #[test]
    fn both_specs_conflict_at_clap_level() {
        let args = [
            "hp_tr_finder", "ref", "t.fa", "--unitAndRepeats", "1-2", "--motifAndRepeats", "A-2",
        ];
        let Err(err) = Cli::try_parse_from(args) else {
            panic!("expected a clap conflict error");
        };
        // clap group="spec": the two options cannot be used together.
        assert!(err.to_string().contains("cannot be used with"));
    }

    #[test]
    fn invalid_unit_spec_items_are_rejected() {
        // malformed item: no '-'
        let cli = cli_from(&["ref", "t.fa", "--unitAndRepeats", "1-2,zz"]);
        let err = build_finder_regs(&cli).unwrap_err();
        assert!(err.contains("expected N-N"), "unexpected error: {err}");

        // non-numeric min repeats
        let cli = cli_from(&["ref", "t.fa", "--unitAndRepeats", "1-x"]);
        let err = build_finder_regs(&cli).unwrap_err();
        assert!(err.contains("invalid min repeats"), "unexpected error: {err}");

        // zero unit size
        let cli = cli_from(&["ref", "t.fa", "--unitAndRepeats", "0-2"]);
        let err = build_finder_regs(&cli).unwrap_err();
        assert!(err.contains("unit size must be >= 1"), "unexpected error: {err}");
    }

    #[test]
    fn invalid_motif_spec_items_are_rejected() {
        // malformed item: no '-'
        let cli = cli_from(&["ref", "t.fa", "--motifAndRepeats", "A2"]);
        let err = build_finder_regs(&cli).unwrap_err();
        assert!(err.contains("invalid motif-repeats spec"), "unexpected error: {err}");

        // motif that yields an invalid regex: '[' alone is an unclosed class.
        let cli = cli_from(&["ref", "t.fa", "--motifAndRepeats", "[-2"]);
        let err = build_finder_regs(&cli).unwrap_err();
        assert!(err.contains("invalid regex"), "unexpected error: {err}");
    }

    #[test]
    fn empty_motif_and_zero_min_repeats_are_rejected() {
        // `--motifAndRepeats=-2`：motif 为空（'-' 被当作分隔符），不得 panic
        let cli = cli_from(&["ref", "t.fa", "--motifAndRepeats=-2"]);
        let err = build_finder_regs(&cli).unwrap_err();
        assert!(err.contains("motif must not be empty"), "unexpected error: {err}");

        // min repeats 为 0 会产生 {0,} 的零长度注释，两种模式都应拒绝
        let cli = cli_from(&["ref", "t.fa", "--unitAndRepeats", "1-0"]);
        let err = build_finder_regs(&cli).unwrap_err();
        assert!(err.contains("min repeats must be >= 1"), "unexpected error: {err}");

        let cli = cli_from(&["ref", "t.fa", "--motifAndRepeats", "AT-0"]);
        let err = build_finder_regs(&cli).unwrap_err();
        assert!(err.contains("min repeats must be >= 1"), "unexpected error: {err}");
    }

    /// 写一个临时 FASTA，跑完整 CLI 流程，返回去除 `#` 头后的 GFF 行。
    fn run_on_seq(tmp_tag: &str, seq: &str, scenario: &str, args: &[&str]) -> Vec<String> {
        let dir = std::env::temp_dir();
        let stem = format!("hp_tr_finder_test_{}_{tmp_tag}", std::process::id());
        let fa = dir.join(format!("{stem}.fa"));
        let out = dir.join(format!("{stem}.gff"));
        std::fs::write(&fa, format!(">t1\n{seq}\n")).unwrap();

        let mut argv = vec![scenario, fa.to_str().unwrap(), "-o", out.to_str().unwrap()];
        argv.extend_from_slice(args);
        let cli = cli_from(&argv);
        super::tr_finder(&cli).unwrap();

        let content = std::fs::read_to_string(&out).unwrap();
        let _ = std::fs::remove_file(&fa);
        let _ = std::fs::remove_file(&out);
        content
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(str::to_string)
            .collect()
    }

    #[test]
    fn called_scenario_extends_and_merges_regions() {
        // (ATG)6 at [10,28), (TGCA)5 at [34,54)
        let seq = "ACGTACGTAAATGATGATGATGATGATGNNNNNNTGCATGCATGCATGCATGCA";

        // extLeft/extRight 3：两端各扩 3bp，[7,31) 与 [31,57) 相邻合并为 [7,57)
        let lines = run_on_seq(
            "called_ext",
            seq,
            "called",
            &["--unitAndRepeats", "2-5,3-4,4-4", "--extLeft", "3", "--extRight", "3"],
        );
        assert_eq!(lines.len(), 1, "unexpected output: {lines:?}");
        let cols: Vec<_> = lines[0].split('\t').collect();
        // 右端扩区 clamp 到序列末尾（本例 seq.len() == 54）
        assert_eq!((cols[3], cols[4]), ("8", "54"), "unexpected line: {}", lines[0]);
        assert_eq!(cols[8], &seq[7..54]);

        // extRatio 0.5 + extMax 4：18bp/20bp 区域各扩 9/10bp，被 extMax 封顶到 4
        let lines = run_on_seq(
            "called_ratio",
            seq,
            "called",
            &[
                "--unitAndRepeats",
                "2-5,3-4,4-4",
                "--extRatio",
                "0.5",
                "--extMax",
                "4",
            ],
        );
        assert_eq!(lines.len(), 1, "unexpected output: {lines:?}");
        let cols: Vec<_> = lines[0].split('\t').collect();
        assert_eq!((cols[3], cols[4]), ("7", "54"), "unexpected line: {}", lines[0]);
        assert_eq!(cols[8], &seq[6..54]);
    }

    #[test]
    fn ref_scenario_dedups_rotations_and_keeps_simplified_motif() {
        // (ATG)6 的旋转 (TGA)/(GAT) 应被去除，只留一条
        let lines = run_on_seq("ref_dedup", "ACGTACGTAAATGATGATGATGATGATG", "ref", &[
            "--unitAndRepeats",
            "2-5,3-4,4-4",
        ]);
        assert_eq!(lines.len(), 1, "unexpected output: {lines:?}");
        assert!(lines[0].ends_with("Note=(ATG)6,ATGATGATGATGATGATG"), "{}", lines[0]);

        // 显式非最简 motif 化简后按最简形式报告
        let lines = run_on_seq("ref_simplify", "NNNNNNNNNNGTGTGTGTGTGT", "ref", &[
            "--motifAndRepeats",
            "GTGT-3",
        ]);
        assert_eq!(lines.len(), 1, "unexpected output: {lines:?}");
        assert!(lines[0].ends_with("Note=(GT)6,GTGTGTGTGTGT"), "{}", lines[0]);
    }

    #[test]
    fn scenario_value_enum_parses() {
        assert!(matches!(cli_from(&["ref", "t.fa"]).scenario, Scenario::Ref));
        assert!(matches!(
            cli_from(&["called", "t.fa"]).scenario,
            Scenario::Called
        ));
    }
}
