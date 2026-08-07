use clap::{Parser, ValueEnum};
use hp_tr_finder::{all_seq_hp_tr_finder, Region2Motif, UnitAndRepeats};
use regex::Regex;
use std::{
    collections::HashMap,
    error::Error,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
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
/// `UnitAndRepeats` (which enumerates non-degenerate motifs of that size).
/// motif mode: a single map built directly from the explicit `motif-N`
/// pairs, mirroring tr-finder's `build_tr_finder_regex_through_modif_min_repeats`.
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
            regs.push(UnitAndRepeats::new(unit, min_repeats).build_finder_regrex());
        }
        Ok(regs)
    } else if let Some(spec) = &cli.motif_and_min_repeats {
        let mut regs = HashMap::new();
        for item in spec.trim().split(',') {
            let item = item.trim();
            let (motif, min_repeats) = item
                .split_once('-')
                .ok_or_else(|| format!("invalid motif-repeats spec: '{item}', expected N-N"))?;
            let min_repeats: usize = min_repeats
                .parse()
                .map_err(|_| format!("invalid min repeats: '{min_repeats}'"))?;
            let regex_str = format!("({motif}){{{min_repeats},}}");
            let reg = Regex::new(&regex_str)
                .map_err(|e| format!("invalid regex from motif '{motif}': {e}"))?;
            regs.insert(motif.to_string(), reg);
        }
        Ok(vec![regs])
    } else {
        Err("must specify --unitAndRepeats or --motifAndRepeats".to_string())
    }
}

/// Parse a FASTA file into (id, seq) records in file order. Duplicate ids are
/// merged by concatenation (scaffold-style FASTA can repeat a contig name);
/// the FASTA header's first whitespace-separated token is the id.
fn parse_fasta(path: &str) -> Result<Vec<(String, String)>, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut records: Vec<(String, String)> = Vec::new();
    let mut current: Option<(String, String)> = None;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if let Some(record) = current.take() {
                records.push(record);
            }
            let id = line[1..].split_whitespace().next().unwrap_or("").to_string();
            current = Some((id, String::new()));
        } else if let Some((_, seq)) = current.as_mut() {
            seq.push_str(line.trim());
        }
    }
    if let Some(record) = current.take() {
        records.push(record);
    }

    Ok(records)
}

/// Merge duplicate ids, concatenating their sequences.
fn merge_records(records: Vec<(String, String)>) -> HashMap<String, String> {
    let mut seqs: HashMap<String, String> = HashMap::new();
    for (id, seq) in records {
        if let Some(prev) = seqs.get_mut(&id) {
            prev.push_str(&seq);
        } else {
            seqs.insert(id, seq);
        }
    }
    seqs
}

fn tr_finder(cli: &Cli) -> Result<(), Box<dyn Error>> {
    if cli.scenario != Scenario::Ref {
        return Err("called scenario not implemented yet".into());
    }

    let seqs = merge_records(parse_fasta(&cli.inp)?);
    let regs = build_finder_regs(cli).map_err(|e| -> Box<dyn Error> { e.into() })?;

    let res: HashMap<String, Region2Motif<Arc<String>>> = all_seq_hp_tr_finder(&regs, &seqs);

    let cmd_line = std::env::args().collect::<Vec<_>>().join(" ");
    let mut writer = BufWriter::new(File::create(cli.get_oup_filepath())?);
    writeln!(&mut writer, "##gff-version 3")?;
    writeln!(&mut writer, "# {cmd_line}")?;

    // Iterate in FASTA order, not HashMap order.
    for (name, seq) in &seqs {
        let mut annotations: Vec<_> = res.get(name).into_iter().flat_map(|r| r.iter()).collect();
        annotations.sort_by_key(|&(&(start, _), _)| start);

        for (&(start, end), pattern) in annotations {
            // pattern is `(motif)N`; DNA motifs never contain parentheses.
            let motif = &pattern[pattern.find('(').unwrap() + 1..pattern.rfind(')').unwrap()];
            let copies = (end - start) / motif.len();
            writeln!(
                &mut writer,
                "{name}\t.\t.\t{}\t{end}\t.\t.\t.\tNote=({motif}){copies},{}",
                start + 1, // GFF3 start is 1-based; end stays 0-based like tr-finder
                &seq[start..end]
            )?;
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
    fn scenario_value_enum_parses() {
        assert!(matches!(cli_from(&["ref", "t.fa"]).scenario, Scenario::Ref));
        assert!(matches!(
            cli_from(&["called", "t.fa"]).scenario,
            Scenario::Called
        ));
    }
}
