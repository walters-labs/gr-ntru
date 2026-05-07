use std::env;
use std::time::Instant;

use ntru_group_algebra::{FiniteGroup, GroupAlgebraNtru, NtruError, TrialSummary};

#[derive(Clone, Debug)]
struct Experiment {
    name: &'static str,
    group: &'static str,
    n: usize,
    p: u64,
    q: u64,
    d: usize,
    trials: usize,
}

const QUICK: &[Experiment] = &[
    Experiment {
        name: "classical cyclic C7",
        group: "cyclic",
        n: 7,
        p: 3,
        q: 41,
        d: 2,
        trials: 8,
    },
    Experiment {
        name: "symmetric S3",
        group: "symmetric",
        n: 3,
        p: 5,
        q: 67,
        d: 2,
        trials: 8,
    },
    Experiment {
        name: "dihedral D8",
        group: "dihedral",
        n: 8,
        p: 3,
        q: 97,
        d: 2,
        trials: 4,
    },
];

const FFT: &[Experiment] = &[
    Experiment {
        name: "dihedral D64 over FFT prime",
        group: "dihedral",
        n: 64,
        p: 3,
        q: 2_013_265_921,
        d: 16,
        trials: 2,
    },
    Experiment {
        name: "symmetric S5 over prime fields",
        group: "symmetric",
        n: 5,
        p: 7,
        q: 4099,
        d: 15,
        trials: 2,
    },
];

fn main() -> Result<(), NtruError> {
    let mut profile = String::from("quick");
    let mut seed = 20260506_u64;
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--profile" => {
                profile = args.next().unwrap_or_else(|| "quick".to_string());
            }
            "--seed" => {
                if let Some(value) = args.next() {
                    seed = value.parse().unwrap_or(seed);
                }
            }
            "--help" | "-h" => {
                print_help();
                return Ok(());
            }
            _ => {}
        }
    }

    let experiments = match profile.as_str() {
        "quick" => QUICK,
        "fft" => FFT,
        _ => QUICK,
    };

    let mut total_successes = 0;
    let mut total_trials = 0;
    for (index, experiment) in experiments.iter().enumerate() {
        let summary = run_experiment(experiment, seed + index as u64)?;
        total_successes += summary.successes;
        total_trials += summary.completed_trials;
    }
    println!("\nTotal successful decryptions: {total_successes}/{total_trials}");
    Ok(())
}

fn print_help() {
    println!("ntru-group-algebra --profile quick|fft [--seed N]");
}

fn run_experiment(experiment: &Experiment, seed: u64) -> Result<TrialSummary, NtruError> {
    let group = match experiment.group {
        "cyclic" => FiniteGroup::cyclic(experiment.n)?,
        "dihedral" => FiniteGroup::dihedral(experiment.n)?,
        "symmetric" => FiniteGroup::symmetric(experiment.n)?,
        _ => unreachable!(),
    };
    let mut scheme = GroupAlgebraNtru::new(group, experiment.p, experiment.q, experiment.d)?;
    let start = Instant::now();
    let summary = scheme.run_trials(experiment.trials, 1_000, seed);
    print_summary(experiment, &scheme, &summary, start.elapsed().as_secs_f64());
    Ok(summary)
}

fn print_summary(
    experiment: &Experiment,
    scheme: &GroupAlgebraNtru,
    summary: &TrialSummary,
    elapsed: f64,
) {
    println!("\n{}", experiment.name);
    println!(
        "  |G|={}, p={}, q={}, d={}",
        scheme.group().order(),
        scheme.p(),
        scheme.q(),
        scheme.d()
    );
    println!(
        "  decryptions: {}/{} successful ({} keygen failures)",
        summary.successes, summary.completed_trials, summary.keygen_failures
    );
    println!(
        "  no-wrap checks: {}/{} matched p*g*r + f*m",
        summary.no_wraps, summary.completed_trials
    );
    println!("  elapsed: {elapsed:.2}s");
    if let Some(avg) = summary.avg_key_attempts {
        println!(
            "  key search attempts: avg {:.2}, max {}",
            avg,
            summary.max_key_attempts.unwrap_or(0)
        );
    }
    if !summary.backend_stats.counts().is_empty() {
        println!("  backend calls: {}", summary.backend_stats.format());
    }
}
