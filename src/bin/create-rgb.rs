use anyhow::Result;
use log::*;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;
use simplelog::*;
use std::path::PathBuf;
use structopt::StructOpt;

use rgb::forward;
use rgb::ciff;
use rgb::output;

#[derive(Debug, StructOpt)]
#[structopt(
    name = "create-rgb",
    about = "Reorders documents using recursive graph bisection and ciff files."
)]
struct Opt {
    /// Input file ciff file
    #[structopt(short, long, parse(from_os_str))]
    input: PathBuf,

    /// Output ciff file
    #[structopt(short, long, parse(from_os_str))]
    output: PathBuf,

    /// Minimum number of occurrences to consider
    #[structopt(short, long, default_value = "4096")]
    min_len: usize,

    /// Maximum length to consider in percentage of documents in the index
    #[structopt(short, long, default_value = "0.1")]
    cutoff_frequency: f32,

    /// Min partition size
    #[structopt(short, long, default_value = "16")]
    recursion_stop: usize,

    /// Swap iterations
    #[structopt(short, long, default_value = "20")]
    swap_iterations: usize,

    /// Show loggap cost
    #[structopt(short, long)]
    loggap: bool,

    /// Sort leaf by identifier
    #[structopt(short, long)]
    sort_leaf: bool,
    
}

fn compute_loggapsum<P: AsRef<std::path::Path>>(file_path: P) -> (f64, usize) {
    let ciff_file = std::fs::File::open(file_path).expect("can't open ciff file");
    let mut ciff_file = std::io::BufReader::new(ciff_file);
    let mut ciff_reader = ciff::Reader::new(&mut ciff_file).expect("can't open ciff file");

    let log_gaps: Vec<f64> = (0..256).map(|i| (i as f64).log2()).collect();

    let mut log_sum = 0.0;
    let mut num_postings: usize = 0;
    while let Some(ciff::CiffRecord::PostingsList(plist)) = ciff_reader.next() {
        let postings = plist.get_postings();
        log_sum += ((postings.first().unwrap().get_docid() + 1) as f64).log2();
        num_postings += 1;
        for posting in postings.iter().skip(1) {
            let did_gap = posting.get_docid();
            if did_gap < 256 {
                log_sum += unsafe { log_gaps.get_unchecked(did_gap as usize) };
            } else {
                log_sum += (did_gap as f64).log2();
            }
            num_postings += 1;
        }
    }
    (log_sum, num_postings)
}

fn main() -> Result<()> {
    CombinedLogger::init(vec![
        TermLogger::new(LevelFilter::Info, Config::default(), TerminalMode::Mixed),
        WriteLogger::new(
            LevelFilter::Info,
            Config::default(),
            std::fs::File::create(format!(
                "create-rgb-{}.log",
                std::process::id()
            ))
            .expect("can't create file log"),
        ),
    ])
    .unwrap();

    let opt = Opt::from_args();
    info!("{:?}", opt);

    let start_fwd = std::time::Instant::now();
    info!("(1) create forward index from ciff");
    let forward::Forward { mut docs, uniq_terms } = forward::from_ciff(&opt.input, opt.min_len, opt.cutoff_frequency)?;

    // move all empty docs to the end and exclude them from reordering
    info!("(2) sort empty docs to the back");
    docs.sort_by(|a, b| b.terms.len().cmp(&a.terms.len()));
    let num_non_empty = docs
        .iter()
        .position(|d| d.terms.len() == 0)
        .unwrap_or(docs.len());
    let fwd_time = start_fwd.elapsed().as_secs_f32();
    info!("fwd duration: {:.2} secs", fwd_time);
    info!("docs {} non_empty {}", docs.len(), num_non_empty);

    let pb = indicatif::ProgressBar::new(num_non_empty as u64);
    pb.set_draw_delta(1000);
    let desc =
        "RGB: [{elapsed_precise}] [{bar:40.cyan/blue}] ({pos}/{len}, ETA {eta}, SPEED: {per_sec})";
    pb.set_style(indicatif::ProgressStyle::default_bar().template(desc));

    info!("(3) perform graph bisection");
    let start_rgb = std::time::Instant::now();
    rgb::recursive_graph_bisection(
        &mut docs[..num_non_empty],
        uniq_terms,
        opt.swap_iterations,
        opt.recursion_stop,
        pb.clone(),
        opt.sort_leaf,
    );
    pb.finish_and_clear();
    let rgb_time = start_rgb.elapsed().as_secs_f32();
    info!("rgb duration: {:.2} secs", rgb_time);

    // now we can clear some space
    info!("(4) clear forward index");
    docs.par_iter_mut().for_each(|doc| {
        doc.terms.truncate(0);
        doc.terms.shrink_to_fit();
    });

    info!("(5) write new ciff file");
    let start_write = std::time::Instant::now();
    output::rewrite_ciff(&docs, &opt.input, &opt.output)?;
    let write_time = start_write.elapsed().as_secs_f32();
    info!("write duration: {:.2} secs", write_time);


    let all_done_time = start_fwd.elapsed().as_secs_f32();

    if opt.loggap {
        info!("(6) compute loggap cost");
        let (before_log_sum, num_postings) = compute_loggapsum(&opt.input);
        let before_bpi = before_log_sum / num_postings as f64;
        info!("\tbefore reorder: {:.3} BPI",before_bpi);
        let (after_log_sum, num_postings) = compute_loggapsum(&opt.output);
        let after_bpi = after_log_sum / num_postings as f64;
        info!("\t after reorder: {:.3} BPI",after_bpi);
    }

    info!("ALL DONE! duration: {:.2} secs", all_done_time);

    Ok(())
}
