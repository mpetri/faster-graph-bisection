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
    output_ciff: Option<PathBuf>,

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
    #[structopt(long)]
    sort_leaf: bool,

    /// Maximum depth
    #[structopt(long, default_value = "100")]
    max_depth: usize,

    /// Read forward index
    #[structopt(long, parse(from_os_str))]
    input_fidx: Option<PathBuf>,

    /// Output forward index
    #[structopt(long, parse(from_os_str))]
    output_fidx: Option<PathBuf>,

    /// Dump the document map
    #[structopt(long, parse(from_os_str))]
    output_mapping: Option<PathBuf>,
    
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

fn validate_gain() {

    let gain_func: Option<&'static str> = std::option_env!("GAIN");
    if gain_func.is_none() {
        log::info!("Error: A gain function needs to be passed at compile time via the environment variable `GAIN` -- Please recompile...");
        std::process::exit(1);
    }
    let gain_func = gain_func.unwrap();
    let gain_types = vec!["default", "approx_1", "approx_2"];
    if gain_types.iter().any(|&i| i == gain_func) {
        log::info!("Using the `{}` gain function.", gain_func);
    } else {
        log::info!("Error: Couldn't match the gain function.");
        std::process::exit(1); 
    }
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

    // validate the compile-time gain function
    validate_gain();

    let opt = Opt::from_args();
    info!("{:?}", opt);

    // Sanity check output options. We want to at least dump the map, dump the ciff, dump the fidx, or compute loggap on the fly...
    if opt.output_fidx.is_none() && opt.output_ciff.is_none() && opt.output_mapping.is_none() && !opt.loggap{
        info!("Error: Nothing will be output. Check your options and try again.");
        std::process::exit(1);
    }

    // Check that we're not trying to both read and write a pre-existing forward index
    let start_fwd = std::time::Instant::now();
    let forward::Forward { mut docs, uniq_terms } = match opt.input_fidx {
    	Some(index) => {info!("(1) reading forward index from file"); forward::from_file(index)?},
    	None => {info!("(1) building forward index"); forward::from_ciff(&opt.input, opt.min_len, opt.cutoff_frequency, opt.output_fidx.as_ref())?},
    };

    // move all empty docs to the end and exclude them from reordering
    info!("(2) sort empty docs to the back");
    docs.sort_by(|a, b| b.postings.len().cmp(&a.postings.len()));
    let num_non_empty = docs
        .iter()
        .position(|d| d.postings.len() == 0)
        .unwrap_or(docs.len());
    let fwd_time = start_fwd.elapsed().as_secs_f32();
    info!("fwd duration: {:.2} secs", fwd_time);
    info!("docs {} non_empty {}", docs.len(), num_non_empty);
    info!("put docs back into default order...");
    docs.sort_by(|a, b| a.org_id.cmp(&b.org_id));

    info!("(3) perform graph bisection");
    let start_rgb = std::time::Instant::now();
    let depth = 1;
    rgb::recursive_graph_bisection(
        &mut docs[..num_non_empty],
        uniq_terms,
        opt.swap_iterations,
        opt.recursion_stop,
        opt.max_depth,
        depth,
        opt.sort_leaf,
    );
    let rgb_time = start_rgb.elapsed().as_secs_f32();
    info!("rgb duration: {:.2} secs", rgb_time);

    // now we can clear some space
    info!("(4) clear forward index");
    docs.par_iter_mut().for_each(|doc| {
        doc.postings.truncate(0);
        doc.postings.shrink_to_fit();
    });

    info!("(5) starting output operations...");

 
    if let Some(output_map) = opt.output_mapping {
        info!(" --> (5.1) output the plain-text mapping file");
        output::dump_order(&docs, output_map);
    }

    if let Some(output_ciff) = opt.output_ciff {
        info!(" --> (5.2) write new ciff file");
        let start_write = std::time::Instant::now();
        output::rewrite_ciff(&docs, &opt.input, &output_ciff)?;
        let write_time = start_write.elapsed().as_secs_f32();
        info!("write duration: {:.2} secs", write_time);
        
        if opt.loggap {
            info!("(6) compute loggap cost");
            let (before_log_sum, num_postings) = compute_loggapsum(&opt.input);
            let before_bpi = before_log_sum / num_postings as f64;
            info!("\tbefore reorder: {:.3} BPI",before_bpi);
            let (after_log_sum, num_postings) = compute_loggapsum(&output_ciff);
            let after_bpi = after_log_sum / num_postings as f64;
            info!("\t after reorder: {:.3} BPI",after_bpi);
        }
    } else {
        if opt.loggap {
            info!(" --> (5.2a) Computing LogGap on remapped CIFF file, but not writing it...");
            let start_write = std::time::Instant::now();
            output::remap_ciff(&docs, &opt.input)?;
            let write_time = start_write.elapsed().as_secs_f32();
            info!("Remap + LogGap duration: {:.2} secs", write_time);
        }
    }
    

    let all_done_time = start_fwd.elapsed().as_secs_f32();

    info!("ALL DONE! duration: {:.2} secs", all_done_time);
    
    Ok(())
}
