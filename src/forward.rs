use log::*;
use thiserror::Error;

use crate::ciff;

use std::fs::File;
use bincode::{serialize_into, deserialize_from};
use std::io::{BufWriter};
use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize, Default, PartialEq, Clone, Debug)]
pub struct Doc {
    pub terms: Vec<u32>,
    pub org_id: u32,
    pub gain: f32,
    pub leaf_id: i32,
}

#[derive(Error, Debug)]
pub enum Error {
    #[error("forward I/O error")]
    ReadError(#[from] std::io::Error),
    #[error("open ciff file")]
    CiffOpenError,
}

#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct Forward {
    pub docs: Vec<Doc>,
    pub uniq_terms: usize,
}


pub fn from_file<P: AsRef<std::path::Path>>(
    file_path: P
) -> Result<Forward, bincode::Error> {
    let in_file = std::io::BufReader::new(std::fs::File::open(file_path)?);
    deserialize_from(in_file)
}

pub fn from_ciff<P: AsRef<std::path::Path>>(
    file_path: P,
    min_len: usize,
    cutoff_frequency: f32,
    output_path: Option<P>,
) -> Result<Forward, Error> {
    let ciff_file = std::fs::File::open(file_path)?;
    let mut ciff_file = std::io::BufReader::new(ciff_file);
    let mut ciff_reader =
        ciff::Reader::new(&mut ciff_file).map_err(|_e| Error::CiffOpenError)?;
    let num_docs = ciff_reader.num_docs();
    let mut docs = Vec::with_capacity(num_docs as usize);
    for doc_id in 0..num_docs {
        docs.push(Doc {
            terms: Vec::with_capacity(256), // initial estimate for uniq terms in doc
            org_id: doc_id as u32,
            gain: 0.0,
            leaf_id: -1,
        });
    }
    let cutoff_len: usize = (num_docs as f32 * cutoff_frequency).ceil() as usize;

    let num_msgs = ciff_reader.num_postings_lists() as u64;
    let pb_plist = indicatif::ProgressBar::new(num_msgs);
    pb_plist.set_draw_delta(1000);
    pb_plist.set_style(indicatif::ProgressStyle::default_bar().template(
        "create_fwd: {spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] ({pos}/{len}, ETA {eta}, SPEED: {per_sec})",
    ));

    let mut term_id = 0;
    let mut uniq_terms = 0;
    let mut total_terms = 0;
    let mut frequent_terms = 0;
    let mut infrequent_terms = 0;
    while let Some(ciff::CiffRecord::PostingsList(plist)) = ciff_reader.next() {
        total_terms += 1;
        pb_plist.inc(1);
        if plist.postings.len() < min_len {
            infrequent_terms += 1;
            continue;
        }
        if plist.postings.len() >= cutoff_len {
            frequent_terms += 1;
            continue;
        }
        let postings = plist.get_postings();
        let mut doc_id: usize = 0;
        for posting in postings {
            doc_id += posting.get_docid() as usize;
            docs[doc_id].terms.push(term_id);
        }
        term_id += 1;
        uniq_terms += 1;
    }
    pb_plist.finish_and_clear();
    for doc in docs.iter_mut() {
        doc.terms.shrink_to_fit();
    }
    info!("forward index stats:");
    info!("\ttotal terms: {}", total_terms);
    info!("\tdiscarded frequent terms: {}", frequent_terms);
    info!("\tdiscarded infrequent terms: {}", infrequent_terms);
    info!("\tremaining terms: {}", term_id + 1);

    let forward_idx = Forward { docs, uniq_terms };

    if let Some(path) = output_path {
        info!("Saving forward index to file: {:?}", path.as_ref().to_str());
        let mut of = BufWriter::new(File::create(path).unwrap());
        serialize_into(&mut of, &forward_idx).unwrap(); 
    }

    Ok(forward_idx)
}
