
use crate::ciff;
use crate::forward::Doc;
use crate::ciff::proto::*;
use std::io::Write;

pub fn dump_order<P: AsRef<std::path::Path>>(docs: &[Doc], output_map: P) {
    let mut out = std::fs::File::create(output_map).expect("can not open output mapping file");
    for (new_id, doc) in docs.iter().enumerate() {
        std::writeln!(out, "{} {} {}", doc.org_id, new_id, doc.leaf_id).expect("can not write to output mapping file");
    }
}

pub fn remap_ciff<P: AsRef<std::path::Path>>(docs: &[Doc], input_ciff: P) -> Result<(), protobuf::ProtobufError> {
    // (1) create the reverse mapping
    let mut document_mapping: Vec<u32> = vec![0; docs.len()];
    for (new_id, doc) in docs.iter().enumerate() {
        document_mapping[doc.org_id as usize] = new_id as u32;
    }
    let ciff_file = std::fs::File::open(input_ciff).expect("can not open input ciff");
    let mut ciff_file = std::io::BufReader::new(ciff_file);
    let mut ciff_reader = ciff::Reader::new(&mut ciff_file).expect("can not open input ciff");

    let log_gaps: Vec<f64> = (0..256).map(|i| (i as f64).log2()).collect();
    let mut log_sum = 0.0;
    let mut num_postings: usize = 0;
    
    let mut record = ciff_reader.next();
    while let Some(ciff::CiffRecord::PostingsList(plist)) = record {
        let mut mapped_list = PostingsList::default();
        mapped_list.set_term(plist.get_term().to_string());
        mapped_list.set_df(plist.get_df());
        mapped_list.set_cf(plist.get_cf());
        let postings = plist.get_postings();
        let mut doc_id: usize = 0;
        for posting in postings {
            let mut mapped_posting = Posting::default();
            doc_id += posting.get_docid() as usize;
            let mapped_id = document_mapping[doc_id as usize];
            let tf = posting.get_tf();
            mapped_posting.set_docid(mapped_id as i32);
            mapped_posting.set_tf(tf);
            mapped_list.postings.push(mapped_posting);
        }
        mapped_list.postings.sort_by_key(|p| p.get_docid());
        let mut prev: i32 = 0;
        let mut first = true;
        for posting in mapped_list.postings.iter_mut() {
            let absolute_did = posting.get_docid();
            let dgap = absolute_did - prev;
            prev = absolute_did;
            if first {
                log_sum += ((absolute_did + 1) as f64).log2();
                first = false;
            } else {
                if dgap < 256 {
                    log_sum += unsafe { log_gaps.get_unchecked(dgap as usize) };
                } else {
                    log_sum += (dgap as f64).log2();
                }
            }
            num_postings += 1;
        }
        record = ciff_reader.next();
    }
    let after_bpi = log_sum / num_postings as f64;
    log::info!("LogGap after reorder: {:.3} BPI",after_bpi);
    
    Ok(())
}

pub fn rewrite_ciff<P: AsRef<std::path::Path>>(docs: &[Doc], input_ciff: P, output_ciff: P) -> Result<(), protobuf::ProtobufError> {
    // (1) create the reverse mapping
    let mut document_mapping: Vec<u32> = vec![0; docs.len()];
    for (new_id, doc) in docs.iter().enumerate() {
        document_mapping[doc.org_id as usize] = new_id as u32;
    }
    // (2) read the original ciff and remap ids, write the new ciff
    log::info!("writing to ciff file: {}",&output_ciff.as_ref().display());
    let mut ciff_outfile = std::fs::File::create(output_ciff).expect("can not create output ciff");
    let ciff_file = std::fs::File::open(input_ciff).expect("can not open input ciff");
    let mut ciff_file = std::io::BufReader::new(ciff_file);
    let mut ciff_reader = ciff::Reader::new(&mut ciff_file).expect("can not open input ciff");
    let num_msgs = ciff_reader.num_postings_lists() as u64;
    let pb_plist = indicatif::ProgressBar::new(num_msgs);
    pb_plist.set_draw_delta(1000);
    pb_plist.set_style(indicatif::ProgressStyle::default_bar().template(
        "plist_storage: {spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] ({pos}/{len}, ETA {eta}, SPEED: {per_sec})",
    ));
    let mut ciff_writer = ciff::Writer::from_reader(&mut ciff_outfile,&ciff_reader,"recursive graph bisection")?;

    struct CiffDoc {
        id: i32,
        external_id: String,
        length: i32,
    };

    let mut record = ciff_reader.next();
    while let Some(ciff::CiffRecord::PostingsList(plist)) = record {
        let mut mapped_list = PostingsList::default();
        mapped_list.set_term(plist.get_term().to_string());
        mapped_list.set_df(plist.get_df());
        mapped_list.set_cf(plist.get_cf());
        let postings = plist.get_postings();
        let mut doc_id: usize = 0;
        for posting in postings {
            let mut mapped_posting = Posting::default();
            doc_id += posting.get_docid() as usize;
            let mapped_id = document_mapping[doc_id as usize];
            let tf = posting.get_tf();
            mapped_posting.set_docid(mapped_id as i32);
            mapped_posting.set_tf(tf);
            mapped_list.postings.push(mapped_posting);
        }
        mapped_list.postings.sort_by_key(|p| p.get_docid());
        let mut prev: i32 = 0;
        for posting in mapped_list.postings.iter_mut() {
            let absolute_did = posting.get_docid();
            let dgap = absolute_did - prev;
            prev = absolute_did;
            posting.set_docid(dgap);
        }
        ciff_writer.write_postingslist(mapped_list)?;

        record = ciff_reader.next();

        pb_plist.inc(1);
    }
    pb_plist.finish_and_clear();

    let num_msgs = 2 * ciff_reader.num_docs() as u64;
    let pb_docs = indicatif::ProgressBar::new(num_msgs);
    pb_docs.set_draw_delta(1000);
    pb_docs.set_style(indicatif::ProgressStyle::default_bar().template(
        "docmeta_storage: {spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] ({pos}/{len}, ETA {eta}, SPEED: {per_sec})",
    ));

    let mut doc_records = Vec::with_capacity(docs.len());
    while let Some(ciff::CiffRecord::Document {
        doc_id,
        external_id,
        length,
    }) = record
    {
        let mapped_id = document_mapping[doc_id as usize];
        doc_records.push(CiffDoc {
            id: mapped_id as i32,
            external_id,
            length: length as i32,
        });

        record = ciff_reader.next();
        pb_plist.inc(1);
    }

    // (3) sort the doc records and write them in correct order
    doc_records.sort_by_key(|d| d.id);
    for doc in doc_records {
        ciff_writer.write_document(doc.id, doc.external_id, doc.length)?;
        pb_plist.inc(1);
    }
    pb_plist.finish_and_clear();

    ciff_writer.flush()?;

    Ok(())
}
