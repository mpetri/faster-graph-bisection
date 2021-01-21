use crate::forward::Doc;
use rand::seq::SliceRandom;
use rand::thread_rng;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSliceMut;
use std::cmp::Ordering::Equal;

const LOG_PRECOMP_LIMIT: i32 = 4096;

lazy_static::lazy_static! {
    static ref LOG2: Vec<f32> = (0..LOG_PRECOMP_LIMIT).map(|i| (i as f32).log2()).collect();
}

#[inline]
fn cmp_log2(x: i32) -> f32 {
    if x < LOG_PRECOMP_LIMIT {
        LOG2[x as usize]
    } else {
        (x as f32).log2()
    }
}

fn swap_documents(
    left: &mut [Doc],
    right: &mut [Doc],
    left_degs: &mut [i32],
    right_degs: &mut [i32],
) -> usize {
    let mut num_swaps = 0;
    for (l, r) in left.iter_mut().zip(right.iter_mut()) {
        if l.gain + r.gain > 0.0 {
            for term in &l.terms {
                left_degs[*term as usize] -= 1;
                right_degs[*term as usize] += 1;
            }
            for term in &r.terms {
                left_degs[*term as usize] += 1;
                right_degs[*term as usize] -= 1;
            }
            std::mem::swap(l, r);
            num_swaps += 1;
        } else {
            break;
        }
    }
    num_swaps
}

fn compute_degrees(docs: &[Doc], num_terms: usize) -> Vec<i32> {
    let mut degrees = vec![0; num_terms];
    for doc in docs {
        for term in &doc.terms {
            degrees[*term as usize] += 1;
        }
    }
    degrees
}

fn expb(log_from: f32, log_to: f32, deg1: i32, deg2: i32) -> f32 {
    let d1f = deg1 as f32;
    let d2f = deg2 as f32;

    let a = d1f * log_from;
    let b = d1f * cmp_log2(deg1 + 1);

    let c = d2f * log_to;
    let d = d2f * cmp_log2(deg2 + 1);
    a - b + c - d
}

fn compute_move_gains(from: &mut [Doc], log2_from: f32, log2_to: f32, fdeg: &[i32], tdeg: &[i32]) {
    from.par_iter_mut().for_each(|doc| {
        let mut doc_gain = 0.0f32;
        for term in &doc.terms {
            let from_deg = fdeg[*term as usize];
            let to_deg = tdeg[*term as usize];
            let term_gain = expb(log2_from, log2_to, from_deg, to_deg)
                - expb(log2_from, log2_to, from_deg - 1, to_deg + 1);
            doc_gain += term_gain;
        }
        doc.gain = doc_gain;
    });
}

fn compute_gains(mut left: &mut [Doc], mut right: &mut [Doc], ldeg: &[i32], rdeg: &[i32]) {
    let log2_left = (left.len() as f32).log2();
    let log2_right = (right.len() as f32).log2();
    rayon::scope(|s| {
        s.spawn(|_| {
            compute_move_gains(&mut left, log2_left, log2_right, &ldeg, &rdeg);
        });
        s.spawn(|_| {
            compute_move_gains(&mut right, log2_right, log2_left, &rdeg, &ldeg);
        });
    });
}

fn process_partitions(
    mut left: &mut [Doc],
    mut right: &mut [Doc],
    num_terms: usize,
    iterations: usize,
) {
    // compute degrees in left and right partition for each term
    let (mut left_deg, mut right_deg) = rayon::join(
        || compute_degrees(&left, num_terms),
        || compute_degrees(&right, num_terms),
    );

    for _iter in 0..iterations {
        compute_gains(&mut left, &mut right, &left_deg[..], &right_deg[..]);

        // sort by decreasing move gain
        left.par_sort_by(|a, b| b.gain.partial_cmp(&a.gain).unwrap_or(Equal));
        right.par_sort_by(|a, b| b.gain.partial_cmp(&a.gain).unwrap_or(Equal));
        // swap docs around while swap gain > 0
        let nswaps = swap_documents(&mut left, &mut right, &mut left_deg[..], &mut right_deg[..]);
        if nswaps == 0 {
            break;
        }
    }
}

pub fn recursive_graph_bisection(
    docs: &mut [Doc],
    num_terms: usize,
    iterations: usize,
    min_partition_size: usize,
    max_depth: usize,
    depth: usize,
    progress: indicatif::ProgressBar,
    sort_leaf: bool,
    slice_id: i32,
) {
    // recursion end?
    if docs.len() <= min_partition_size || depth > max_depth {
        // Sort leaf by input identifier
        if sort_leaf {
            docs.sort_by(|a, b| a.org_id.cmp(&b.org_id));
        }
        // Set up leaf identifiers
        for doc in docs.iter_mut() {
            doc.leaf_id = slice_id;
        }
        progress.inc(docs.len() as u64);
        return;
    }
    // (1) random shuffle first
    let mut rng = thread_rng();
    docs.shuffle(&mut rng);

    // (2) split into two halfs
    let (mut left, mut right) = docs.split_at_mut(docs.len() / 2);

    // (3) swap around docs between the two halves based on move gains
    process_partitions(left, right, num_terms, iterations);

    // (4) recurse left and right
    rayon::scope(|s| {
        let progress_left = progress.clone();
        let progress_right = progress.clone();
        s.spawn(|_| {
            recursive_graph_bisection(&mut left, num_terms, iterations, min_partition_size, max_depth, depth+1, progress_left, sort_leaf, slice_id * 2);
        });
        s.spawn(|_| {
            recursive_graph_bisection(&mut right, num_terms, iterations, min_partition_size, max_depth, depth+1, progress_right, sort_leaf, slice_id * 2 + 1);
        });
    });
}
