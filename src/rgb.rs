use crate::forward::Doc;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSliceMut;
use std::cmp::Ordering::Equal;


// log2(e) approximation 
const LOG2_E: f32 = 1.44269504089;

// The following pre-computes a table of log2 values to improve
// efficiency. Note that log2(0) is defined to be 0 for convenience.
const LOG_PRECOMP_LIMIT: i32 = 4096;

lazy_static::lazy_static! {
    static ref LOG2: Vec<f32> = (0..LOG_PRECOMP_LIMIT).map(|i|
        match i {
            0 => 0.0,
            _ => (i as f32).log2()
        }).collect();
}

#[inline]
fn cmp_log2(x: i32) -> f32 {
    if x < LOG_PRECOMP_LIMIT {
        LOG2[x as usize]
    } else {
        (x as f32).log2()
    }
}

// This is the default document swap method. It assumes positive gains
// on both left and right indicate a desire to swap, and so swaps while
// the sum of these gains is > 0
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

// Computes the sum of the term degrees across a slice of documents
fn compute_degrees(docs: &[Doc], num_terms: usize) -> Vec<i32> {
    let mut degrees = vec![0; num_terms];
    for doc in docs {
        for term in &doc.terms {
            degrees[*term as usize] += 1;
        }
    }
    degrees
}

// This is the original cost function: Asymmetric
fn expb(log_from: f32, log_to: f32, deg1: i32, deg2: i32) -> f32 {
    let d1f = deg1 as f32;
    let d2f = deg2 as f32;

    let a = d1f * log_from;
    let b = d1f * cmp_log2(deg1 + 1);

    let c = d2f * log_to;
    let d = d2f * cmp_log2(deg2 + 1);
    a - b + c - d
}

// This is the first approximation: Asymmetric
fn approx_one_a(_log_to: f32, _log_from: f32, deg_to: i32, deg_from: i32) -> f32 {
    cmp_log2(deg_to + 2) - cmp_log2(deg_from) - LOG2_E / (1.0 + deg_to as f32)
}

// This is the second approximation: Symmetric
fn approx_two_s(_log_to: f32, _log_from: f32, deg_to: i32, deg_from: i32) -> f32 {
    cmp_log2(deg_to) - cmp_log2(deg_from)
}

// This function will compute gains using the baseline approach
fn compute_move_gains_default(
    from: &mut [Doc],
    log2_from: f32,
    log2_to: f32,
    fdeg: &[i32],
    tdeg: &[i32],
) {
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

// Computes gains using the first approximation, saving two log calls
fn compute_move_gains_a1(
    from: &mut [Doc],
    log2_from: f32,
    log2_to: f32,
    fdeg: &[i32],
    tdeg: &[i32],
) {
    from.par_iter_mut().for_each(|doc| {
        let mut doc_gain = 0.0f32;
        for term in &doc.terms {
            let from_deg = fdeg[*term as usize];
            let to_deg = tdeg[*term as usize];
            let term_gain = approx_one_a(log2_to, log2_from, to_deg, from_deg);
            doc_gain += term_gain;
        }
        doc.gain = doc_gain;
    });
}

// Computes gains using the second approximation, saving four log calls
fn compute_move_gains_a2(
    from: &mut [Doc],
    log2_from: f32,
    log2_to: f32,
    fdeg: &[i32],
    tdeg: &[i32],
) {
    from.par_iter_mut().for_each(|doc| {
        let mut doc_gain = 0.0f32;
        for term in &doc.terms {
            let from_deg = fdeg[*term as usize];
            let to_deg = tdeg[*term as usize];
            let term_gain = approx_two_s(log2_to, log2_from, to_deg, from_deg);
            doc_gain += term_gain;
        }
        doc.gain = doc_gain;
    });
}

// This function is a wrapper to the correct gain function, which is specified at compile-time
// using the `GAIN` environment variable
fn compute_gains(mut left: &mut [Doc], mut right: &mut [Doc], ldeg: &[i32], rdeg: &[i32]) {
    let gain_func: Option<&'static str> = std::option_env!("GAIN");
    let gain_func = gain_func.unwrap();

    let log2_left = 0.0;
    let log2_right = 0.0;

    // (0) -- Default gain
    match gain_func {
        "default" => {
            let log2_left = (left.len() as f32).log2();
            let log2_right = (right.len() as f32).log2();
            rayon::scope(|s| {
                s.spawn(|_| {
                    compute_move_gains_default(&mut left, log2_left, log2_right, &ldeg, &rdeg);
                });
                s.spawn(|_| {
                    compute_move_gains_default(&mut right, log2_right, log2_left, &rdeg, &ldeg);
                });
            });
        }
        // (1) -- First approximation
        "approx_1" => {
            rayon::scope(|s| {
                s.spawn(|_| {
                    compute_move_gains_a1(&mut left, log2_left, log2_right, &ldeg, &rdeg);
                });
                s.spawn(|_| {
                    compute_move_gains_a1(&mut right, log2_right, log2_left, &rdeg, &ldeg);
                });
            });
        }

        // (2) -- Second approximation: Note that the right hand call needs to have the parameters
        // reversed for rdeg and ldeg if we opt for the `sorting` computation mode instead of the
        // `swapping` mode
        "approx_2" => {
            rayon::scope(|s| {
                s.spawn(|_| {
                    compute_move_gains_a2(&mut left, log2_left, log2_right, &ldeg, &rdeg);
                });
                s.spawn(|_| {
                    compute_move_gains_a2(&mut right, log2_right, log2_left, &rdeg, &ldeg);
                });
            });
        }

        // Should be unreachable...
        _ => {
            log::info!("Error: Couldn't match the gain function.");
        }
    }
}

// The heavy lifting -- core logic for the BP process
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

    // for each iteration, compute the gains, swap stuff, and check if we can terminate
    for _iter in 0..iterations {
        compute_gains(&mut left, &mut right, &left_deg[..], &right_deg[..]);

        // parallel sort by decreasing move gain
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
    sort_leaf: bool,
) {
    // recursion end?
    if docs.len() <= min_partition_size || depth > max_depth {
        // Sort leaf by input identifier
        if sort_leaf {
            docs.sort_by(|a, b| a.org_id.cmp(&b.org_id));
        }
        return;
    }

    // (1) split into two halfs
    let (mut left, mut right) = docs.split_at_mut(docs.len() / 2);

    // (2) swap around docs between the two halves based on move gains
    process_partitions(left, right, num_terms, iterations);

    // (3) recurse left and right
    rayon::scope(|s| {
        s.spawn(|_| {
            recursive_graph_bisection(&mut left, num_terms, iterations, min_partition_size, max_depth, depth+1, sort_leaf);
        });
        s.spawn(|_| {
            recursive_graph_bisection(&mut right, num_terms, iterations, min_partition_size, max_depth, depth+1, sort_leaf);
        });
    });
}
