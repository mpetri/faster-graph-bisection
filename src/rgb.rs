use crate::forward::Doc;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSliceMut;
use std::cmp::Ordering::Equal;
use std::cmp::{min, max};

// XXX uses quickselect by default. If QUICKSELECT is not enabled,
// the sort-based method will be used. If PSORT is true, the sorts
// will be invoked in parallel. Otherwise, sequential sorts will be
// used.
const QUICKSELECT:bool = true;
const PSORT:bool = true;

// XXX does not use cooling by default. But if you enable this flag,
// cooling will be activated.
const COOLING:bool = false;

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

// Wrapper to various partitioning methods
// Default is to use floydrivest
// can specify `cooling` parameters to the `not_partitioned` call
// to allow higher tolerances on swaps
fn partition_quickselect(docs: &mut [Doc], k: usize, tolerance: f32) {

    // (1) Use quickselect
    //_quickselect(docs, k, 0, docs.len() - 1);
    
    // (2) Use dutch-flag quickselect
    //_dutch_flag_quickselect(docs, k, 0, docs.len() - 1);

    // (3) Use Floyd-Rivest with partitioned check
    if not_partitioned(docs, k, tolerance) {
        floydrivest(docs, k, 0, docs.len() - 1);
    }
}


// Checks if we have a partition already or not. A partition means that
// all values to the left of the median are stricly < median, and all values
// to the right are > median (based on document gains). Note that ties with
// the median are considered to be inconsequential to the computation
// Note: if tolerance is set, this will allow errors -- the higher the
// bound, the higher the tolerance for out-of-place values
fn not_partitioned(
    docs: &mut [Doc],
    median_idx: usize,
    tolerance: f32,
) -> bool {
    let median_val = docs[median_idx].gain;
    let (left, right) = docs.split_at(median_idx);
    left.iter().any(|d| (d.gain - tolerance) > median_val) || right.iter().any(|d| (d.gain + tolerance) < median_val)
}

// This function swaps two documents in the global document slice, indexed by `left` and `right`
// If the swap occurs across the median, we also need to flag that we need to update the degrees. 
// Otherwise, a simple pointer swap is all that is required
fn swap_docs(docs: &mut [Doc], left: usize, right: usize, median_idx: usize)  {

    // If the swap occurs entirely on one half, we don't need to update the degrees
    if (left < median_idx && right < median_idx) || (left >= median_idx && right >= median_idx) {
        docs.swap(left, right);
    } else {
        // Otherwise, we need to do the update...
        docs[left].leaf_id += 1; // Moved left2right --> increment by 1
        docs[right].leaf_id -= 1; // Moved right2left --> decrement by 1
        docs.swap(left, right);
    }
}


// Floyd Rivest algorithm to partition the median document by gains.
// Results in a modified document slice such that the median value is in its correct position, everything to the
// left has a gain < median, and everything to the right has a gain > median. Updates the degrees during the
// swapping process.
// Based on https://github.com/huonw/order-stat/blob/master/src/floyd_rivest.rs
// and https://github.com/frjnn/floydrivest/blob/master/src/lib.rs
fn floydrivest(docs: &mut [Doc], k: usize, in_left: usize, in_right: usize)
{
    let mut i: usize;
    let mut j: usize;
    let mut t_idx: usize;
    let mut left = in_left;
    let mut right = in_right;

    // Constants
    const MAX_PARTITION: usize = 600; // Arbitrary, used in the original paper
    const C: f32 = 0.5; // Arbitrary, used in the original paper

    while right > left {

        // If the partition we are looking at has more than MAX_PARTITION elements
        if right - left > MAX_PARTITION {

            // Use recursion on a sample of size s to get an estimate
            // for the (k - left + 1 )-th smallest elementh into a[k],
            // biased slightly so that the (k - left + 1)-th element is expected
            // to lie in the smallest set after partitioning
            let no_elements: f32 = (right - left + 1) as f32;
            let i: f32 = (k - left + 1) as f32;
            let z: f32 = no_elements.ln();
            let s: f32 = C * (z * (2.0 / 3.0)).exp();
            let sn: f32 = s / no_elements;
            let sd: f32 = C * (z * s * (1.0 - sn)).sqrt() * (i - no_elements * 0.5).signum();

            let isn: f32 = i * s / no_elements;
            let inner: f32 = k as f32 - isn + sd;
            let ll: usize = max(left, inner as usize);
            let rr: usize = min(right, (inner + s) as usize);
            floydrivest(docs, k, ll, rr);
        }

        // The following code partitions a[l : r] about t, it is similar to Hoare's
        // algorithm but it'll run faster on most machines since the subscript range
        // checking on i and j has been removed.
        
        i = left + 1;
        j = right - 1;
        
        swap_docs(docs, left, k, k);

        // if left doc > right doc, swap them
        if docs[left].gain >= docs[right].gain {
            swap_docs(docs, left, right, k);
            t_idx = right;
        } else {
            t_idx = left;
        }

        // Move left pointer up
        while docs[i].gain < docs[t_idx].gain {
            i+=1 ;
        }
        
        // Move right pointer down
        while docs[j].gain > docs[t_idx].gain {
            j-=1;
        }
      

        while i < j { 
            swap_docs(docs, i, j, k);
            i += 1;
            j -= 1;
            
            // Move left pointer up
            while docs[i].gain < docs[t_idx].gain {
                i+=1 ;
            }
        
            // Move right pointer down
            while docs[j].gain > docs[t_idx].gain {
                j-=1;
            }
        }
       
        if left == t_idx {
            swap_docs(docs, left, j, k);
        } else {
            j += 1;
            swap_docs(docs, j, right, k);
        }
        if j <= k {
            left = j + 1;
        }
        if k <= j {
                right = j.saturating_sub(1);
        }
    }
}

// Uses quickselect algorithm with Hoare's partitioning to partition the median document by gains
// Results in a modified document slice such that the median value is in its correct position, everything to the
// left has a gain < median, and everything to the right has a gain > median. Updates the degrees during the
// swapping process.
fn _quickselect(docs: &mut [Doc], k: usize, in_left: usize, in_right: usize)
{
    let mut left = in_left;
    let mut right = in_right;

    loop {

        if left == right {
            break;
        }

        let original_low = left;
        let original_high = right;
        let pivot = (left + right + 1) / 2;
        let pivot_gain = docs[pivot].gain;
        
        // Hoare partitioning 
        loop {
            // Move left pointer up
            while docs[left].gain < pivot_gain {
                left+=1;
            }
            // Move right pointer down
            while docs[right].gain > pivot_gain {
                right-=1;
            }
            // Pointers crossed, we can exit here and check the partitions
            if left >= right {
                break;
            }
            swap_docs(docs, left, right, k);
            left += 1;
            right -= 1;
        }
    
        // We'll look at docs[low, left-1]
        if k < left {
            right = left - 1;
            left = original_low;
        }
        // We'll look at docs[left, max]
        else {
            right = original_high;
        }
    
    } 
} 

// Uses a modified quickselect algorithm to with dutch flag partitioning
// Results in a modified document slice such that the median value is in its correct position, everything to the
// left has a gain < median, and everything to the right has a gain > median. Updates the degrees during the
// swapping process.
fn _dutch_flag_quickselect(docs: &mut [Doc], k: usize, in_left: usize, in_right: usize)
{

    let mut bot = 0;
    let mut size = in_right;

    while size > 1 {
        
        let pivot = k;
        let pivot_gain = docs[pivot].gain;
        
        let mut nxt = in_left;
        let mut feq = in_left;
        let mut fgt = size + bot;

        while nxt < fgt {
            
            if docs[nxt].gain < pivot_gain {
                swap_docs(docs, feq, nxt, k);
                feq +=1;
                nxt +=1;
            } else if docs[nxt].gain > pivot_gain {
                fgt -= 1;
                swap_docs(docs, nxt, fgt, k);
            } else {
                nxt += 1;
            }
            
        }

        if docs[k].gain == docs[feq].gain {
            break;
        } else {
            if k < feq {
                size = feq - bot
            } else {
                size = size - (fgt-bot);
                bot = fgt;
            }
        }
    } 
} 


// This method will rip through the documents vector
// and update the degrees of documents which swapped
fn fix_degrees(
    docs: &[Doc],
    left_degs: &mut [i32],
    right_degs: &mut [i32],
) -> usize {
    let mut num_swaps = 0;
    for doc in docs.iter() {
        // Doc went right to left
        if doc.leaf_id == -1 {
            for (term, _) in &doc.postings {
                left_degs[*term as usize] += 1;
                right_degs[*term as usize] -= 1;
            }
            num_swaps += 1;
        } 
        // Moved left to right
        else if doc.leaf_id == 1 {
            for (term, _) in &doc.postings {
                left_degs[*term as usize] -= 1;
                right_degs[*term as usize] += 1;
            }
            num_swaps += 1;
        } 
    }
    num_swaps
}


// This is an alternative swap method, which assumes negative
// numbers want to go left, and positive numbers want to go
// right. It assumes that left comes in sorted descending (so documents
// wanting to move right appear first), and right comes in sorted
// ascending (so documents wanting to move left appear first).
// Thus, a swap will take place if (d_left - d_right) > tolerance.
fn swap_documents(
    left: &mut [Doc],
    right: &mut [Doc],
    left_degs: &mut [i32],
    right_degs: &mut [i32],
    tolerance: f32,
) -> usize {
    let mut num_swaps = 0;
    for (l, r) in left.iter_mut().zip(right.iter_mut()) {
        if l.gain - r.gain > tolerance {
            for (term, _) in &l.postings {
                left_degs[*term as usize] -= 1;
                right_degs[*term as usize] += 1;
            }
            for (term, _) in &r.postings {
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
        for (term, _) in &doc.postings {
            degrees[*term as usize] += 1;
        }
    }
    degrees
}

// Two convenience functions used to compute the left/right degrees of
// a document slice in parallel
fn compute_degrees_l(docs: &[Doc], num_terms: usize) -> Vec<i32> {
    let (left, _) = docs.split_at(docs.len() / 2);
    compute_degrees(left, num_terms)    
}

fn compute_degrees_r(docs: &[Doc], num_terms: usize) -> Vec<i32> { 
    let (_, right) = docs.split_at(docs.len() / 2);
    compute_degrees(right, num_terms)    
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
fn compute_move_gains_default_l2r(
    from: &mut [Doc],
    log2_from: f32,
    log2_to: f32,
    fdeg: &[i32],
    tdeg: &[i32],
) {
    from.par_iter_mut().for_each(|doc| {
        let mut doc_gain = 0.0f32;
        for (term, _) in &doc.postings {
            let from_deg = fdeg[*term as usize];
            let to_deg = tdeg[*term as usize];
            let term_gain = expb(log2_from, log2_to, from_deg, to_deg)
                - expb(log2_from, log2_to, from_deg - 1, to_deg + 1);
            doc_gain += term_gain;
        }
        doc.gain = doc_gain;
        doc.leaf_id = 0;
    });
}

fn compute_move_gains_default_r2l(
    from: &mut [Doc],
    log2_from: f32,
    log2_to: f32,
    fdeg: &[i32],
    tdeg: &[i32],
) {
    from.par_iter_mut().for_each(|doc| {
        let mut doc_gain = 0.0f32;
        for (term, _) in &doc.postings {
            let from_deg = fdeg[*term as usize];
            let to_deg = tdeg[*term as usize];
            let term_gain = expb(log2_from, log2_to, from_deg, to_deg)
                - expb(log2_from, log2_to, from_deg - 1, to_deg + 1);
            doc_gain -= term_gain;
        }
        doc.gain = doc_gain;
        doc.leaf_id = 0;
    });
}



// Computes gains using the first approximation, saving two log calls
fn compute_move_gains_a1_l2r(
    from: &mut [Doc],
    log2_from: f32,
    log2_to: f32,
    fdeg: &[i32],
    tdeg: &[i32],
) {
    from.par_iter_mut().for_each(|doc| {
        let mut doc_gain = 0.0f32;
        for (term, _) in &doc.postings {
            let from_deg = fdeg[*term as usize];
            let to_deg = tdeg[*term as usize];
            let term_gain = approx_one_a(log2_to, log2_from, to_deg, from_deg);
            doc_gain += term_gain;
        }
        doc.gain = doc_gain;
        doc.leaf_id = 0;
    });
}
// Computes gains using the first approximation, saving two log calls
fn compute_move_gains_a1_r2l(
    from: &mut [Doc],
    log2_from: f32,
    log2_to: f32,
    fdeg: &[i32],
    tdeg: &[i32],
) {
    from.par_iter_mut().for_each(|doc| {
        let mut doc_gain = 0.0f32;
        for (term, _) in &doc.postings {
            let from_deg = fdeg[*term as usize];
            let to_deg = tdeg[*term as usize];
            let term_gain = approx_one_a(log2_to, log2_from, to_deg, from_deg);
            doc_gain -= term_gain; // Note the negative sign here
        }
        doc.gain = doc_gain;
        doc.leaf_id = 0;
    });
}



// Computes gains using the second approximation, saving four log calls
// Since it's symmetric, we don't need an l2r or r2l function
fn compute_move_gains_a2(
    from: &mut [Doc],
    log2_from: f32,
    log2_to: f32,
    fdeg: &[i32],
    tdeg: &[i32],
) {
    from.par_iter_mut().for_each(|doc| {
        let mut doc_gain = 0.0f32;
        for (term, _) in &doc.postings {
            let from_deg = fdeg[*term as usize];
            let to_deg = tdeg[*term as usize];
            let term_gain = approx_two_s(log2_to, log2_from, to_deg, from_deg);
            doc_gain += term_gain;
        }
        doc.gain = doc_gain;
        doc.leaf_id = 0; 
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
                    compute_move_gains_default_l2r(&mut left, log2_left, log2_right, &ldeg, &rdeg);
                });
                s.spawn(|_| {
                    compute_move_gains_default_r2l(&mut right, log2_right, log2_left, &rdeg, &ldeg);
                });
            });
        }
        // (1) -- First approximation
        "approx_1" => {
            rayon::scope(|s| {
                s.spawn(|_| {
                    compute_move_gains_a1_l2r(&mut left, log2_left, log2_right, &ldeg, &rdeg);
                });
                s.spawn(|_| {
                    compute_move_gains_a1_r2l(&mut right, log2_right, log2_left, &rdeg, &ldeg);
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
                    compute_move_gains_a2(&mut right, log2_right, log2_left, &ldeg, &rdeg);
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
    mut docs: &mut [Doc],
    num_terms: usize,
    iterations: usize,
) {

    // compute degrees in left and right partition for each term
    let (mut left_deg, mut right_deg) = rayon::join(
        || compute_degrees_l(&docs, num_terms),
        || compute_degrees_r(&docs, num_terms), 
    );

    for _iter in 0..iterations {
    
        // Minselect method
        if QUICKSELECT { 
            
            // Split in half and compute gains
            let (mut left, mut right) = docs.split_at_mut(docs.len() / 2);
            compute_gains(&mut left, &mut right, &left_deg[..], &right_deg[..]);

            // Use quickselect to partition the global doc slice into < median and > median
            let median_idx = docs.len() / 2;

            if COOLING {
                partition_quickselect(&mut docs, median_idx, (_iter as f32) * 0.5); // Simulated annealing
            } else {
                partition_quickselect(&mut docs, median_idx, 0.0);
            }

            // Go through swapped documents and fix the term degrees
            let nswaps = fix_degrees(docs, &mut left_deg[..], &mut right_deg[..]);
            if nswaps == 0 {
                break;
            }

        } 
        
        // Regular sort+swap method. Note that this uses a parallel sort if PSORT is set to `true`.
        else {

            // Split in half and compute gains
            let (mut left, mut right) = docs.split_at_mut(docs.len() / 2);
            compute_gains(&mut left, &mut right, &left_deg[..], &right_deg[..]);

            // Use parallel sorts on each half to re-arrange documents by their gains
            if PSORT {
                left.par_sort_by(|a, b| b.gain.partial_cmp(&a.gain).unwrap_or(Equal)); // Sort gains high to low
                right.par_sort_by(|a, b| a.gain.partial_cmp(&b.gain).unwrap_or(Equal)); // Sort gains low to high
            } else {
                left.sort_by(|a, b| b.gain.partial_cmp(&a.gain).unwrap_or(Equal)); // Sort gains high to low
                right.sort_by(|a, b| a.gain.partial_cmp(&b.gain).unwrap_or(Equal)); // Sort gains low to high
            }

            // Go through and swap documents between partitions
            let nswaps = if COOLING {
                swap_documents(&mut left, &mut right, &mut left_deg[..], &mut right_deg[..], _iter as f32) // Simulated annealing
            } else {
                swap_documents(&mut left, &mut right, &mut left_deg[..], &mut right_deg[..], 0.0)
            };

            if nswaps == 0 {
                break;
            }
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

    // (1) swap around docs between the two halves based on move gains
    process_partitions(docs, num_terms, iterations);

    let (mut left, mut right) = docs.split_at_mut(docs.len() / 2);
 
    // (2) recurse left and right
    rayon::scope(|s| {
        s.spawn(|_| {
            recursive_graph_bisection(&mut left, num_terms, iterations, min_partition_size, max_depth, depth+1, sort_leaf);
        });
        s.spawn(|_| {
            recursive_graph_bisection(&mut right, num_terms, iterations, min_partition_size, max_depth, depth+1, sort_leaf);
        });
    });
}
