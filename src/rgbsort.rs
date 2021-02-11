use crate::forward::Doc;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSliceMut;
use std::cmp::Ordering::Equal;
use std::cmp::{max, min, Ordering};

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

// Wrapper for floydrivest
fn partition_floyd_rivest(docs: &mut [Doc], nth_el: usize,)  {
    let mut cmp = |a: &Doc, b: &Doc| a.gain.partial_cmp(&b.gain).unwrap_or(Equal);
    if not_partitioned(docs, nth_el) {
        floydrivest(docs, nth_el, 0, docs.len() - 1, &mut cmp, nth_el);
    }
}

fn partition_quickselect(docs: &mut [Doc], nth_el: usize) {
    quickselect(docs, nth_el, 0, docs.len() - 1, nth_el);
}

// Checks if we have a partition already or not. A partition means that
// all values to the left of the median are stricly < median, and all values
// to the right are > median (based on document gains). Note that ties with
// the median are considered to be inconsequential to the computation
fn not_partitioned(
    docs: &mut [Doc],
    median_idx: usize,
) -> bool {
    let median_val = docs[median_idx].gain;
    let (left, right) = docs.split_at(median_idx);
    left.iter().any(|d| d.gain > median_val) || right.iter().any(|d| d.gain < median_val)
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


// Heavy lifting for median partition. Uses modified Floyd Rivest algorithm to partition the median document by gains
// Results in a modified document slice such that the median value is in its correct position, everything to the
// left has a gain < median, and everything to the right has a gain > median. Updates the degrees during the
// swapping process.
// Based on https://github.com/huonw/order-stat/blob/master/src/floyd_rivest.rs
// and https://github.com/frjnn/floydrivest/blob/master/src/lib.rs
fn floydrivest<F>(docs: &mut [Doc], nth_el: usize, mut left: usize, mut right: usize, cmp: &mut F, median_idx: usize)
where
    F: FnMut(&Doc, &Doc) -> Ordering,
{
    let mut i: usize;
    let mut j: usize;
    let mut t_idx: usize;
    while right > left {
        if right - left > 600 {
            // Use recursion on a sample of size s to get an estimate
            // for the (nth_el - left + 1 )-th smallest elementh into a[nth_el],
            // biased slightly so that the (nth_el - left + 1)-th element is expected
            // to lie in the smallest set after partitioning.
            let n: f64 = (right - left + 1) as f64;
            let i: f64 = (nth_el - left + 1) as f64;
            let z: f64 = n.ln();
            let s: f64 = 0.5 * (z * (2.0 / 3.0)).exp();
            let sn: f64 = s / n;
            let sd: f64 = 0.5 * (z * s * (1.0 - sn)).sqrt() * (i - n * 0.5).signum();

            let isn: f64 = i * s / n;
            let inner: f64 = nth_el as f64 - isn + sd;
            let ll: usize = max(left, inner as usize);
            let rr: usize = min(right, (inner + s) as usize);
            floydrivest(docs, nth_el, ll, rr, cmp, median_idx);
        }
        // The following code partitions a[l : r] about t, it is similar to Hoare's
        // algorithm but it'll run faster on most machines since the subscript range
        // checking on i and j has been removed.
        
        i = left + 1;
        j = right - 1;
        
        swap_docs(docs, left, nth_el, median_idx);

        if cmp(&docs[left], &docs[right]) != Ordering::Less {
            swap_docs(docs, left, right, median_idx);
            t_idx = right;
        } else {
            t_idx = left;
        }


        while cmp(&docs[i], &docs[t_idx]) == Ordering::Less { i +=1 }
        while cmp(&docs[j], &docs[t_idx]) == Ordering::Greater { j -= 1 }
       
        while i < j { 
            swap_docs(docs, i, j, median_idx);
            i += 1;
            j -= 1;
            while cmp(&docs[i], &docs[t_idx]) == Ordering::Less { i +=1 }
            while cmp(&docs[j], &docs[t_idx]) == Ordering::Greater { j -= 1 }
        }
       
        if left == t_idx {
            swap_docs(docs, left, j, median_idx);
        } else {
            j += 1;
            swap_docs(docs, j, right, median_idx);
        }
        if j <= nth_el {
            left = j + 1;
        }
        if nth_el <= j {
                right = j.saturating_sub(1);
        }
    }
}

// Heavy lifting for median partition. Uses modified quickselect algorithm to partition the median document by gains
// Results in a modified document slice such that the median value is in its correct position, everything to the
// left has a gain < median, and everything to the right has a gain > median. Updates the degrees during the
// swapping process.
fn quickselect(docs: &mut [Doc], nth_el: usize, mut left: usize, mut right: usize, median_idx: usize)
{

    loop {
        if left == right {
            break;
        }
        let mut pivot = left + (right - left) / 2; // XXX Use the median for now
        
        // lumunto partitioning
        let mut pivot_gain = docs[pivot].gain;
        swap_docs(docs, pivot, right, median_idx);
        let mut store_idx = left;
        for i in left..right {
            if docs[i].gain < pivot_gain {
                swap_docs(docs, store_idx, i, median_idx);
                store_idx += 1;
            }
        }
        swap_docs(docs, right, store_idx, median_idx);
        pivot = store_idx;


        if pivot == nth_el {
            break;
        } else if nth_el < pivot {
            right = pivot - 1;
        } else {
            left = pivot + 1;
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
            for term in &doc.terms {
                left_degs[*term as usize] += 1;
                right_degs[*term as usize] -= 1;
            }
            num_swaps += 1;
        } 
        // Moved left to right
        else if doc.leaf_id == 1 {
            for term in &doc.terms {
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
// Thus, a swap will take place if (d_left - d_right) > 0.0.
fn swap_documents(
    left: &mut [Doc],
    right: &mut [Doc],
    left_degs: &mut [i32],
    right_degs: &mut [i32],
) -> usize {
    let mut num_swaps = 0;
    for (l, r) in left.iter_mut().zip(right.iter_mut()) {
        if l.gain - r.gain > 0.0 {
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
        for term in &doc.terms {
            let from_deg = fdeg[*term as usize];
            let to_deg = tdeg[*term as usize];
            let term_gain = expb(log2_from, log2_to, from_deg, to_deg)
                - expb(log2_from, log2_to, from_deg - 1, to_deg + 1);
            doc_gain += term_gain;
        }
        doc.gain = doc_gain;
        doc.leaf_id = 0; //XXX
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
        for term in &doc.terms {
            let from_deg = fdeg[*term as usize];
            let to_deg = tdeg[*term as usize];
            let term_gain = expb(log2_from, log2_to, from_deg, to_deg)
                - expb(log2_from, log2_to, from_deg - 1, to_deg + 1);
            doc_gain -= term_gain;
        }
        doc.gain = doc_gain;
        doc.leaf_id = 0; //XXX
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
        for term in &doc.terms {
            let from_deg = fdeg[*term as usize];
            let to_deg = tdeg[*term as usize];
            let term_gain = approx_one_a(log2_to, log2_from, to_deg, from_deg);
            doc_gain += term_gain;
        }
        doc.gain = doc_gain;
        doc.leaf_id = 0; //XXX
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
        for term in &doc.terms {
            let from_deg = fdeg[*term as usize];
            let to_deg = tdeg[*term as usize];
            let term_gain = approx_one_a(log2_to, log2_from, to_deg, from_deg);
            doc_gain -= term_gain; // Note the negative sign here
        }
        doc.gain = doc_gain;
        doc.leaf_id = 0; //XXX
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
        for term in &doc.terms {
            let from_deg = fdeg[*term as usize];
            let to_deg = tdeg[*term as usize];
            let term_gain = approx_two_s(log2_to, log2_from, to_deg, from_deg);
            doc_gain += term_gain;
        }
        doc.gain = doc_gain;
        doc.leaf_id = 0; //XXX
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
    depth: usize,
) {
    // compute degrees in left and right partition for each term
    let (mut left_deg, mut right_deg) = rayon::join(
        || compute_degrees_l(&docs, num_terms),
        || compute_degrees_r(&docs, num_terms), 
    );

    for _iter in 0..iterations {
    
        if true { // quickselect

            //let start_gains = std::time::Instant::now();
            // Split in half and compute gains
            {
                let (mut left, mut right) = docs.split_at_mut(docs.len() / 2);
                compute_gains(&mut left, &mut right, &left_deg[..], &right_deg[..]);
            }
            //let gains_time = start_gains.elapsed().as_micros();
            //println!("gains it = {} depth = {} time_micro = {}", _iter, depth, gains_time);
 
            //let start_gains = std::time::Instant::now();
            let median_idx = docs.len() / 2;
            partition_floyd_rivest(&mut docs, median_idx);
            //partition_quickselect(&mut docs, median_idx);
            //for doc in docs.iter() {
            //    println!("{}", doc.gain);
            //}
            //let gains_time = start_gains.elapsed().as_micros();
            //println!("quickselect it = {} depth = {} time_micro = {} total_swaps = {}", _iter, depth, gains_time, nswaps);
            let nswaps = fix_degrees(docs, &mut left_deg[..], &mut right_deg[..]);
            if nswaps == 0 {
                break;
            }
        } else {

            // We do the split here, because if we decide later that we want to do something
            // with the whole document partition, we can refactor easily...
            let start_gains = std::time::Instant::now();
            let (mut left, mut right) = docs.split_at_mut(docs.len() / 2);
            compute_gains(&mut left, &mut right, &left_deg[..], &right_deg[..]);
            let gains_time = start_gains.elapsed().as_micros();
            println!("gains it = {} depth = {} time_micro = {}", _iter, depth, gains_time);
   
            let start_gains = std::time::Instant::now();
            left.par_sort_by(|a, b| b.gain.partial_cmp(&a.gain).unwrap_or(Equal)); // Sort gains high to low
            right.par_sort_by(|a, b| a.gain.partial_cmp(&b.gain).unwrap_or(Equal)); // Sort gains low to high
            let gains_time = start_gains.elapsed().as_micros();
            println!("sort it = {} depth = {} time_micro = {}", _iter, depth, gains_time);
 
            let start_gains = std::time::Instant::now();
            let nswaps = swap_documents(&mut left, &mut right, &mut left_deg[..], &mut right_deg[..]);
            let gains_time = start_gains.elapsed().as_micros();
            println!("swaps it = {} depth = {} time_micro = {} total_swaps = {}", _iter, depth, gains_time, nswaps);
 
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
    process_partitions(docs, num_terms, iterations, depth);

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
