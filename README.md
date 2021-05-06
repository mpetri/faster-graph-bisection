# Faster Recursive Graph Bisection

This repo contains the code corresponding to the SIGIR 2021 short paper 
*Faster Index Reordering with Bipartite Graph Partitioning* by Joel Mackenzie,
Matthias Petri, and Alistair Moffat. 

## Citation Information
If you use this code in your own work or research, please consider citing
our work:
```
@inproceedings{mpm21-sigir,
 title = {Faster Index Reordering with Bipartite Graph Partitioning},
 author = {J. Mackenzie and M. Petri and A. Moffat},
 booktitle = {Proc. SIGIR},
 pages = {To Appear},
 year = {2021},
}
```
The paper can be found at the following DOI: https://doi.org/10.1145/3404835.3462991

## Acknowledgements
This work was built on previous work from Dhulipala et. al:
[Compressing Graphs and Indexes with Recursive Graph Bisection](http://www.kdd.org/kdd2016/papers/files/rpp0883-dhulipalaAemb.pdf), 
[ACM Proceedings](https://dl.acm.org/citation.cfm?id=2939862).

We also used the reproducibility study from Mackenzie et. al:
[Compressing Inverted Indexes with Recursive Graph Bisection: A Reproducibility Study](http://engineering.nyu.edu/~suel/papers/bp-ecir19.pdf),
[Springer Proceedings](https://link.springer.com/chapter/10.1007/978-3-030-15712-8_22).

Our codebase is based on the implementation found in the [PISA](https://github.com/pisa-engine/pisa/) search engine, which
corresponds to the reproducibility study discussed above. The codebase works with the
[Common Index File Format](https://github.com/osirrc/ciff), an open-source index exchange format for information
retrieval experimentation.


## Building the code
You can build the code using Cargo:
```
cargo build --release
```

However, if you follow the command above, running the code will give an error:
```
./target/release/create-rgb
03:09:14 [INFO] Error: A gain function needs to be passed at compile time via the environment variable `GAIN` -- Please recompile...
```

The explanation is that, since we experimented with three different gain functions, the desired gain function must be passed in
at compile time via an environment variable. The valid options are `default`, `approx_1`, or `approx_2`. So, recompile as such:
```
GAIN=approx_1 cargo build --release
```

## Running the code
You will need, at bare minimum, a CIFF index corresponding to whatever data you wish to reorder. Some pre-generated CIFF files
can be found [here](https://github.com/osirrc/ciff#reference-lucene-indexes).

For our following example, let's grab the Robust04 CIFF file.
```
mkdir example
cd example
wget https://www.dropbox.com/s/rph6udiqs2k7bfo/robust04-complete-20200306.ciff.gz?dl=0 -O robust-ciff.gz
gunzip robust-ciff.gz
mv robust-ciff robust.ciff
```

Now, let's run the BP algorithm, output a reordered CIFF file, and compute the loggap improvement.
```
../target/release/create-rgb --input robust.ciff --output-ciff robust-reordered.ciff --loggap 
03:28:02 [INFO] Using the `approx_1` gain function.
03:28:02 [INFO] Opt { input: "robust.ciff", output_ciff: Some("robust-reordered.ciff"), min_len: 4096, cutoff_frequency: 0.1, recursion_stop: 16, swap_iterations: 20, loggap: true, sort_leaf: false, max_depth: 100, input_fidx: None, output_fidx: None, output_mapping: None }
03:28:02 [INFO] (1) building forward index
create_fwd: ⠒ [00:00:02] [████████████████████░░░░░░░░░░░░░░░░░░░░] (482000/923436, ETA 2s, SPEED: 205002/s)                                                                                crecreate_fwd: ⠙ [00:00:06] [████████████████████████████████░░░░░░░░] (748000/923436, ETA 2s, SPEED: 122324/s)                                                                             03:28:10 [INFO] forward index stats:
03:28:10 [INFO] 	total terms: 923436
03:28:10 [INFO] 	discarded frequent terms: 384
03:28:10 [INFO] 	discarded infrequent terms: 920126
03:28:10 [INFO] 	remaining terms: 2927
03:28:10 [INFO] (2) sort empty docs to the back
03:28:11 [INFO] fwd duration: 8.87 secs
03:28:11 [INFO] docs 528030 non_empty 527908
03:28:11 [INFO] put docs back into default order...
03:28:11 [INFO] (3) perform graph bisection
03:28:20 [INFO] rgb duration: 9.57 secs
03:28:20 [INFO] (4) clear forward index
03:28:21 [INFO] (5) starting output operations...
03:28:21 [INFO]  --> (5.2) write new ciff file
03:28:21 [INFO] writing to ciff file: robust-reordered.ciff
03:28:39 [INFO] write duration: 18.80 secs
03:28:39 [INFO] (6) compute loggap cost
03:28:44 [INFO] 	before reorder: 3.975 BPI
03:28:49 [INFO] 	 after reorder: 2.968 BPI
03:28:49 [INFO] ALL DONE! duration: 47.67 secs
```

So, with this configuration:
 - `approx_1` gain function,
 - minimum postings length of 4096, 
 - maximum postings length of 0.1 * N (where N is the number of documents in the collection), 
 - 20 iterations per level, and 
 - the recursion depth fixed by only recursing while there are more than 16 elements within each partition,
we observe the RGB process taking about 10 seconds, improving loggap from 3.975 to 2.968.

Running the same configuration with the `default` gain function takes about 20 seconds, and yields a final 
loggap of 2.989. Similarly, using `approx_2` takes 6 seconds, and yields a loggap of 3.038.

## Settings and Configuration

A full suite of settings can be found using the `--help` flag, and are listed as follows:
```
create-rgb 0.1.0
Reorders documents using recursive graph bisection and ciff files.

USAGE:
    create-rgb [FLAGS] [OPTIONS] --input <input>

FLAGS:
    -h, --help         Prints help information
    -l, --loggap       Show loggap cost
        --sort-leaf    Sort leaf by identifier
    -V, --version      Prints version information

OPTIONS:
    -c, --cutoff-frequency <cutoff-frequency>
            Maximum length to consider in percentage of documents in the index [default: 0.1]

    -i, --input <input>                          Input file ciff file
        --input-fidx <input-fidx>                Read forward index
        --max-depth <max-depth>                  Maximum depth [default: 100]
    -m, --min-len <min-len>                      Minimum number of occurrences to consider [default: 4096]
    -o, --output-ciff <output-ciff>              Output ciff file
        --output-fidx <output-fidx>              Output forward index
        --output-mapping <output-mapping>        Dump the document map
    -r, --recursion-stop <recursion-stop>        Min partition size [default: 16]
    -s, --swap-iterations <swap-iterations>      Swap iterations [default: 20]
```

For example, you can save a forward index using the `--output-fidx` command, and can read a saved forward index
with the `--input-fidx` flag. If you only wish to dump the reordered document map, use the `--output-mapping`
flag. 

##
