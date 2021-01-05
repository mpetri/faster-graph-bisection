# incremental-index-rgb

```
(base) ➜  incremental-index-rgb git:(main) ✗ ./target/release/create-rgb --help                                                                               
create-rgb 0.1.0
Reorders documents using recursive graph bisection and ciff files.

USAGE:
    create-rgb [FLAGS] [OPTIONS] --input <input> --output <output>

FLAGS:
    -h, --help       Prints help information
    -l, --loggap     Show loggap cost
    -V, --version    Prints version information

OPTIONS:
    -c, --cutoff-frequency <cutoff-frequency>
            Maximum length to consider in percentage of documents in the index [default: 0.1]

    -i, --input <input>                          Input file ciff file
    -m, --min-len <min-len>                      Minimum number of occurrences to consider [default: 4096]
    -o, --output <output>                        Output ciff file
    -r, --recursion-stop <recursion-stop>        Min partition size [default: 16]
    -s, --swap-iterations <swap-iterations>      Swap iterations [default: 20]
(base) ➜  incremental-index-rgb git:(main) ✗ ./target/release/create-rgb -i ./test.ciff -o ./bla.ciff -l -m 10
00:58:29 [INFO] Opt { input: "./test.ciff", output: "./bla.ciff", min_len: 10, cutoff_frequency: 0.1, recursion_stop: 16, swap_iterations: 20, loggap: true }
00:58:29 [INFO] (1) create forward index from ciff
00:58:29 [INFO] forward index stats:
00:58:29 [INFO]         total terms: 564915
00:58:29 [INFO]         discarded frequent terms: 569
00:58:29 [INFO]         discarded infrequent terms: 538259
00:58:29 [INFO]         remaining terms: 26088
00:58:29 [INFO] (2) sort empty docs to the back
00:58:29 [INFO] fwd duration: 0.64 secs
00:58:29 [INFO] docs 18661 non_empty 18661
00:58:29 [INFO] (3) perform graph bisection
00:58:35 [INFO] rgb duration: 5.33 secs
00:58:35 [INFO] (4) clear forward index
00:58:35 [INFO] (5) write new ciff file
00:58:35 [INFO] writing to ciff file: ./bla.ciff
00:58:37 [INFO] write duration: 2.16 secs
00:58:37 [INFO] (6) compute loggap cost
00:58:37 [INFO]         before reorder: 3.681 BPI
00:58:38 [INFO]          after reorder: 3.717 BPI
00:58:38 [INFO] ALL DONE! duration: 8.30 secs
```