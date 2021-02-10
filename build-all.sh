for g in default approx_1 approx_2; do
    GAIN="$g" cargo build --release
    cp target/release/create-rgb target/release/create-rgb-$g
done 
cp target/release/create-rgb-default target/release/create-rgb
