// Link libgomp explicitly so we can call `omp_set_num_threads` directly.
// OpenBLAS depends on libgomp transitively, but its symbols are not
// exposed on the link line unless we ask for them.
fn main() {
    println!("cargo:rustc-link-lib=dylib=gomp");
}
