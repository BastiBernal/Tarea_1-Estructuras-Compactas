add_executable(bench_example example.cpp)
target_link_libraries(bench_example PRIVATE my-project bench-lib)
target_compile_options(bench_example PRIVATE -O3)

add_executable(patron_fijo_benchmark patron_fijo_benchmark.cpp)
target_link_libraries(patron_fijo_benchmark PRIVATE my-project bench-lib sdsl divsufsort divsufsort64)
target_compile_options(patron_fijo_benchmark PRIVATE -O3)