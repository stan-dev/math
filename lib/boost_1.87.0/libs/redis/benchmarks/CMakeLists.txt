
add_library(benchmarks_options INTERFACE)
target_link_libraries(benchmarks_options INTERFACE boost_redis_src)
target_link_libraries(benchmarks_options INTERFACE boost_redis_project_options)
target_compile_features(benchmarks_options INTERFACE cxx_std_20)

add_executable(echo_server_client cpp/asio/echo_server_client.cpp)
target_link_libraries(echo_server_client PRIVATE benchmarks_options)

add_executable(echo_server_direct cpp/asio/echo_server_direct.cpp)
target_link_libraries(echo_server_direct PRIVATE benchmarks_options)

# TODO
#=======================================================================

#.PHONY: bench
#bench:
#	pdflatex --jobname=echo-f0 benchmarks/benchmarks.tex
#	pdflatex --jobname=echo-f1 benchmarks/benchmarks.tex
# pdftoppm {input.pdf} {output.file} -png