cmake_minimum_required(VERSION 4.1)
project(bioinf)

set(CMAKE_CXX_STANDARD 20)

add_executable(bioinf main.cpp)

configure_file(reads.fastq reads.fastq COPYONLY)