# iGenVar

[![Build Status](https://github.com/seqan/iGenVar/workflows/iGenVar%20CI/badge.svg)](https://github.com/seqan/iGenVar/actions?query=workflow%3A%22iGenVar+CI%22+branch%3Amaster)

The official repository for the iGenVar project.

## Installation

Instructions:
1. clone this repository: `git clone --recurse-submodules https://github.com/seqan/iGenVar.git`
    or `git clone https://github.com/seqan/iGenVar.git` and fetch the seqan3 submodule after cloning: `git submodule update --recursive --init`
2. create a build directory and visit it: `mkdir build && cd build`
3. run cmake: `cmake ../iGenVar`
4. build the application: `make`
5. optional: build and run the tests: `make test` or `ctest`
6. optional: build the api documentation: `make doc`
7. execute the app: `./bin/iGenVar`

(Built using the [SeqAn3 App Template](https://github.com/seqan/app-template))

We created a small example, which you can use to test run our app:
`./bin/iGenVar -j ./test/data/mini_example.sam -o ./test/data/output.vcf -m cigar_string -m split_read -l 5`
