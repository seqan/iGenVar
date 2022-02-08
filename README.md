# iGenVar

[![Build Status](https://github.com/seqan/iGenVar/workflows/iGenVar%20CI/badge.svg)](https://github.com/seqan/iGenVar/actions?query=workflow%3A%22iGenVar+CI%22+branch%3Amaster)

The official repository for the iGenVar project.

iGenVar is intended to be a caller for all types of genetic variation: SNPs, indels and larger structural variations
(insertions, deletions, inversions, translocations, CNVs, nested SVs).
It uses both Illumina short reads and PacBio long reads for this purpose.

David Heller in the Vingron lab of the MPI-MG and Tim White in the Kehr lab at BIH have both developed an SV caller for
long read sequencing data. Instead of competing with each other, we want to join forces and combine the two tools, SVIM
and SVIRL, into one better and more versatile tool.
On the other hand, there were some tool developments in the Reinert Lab (FU Berlin): Vaquita, a short read SV caller,
SViper, a refinement tool and Vaquita-LR a further development of Vaquita for long reads.
We want to combine these approaches and use the SeqAn3 library as a basis for this new tool.

## Current status:

We can call insertions and deletions from long read data (SVIM methods implemented).
For more information, see the release plan at the bottom of the page.

## Installation

Instructions:

1. clone this repository: `git clone --recurse-submodules https://github.com/seqan/iGenVar.git`
    or `git clone https://github.com/seqan/iGenVar.git` and fetch the seqan3 submodule after cloning:
    `git submodule update --recursive --init`
2. create a build directory and visit it: `mkdir build && cd build`
3. run cmake: `cmake ../iGenVar`
4. build the application: `make`
5. optional: build and run the tests: `make test` or `ctest`
6. optional: build the api documentation: `make doc`
7. execute the app: `./bin/iGenVar`

(Built using the [SeqAn3 App Template](https://github.com/seqan/app-template))

We created small examples, which you can use to test our app:
```bash
./bin/iGenVar -i ./test/data/paired_end_mini_example.sam -j ./test/data/single_end_mini_example.sam \
-o ./test/data/output.vcf --method cigar_string --method split_read --min_var_length 5
```

## Release plan:

<p align="center"><img height="500" src="https://github.com/seqan/iGenVar/blob/863297c128d9fa67a4ab51206d7338dcbdd8ca1b/doc/ReleasePlan.png"></p>

## Dependencies
We support gcc/g++ 9 and higher.
