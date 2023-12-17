# SuStain.jl

[Original SuStain Paper](https://doi.org/10.1038/s41467-018-05892-0)

## Installation

Run the following commands on a Unix-based system with Julia and make preinstalled:

```shell
git clone https://github.com/neuro-stivenr/SuStain.jl
cd SuStain.jl
make
```

This should generate a SuStain.sh file in the repo directory.
You can move this file anywhere on your system, and it will work, as long as you don't move the repo directory.
If you move the repo directory, re-run make.

## Usage

```shell
./SuStain.sh [nthreads] [datapath] [var1,var2,var3] [n_stages] [n_runs] [outpath]
```
