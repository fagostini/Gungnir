# Gungnir
Approximate string matching (Hamming distance) for genomic sequences, running on GPU

## Requirements
* NVIDIA CUDA-capable GPU (the list of capable graphic cards is available [here](https://developer.nvidia.com/cuda-gpus))
* NVIDIA CUDA Toolkit (installation instructions can be found [here](https://docs.nvidia.com/cuda/index.html#installation-guides))

## Clone and compile

Clone the repository:

```bash
git clone git@github.com:fagostini/Gungnir.git
```

Compile the source code:

```bash
cd Gungnir
mkdir bin
make
```

> **_Troubleshoot_**
> 
> In case the compilation fails or there is the need for specific parameters, it might be necessary to modify some hardcoded values in the following files.
>
>_[Makefile](https://github.com/fagostini/Gungnir/blob/master/Makefile)_
>* The CUDA root directory is set to the default installation path (_i.e._, '/usr/local/cuda');
>* The default compute capability is 6.1 (_i.e._, sm_61).
>
>_[main.cu](https://github.com/fagostini/Gungnir/blob/master/main.cu)_
>* The maximum length for the input sequences (MAXLEN) and number of mismatches (MISMATCH) are set to 100 and 9, respectively;
>* The maximum total reference length (MAXREF) is set to 3 billion characters;
>* The block size (BLOCKSIZE) is set to 32 in order to maximise the parallel threads.

## Usage

The previous step will generate an executable file, namely _runHamming_, that can be run using the following synthax:

```bash
./runHamming <query FASTA> <reference FASTA>
```

## Test [optional]

Run the main test:

```bash
make test
```

The output on the console should be 'Test PASSED', while the output file _results.tsv_ should contain the following lines:

```
>Reference  TCCAGCGCCCGAGCCGTCCAGGCGGCCAGCAGGAGCAGTG  1  1  1  1  1  0  0  0  0  0
>HamDist_1  TCCAGCGCGCGAGCCGTCCAGGCGGCCAGCAGGAGCAGTG  1  2  1  1  0  0  0  0  0  0
>HamDist_2  TCCAGCGCGCGAGCCGTCCAGGCGCCCAGCAGGAGCAGTG  1  2  2  0  0  0  0  0  0  0
>HamDist_3  TCCAGCGCGCGAGCCGTCCAGGCGCCCAGCTGGAGCAGTG  1  2  1  1  0  0  0  0  0  0
>HamDist_4  TCCAGCGTGCGAGCCGTCCAGGCGCCCAGCTGGAGCAGTG  1  1  1  1  1  0  0  0  0  0
```

> **_Troubleshoot_**
> 
> If the _results.tsv_ file contains only zeros:
> * Decrease the block size to either 16 or 8
> * Re-compile 
> * Re-run the test.


Run the 1'000 vs 10'000 sequences test:
```
make test10k
```

Run the 10'000 vs 100'000 sequences test:
```
make test100k
```

Run all tests:
```
make testall
```

Profile the execution using _nvprof_ (included in CUDA Toolkit >= 10.0):
```
make profile
```

The command produces some basic profiling results on the console and saves a more detailed version in the '_extra/runHamming-analysis.nvprof_' file, which can be imported into either _nvprof_ or the NVIDIA Visual Profiler (_nvvp_).

## Bugs and Issues

Report bugs as issues on the [GitHub repository](https://github.com/fagostini/Gungnir/issues)

## Author

* [Federico Agostini](https://github.com/fagostini)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
