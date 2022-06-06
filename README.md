# IPUMA Library

Multi-purpose alignment library on the Graphcore Intelligent Processing Unit (IPU).

This project contains the `ipuma-library` which implements implementations of the Smith-Waterman alignment algorithm for the IPU.

## Installation

### Requirements

The following software are required for building and running the alignment
libraries:

```
poplar-sdk >=2.3.0
cmake >= 3.18.4
gcc >= 7.5.0
```

#### Python script dependencies

For some python scripts used in extracting log files numpy and pandas are
required. These are not required for running the alignment binaries themselves!

Preferably execute the following commands in a venv or a conda environment:

```
$ pip install numpy==1.22.2 pandas==1.4.1
```

### Building

```
$ mkdir -p build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release -GNinja ..
$ ninja
```

## Testing

Correctness tests are contained in `ipuma-tests`.

```
$ build/test/ipuma-tests
```

Artificial performance test scenarios on generated data are contained in `ipuma-perf`. These do not fail, but should present GCUPS results that can be compared.

```
$ build/test/ipuma-perf
```

## Running

The `ipusw` binary can be used to perform benchmark measurements on specified datasets on the IPU.

The `cpusw` binary is a wrapper around [StripedSmithWaterman](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library).

### Running paper benchmarks

The `./scripts/run_ipusw_benchmarks_mult.sh` script will execute IPU-based
alignment benchmarks.

The script expects alignment benchmark datasets to be located in `./download`
and will output log files to `./output`.

The script optionally accepts a name for the execution, otherwise `ipuswrun`
will be used:

```
$ ./scripts/run_ipusw_benchmarks_mult.sh <NAME>
```
