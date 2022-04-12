# IPUMA Library

Multi-purpose alignment library on the Graphcore Intelligent Processing Unit (IPU).

This project contains the `ipuma-library` which implements implementations of the Smith-Waterman alignment algorithm for the IPU.

## Installation

### Requirements

poplar-sdk, cmake, gcc

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