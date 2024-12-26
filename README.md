# Scalable Zero-knowledge Proofs for Non-linear Functions in Machine Learning
Source code and the implementation of the USENIX Security'24 paper - [_Scalable Zero-knowledge Proofs for Non-linear Functions in Machine Learning_](https://www.usenix.org/conference/usenixsecurity24/presentation/hao-meng-scalable).

## Introduction
A novel scalable ZK proof framework for non-linear mathematical functions in machine learning based on LUT techniques.

## Contents
This repository consists of the following parts:
- __src__: Implementations for our ZKP protocols.
- __test__: Test script for performance evaluation.

## Installation
***NOTE: Our implementation is built based on [emp-zk](https://github.com/emp-toolkit/emp-zk) and a complete compilation process is required.
For any issues occurred during step 1, please refer to [emp-zk troubleshooting](https://github.com/emp-toolkit/emp-zk/issues).***

1. Install emp-tool, emp-ot and emp-zk.

```bash
python install.py --deps --tool --ot --zk
```

2. Update code for zkMath.

```bash
cp -rf src/* emp-zk/emp-zk/ && rm -rf src \
&& cp -rf test/* emp-zk/test/ && rm -rf test \
&& mv -f CMakeLists.txt emp-zk/
```

3. Recompile the code, then the binary executables can be found at ./bin

```bash
cd emp-zk
cmake .
make -j8
sudo make install
```

## Usage
### Testing on localhost

`./run ./bin/[binary]`

### Testing on two

1. Change the IP address in the test code, or use `tmux` command to initialize another terminal

2. run `./bin/[binary] 1 [port]` on one machine and

   run `./bin/[binary] 2 [port]` on the other machine.

## Disclaimer
This repository is a proof-of-concept prototype.
