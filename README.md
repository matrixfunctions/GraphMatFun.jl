# Computation graphs for matrix functions

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://matrixfunctions.github.io/GraphMatFun.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://matrixfunctions.github.io/GraphMatFun.jl/dev)
[![Build Status](https://travis-ci.com/matrixfunctions/GraphMatFun.jl.svg?token=1D6L5XDzryqshBZ5ijUH&branch=main)](https://travis-ci.com/matrixfunctions/GraphMatFun.jl)
[![codecov](https://codecov.io/gh/matrixfunctions/GraphMatFun.jl/branch/main/graph/badge.svg?token=ZTKNBNMDEZ)](https://codecov.io/gh/matrixfunctions/GraphMatFun.jl)


This package contains functionality to represent, manipulate and optimize algorithms for matrix functions using directed acyclic graphs. It also contains code generate features for other languages (Julia, MATLAB, C/BLAS).

* Features and examples are described in the manuscript [Computation graph for matrix functions](https://arxiv.org/abs/2107.12198).
* Functions and usage is described in the [online package documentation](https://matrixfunctions.github.io/GraphMatFun.jl/dev/).
* Data files for various matrix functions are available in [the package GraphMatFunData](https://github.com/matrixfunctions/GraphMatFunData).


## Installation

The package is registered with the Julia central registry and can be installed with the command:

```
julia> ]
(v1.7) pkg> add GraphMatFun
```

You can now follow the examples in the [online package documentation](https://matrixfunctions.github.io/GraphMatFun.jl/dev/).