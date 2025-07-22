# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

libpapilo is a fork of scipopt/papilo that aims to provide PaPILO (Parallel Presolve for Integer and Linear Optimization) as a **shared library with a C API**. While the original PaPILO is implemented in C++, this fork's goal is to make the functionality accessible from other programming languages through a stable C interface.

**Current Status**: This is a fresh fork with no development yet. The main objectives are:
- Create a shared library (.so/.dylib/.dll) from the existing C++ codebase
- Design and implement a C API that exposes key PaPILO functionality
- Maintain compatibility with the original PaPILO while providing a stable ABI

**Important**: This fork focuses solely on the presolving functionality. Solver integrations (SCIP, Gurobi, HiGHS, etc.) are not needed and should be excluded from the shared library.

The original PaPILO is a C++14-based presolve package for (mixed integer) linear programming problems with support for parallel execution and multiple precision arithmetic, licensed under LGPLv3.

## Build Commands

```bash
# Standard build
mkdir build
cd build
cmake ..
make -j

# With solver integration (e.g., SCIP)
cmake -DSCIP_DIR=PATH_TO_SCIP_BUILD_DIR ..
make -j

# Run unit tests
make test
# or
ctest

# Run a specific test
./test/unit_test "test-name"

# Install
sudo make install
```

## Key Dependencies

- **C++ Standard**: C++14 (required)
- **CMake**: >= 3.11.0
- **Intel TBB**: >= 2020 (for parallelization, strongly recommended)
- **Boost**: >= 1.65 (headers required; iostreams, program_options, serialization for binaries)

## Source Directory Structure

```
src/
├── papilo.cpp              # Main executable entry point
├── papilolib.cpp/.h        # Library interface (initial C API attempt)
├── convMPS.cpp             # MPS file conversion utility
├── duplicates.cpp          # Duplicate detection utility
└── papilo/
    ├── core/               # Core presolving engine
    │   ├── ConstraintMatrix.hpp      # Sparse matrix representation
    │   ├── MatrixBuffer.hpp          # Matrix modification buffer
    │   ├── Presolve.hpp              # Main presolve orchestrator
    │   ├── Problem.hpp               # Problem data structure
    │   ├── ProblemBuilder.hpp        # Problem construction
    │   ├── Reductions.hpp            # Reduction data structures
    │   ├── RowFlags.hpp/ColFlags.hpp # Row/column properties
    │   ├── SingleRow.hpp             # Single row analysis
    │   ├── SparseStorage.hpp         # Sparse data storage
    │   └── postsolve/               # Postsolve components
    │       ├── Postsolve.hpp         # Main postsolve class
    │       ├── PostsolveStorage.hpp  # Storage for postsolve info
    │       └── ReductionType.hpp     # Reduction type definitions
    ├── presolvers/             # 17 presolve methods
    │   ├── CoefficientStrengthening.hpp
    │   ├── ConstraintPropagation.hpp
    │   ├── DominatedCols.hpp
    │   ├── DualFix.hpp
    │   ├── DualInfer.hpp
    │   ├── FixContinuous.hpp
    │   ├── FreeVarSubstitution.hpp
    │   ├── ImplIntDetection.hpp
    │   ├── ParallelColDetection.hpp
    │   ├── ParallelRowDetection.hpp
    │   ├── Probing.hpp
    │   ├── SimpleProbing.hpp
    │   ├── SimpleSubstitution.hpp
    │   ├── SimplifyInequalities.hpp
    │   ├── SingletonCols.hpp
    │   ├── SingletonStuffing.hpp
    │   └── Sparsify.hpp
    ├── interfaces/            # Solver integrations
    │   ├── ScipInterface.hpp
    │   ├── SoplexInterface.hpp
    │   ├── GurobiInterface.hpp
    │   ├── HighsInterface.hpp
    │   └── GlopInterface.hpp
    ├── io/                    # File I/O
    │   ├── MpsParser.hpp      # MPS format parser
    │   ├── MpsWriter.hpp      # MPS format writer
    │   ├── OpbParser.hpp      # OPB format parser
    │   └── SolParser.hpp      # Solution file parser
    ├── misc/                  # Utilities
    │   ├── compress_vector.hpp
    │   ├── fmt.hpp
    │   ├── Hash.hpp
    │   ├── MultiPrecision.hpp
    │   ├── Num.hpp
    │   ├── ParameterSet.hpp
    │   ├── PrimalDualSolValidation.hpp
    │   ├── StableSum.hpp
    │   ├── Timer.hpp
    │   ├── Validation.hpp
    │   └── Vec.hpp
    ├── external/              # Bundled third-party libraries
    │   ├── catch/            # Unit testing framework
    │   ├── fmt/              # Formatting library
    │   ├── lusol/            # LU decomposition (Fortran)
    │   │   └── src/
    │   │       ├── lusol.f90
    │   │       ├── lusol6b.f
    │   │       ├── lusol7b.f
    │   │       ├── lusol8b.f
    │   │       └── lusol_util.f
    │   ├── pdqsort/          # Pattern-defeating quicksort
    │   └── ska/              # Hash map implementations
    └── verification/          # Certificate generation
        └── VeriPB.hpp

## High-Level Architecture

The codebase is template-based to support different numerical precision types (double, quadruple, rational):

1. **Core Engine** (`src/papilo/core/`):
   - `Problem<REAL>`: Problem representation with constraint matrix, bounds, objectives
   - `Presolve<REAL>`: Orchestrates presolving routines and manages the presolve loop
   - `Postsolve<REAL>`: Transforms solutions back to original problem space
   - `postsolve/`: Contains postsolve storage and reduction classes

2. **Presolvers** (`src/papilo/presolvers/`): 17 individual presolve methods that detect and apply reductions:
   - Each presolver inherits from `PresolveMethod<REAL>`
   - Key methods: `execute()`, `initialize()`, `updateResults()`
   - Examples: DualFix, ConstraintPropagation, ImplIntDetection, ParallelRowDetection

3. **Solver Interfaces** (`src/papilo/interfaces/`):
   - Integrations with SCIP, SoPlex, Gurobi, HiGHS, GLOP
   - **Note**: These interfaces are NOT used in this fork - focus is purely on presolving

4. **I/O** (`src/papilo/io/`):
   - MPS format parser/writer
   - OPB (pseudo-boolean) format support
   - Solution file handling

5. **External Dependencies** (`src/papilo/external/`):
   - **LUSOL**: Fortran-based sparse LU factorization solver for linear systems
   - **fmt**: Modern C++ formatting library
   - **Catch**: Unit testing framework
   - **pdqsort**: High-performance sorting algorithm
   - **ska**: Efficient hash map implementations

## Testing Framework

1. **Unit Tests** (`test/`):
   - Uses Catch2 framework
   - Run with: `make test` or `./test/unit_test`
   - Tests cover individual presolvers and core functionality

2. **Performance Testing** (`check/`):
   - Script-based testing framework for instance solving
   - Run locally: `make test` (from root)
   - Cluster support via `check_cluster.sh`

## Common Development Tasks

```bash
# Build and run tests after changes
cd build
make -j && make test

# Run papilo executable
./papilo solve -f instance.mps
./papilo presolve -f instance.mps -v postsolve.info
./papilo postsolve -v postsolve.info -u solution.sol

# Debug a specific presolver
./test/unit_test "dual-fix-happy-path"
```

## Important Implementation Notes

- The presolve loop in `Presolve<REAL>::presolve()` applies presolvers iteratively until no more reductions are found
- Numerical tolerances are managed through `Tolerances<REAL>` class
- Postsolve information is stored in reverse order (LIFO) for correct reconstruction
- Parallel execution uses Intel TBB's parallel algorithms when available
- The `Message` class provides configurable logging levels for debugging

## C API Development Guidelines (To Be Implemented)

When developing the C API for libpapilo, consider:

1. **Handle-based API**: Use opaque pointers to hide C++ implementation details
2. **Error handling**: Return error codes instead of throwing exceptions
3. **Memory management**: Provide explicit create/destroy functions
4. **Numerical types**: Initially focus on double precision, consider multiple precision support later
5. **Key functionality to expose**:
   - Problem creation and modification (constraint matrix, bounds, objectives)
   - Presolve execution with parameter control
   - Retrieving the presolved problem
   - Postsolve functionality to transform solutions back
   - Status and statistics queries
   - Parameter configuration

**Excluded functionality**:
- Solver interfaces (SCIP, Gurobi, etc.) - users will solve with their own choice of solver
- The main executable functionality - this fork provides only the library