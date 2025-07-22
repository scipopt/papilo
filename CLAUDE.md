# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

libpapilo is a fork of scipopt/papilo that aims to provide PaPILO (Parallel Presolve for Integer and Linear Optimization) as a **shared library with a C API**. While the original PaPILO is implemented in C++, this fork's goal is to make the presolving functionality accessible from other programming languages through a stable C interface.

**Key Design Decisions:**
- **Presolving only**: No solver integrations (SCIP, Gurobi, etc.) - users solve with their own choice of solver
- **Double precision only**: Simplifies C API by avoiding C++ template complexity (original supports double/quadruple/rational)
- **Fresh fork**: No development yet - starting from clean slate

The original PaPILO is a C++14-based presolve package for (mixed integer) linear programming problems with support for parallel execution and multiple precision arithmetic, licensed under LGPLv3.

## Development Roadmap

### Phase 1: Environment Setup
- [ ] Remove solver integration dependencies from build system
- [ ] Configure CMake for shared library generation (.so/.dylib/.dll)
- [ ] Set up testing framework for C API validation

### Phase 2: Core C API Design
- [ ] Extend existing `papilolib.h` with presolve-only functions
- [ ] Implement problem construction API (building on existing foundation)
- [ ] Add presolve execution API (individual presolvers + full pipeline)
- [ ] Add presolved problem data extraction API
- [ ] Add manual postsolve API

### Phase 3: Implementation & Testing
- [ ] Implement C++ wrapper layer for double-precision instantiations
- [ ] Create comprehensive test suite matching `test/papilo/presolve/` patterns
- [ ] Add memory management and error handling
- [ ] Performance testing and optimization

### Target C API Interface

The goal is to support the three-stage workflow demonstrated in PaPILO tests:

1. **Construction**: `papilo_problem_create()` → set dimensions → add variables/constraints → add matrix entries
2. **Execution**: `papilo_solver_presolve_only()` → query status and reductions
3. **Validation**: Extract presolved problem → manual postsolve when needed

## Build Commands

```bash
# Standard build
mkdir build && cd build
cmake .. && make -j

# Run tests
make test

# Run specific test
./test/unit_test "dual-fix-happy-path"
```

## Key Dependencies

- **C++ Standard**: C++14
- **CMake**: >= 3.11.0
- **Intel TBB**: >= 2020 (for parallelization)
- **Boost**: >= 1.65 (headers only for this fork)

## Architecture Overview

### Core Components (C++ Implementation)
- **Problem/ProblemBuilder**: Problem representation and construction
- **Presolve**: Main orchestrator with 17 individual presolving methods
- **Reductions**: Tracks all presolving modifications
- **PostsolveStorage**: Enables solution transformation back to original space

### Template Design
Original PaPILO uses `REAL` template parameter for numerical precision. This fork simplifies by using `double` exclusively, avoiding template complexity in the C API.

### Existing C API Foundation
`src/papilolib.h/cpp` provides excellent groundwork:
- ✅ Problem construction functions
- ✅ Parameter management system  
- ❌ Only integrated solving (presolve+solve+postsolve)
- ❌ No presolve-only execution
- ❌ No presolved problem data access

## Detailed Implementation Notes

### Source Structure
```
src/papilo/
├── core/              # Problem, Presolve, PostsolveStorage
├── presolvers/        # 17 presolving methods
├── io/                # MPS/OPB parsers (needed for C API)
├── misc/              # Vec<T>, utilities, numerical types
└── external/          # Third-party libraries (fmt, LUSOL, etc.)
```

### Key Implementation Details
- **Vec<T>**: Custom vector type (currently identical to std::vector, designed for future allocator optimization)
- **LUSOL**: Fortran-based sparse LU solver for advanced presolving techniques
- **Explicit instantiations**: All templates pre-instantiated for double/quad/rational types
- **LIFO postsolve**: Reductions stored in reverse order for correct solution reconstruction

### Test-Driven API Design
Analysis of `test/papilo/presolve/` reveals the essential API patterns:

**Problem Construction Pattern:**
```cpp
ProblemBuilder<double> pb;
pb.reserve(nnz, nrows, ncols);
pb.setColUbAll(bounds) / setObjAll(coeffs) / addEntryAll(triplets);
Problem<double> problem = pb.build();
```

**Presolve Execution Pattern:**
```cpp
PresolveOptions options{};
ProblemUpdate<double> update(problem, postsolve, stats, options, num, msg);
Reductions<double> reductions{};
PresolveStatus status = presolver.execute(problem, update, num, reductions, timer);
```

**Result Validation Pattern:**
```cpp
// Check status: kUnchanged, kReduced, kInfeasible, kUnbounded
// Access reductions: count, individual details (col, row, value, type)
// Extract modified problem data
```

The C API should mirror these patterns while providing C-compatible data structures and error handling.

## Development Guidelines

- **Memory safety**: All C API functions must handle allocation failures gracefully
- **Error handling**: Use return codes, never throw exceptions across C boundary  
- **Opaque handles**: Hide C++ implementation details behind void* handles
- **Resource management**: Provide explicit create/destroy functions for all objects
- **Thread safety**: Document threading requirements (likely requires external synchronization)

## Next Steps

1. Examine `papilolib.cpp` implementation details
2. Design presolve-only extensions to existing API
3. Create proof-of-concept for problem construction → presolve → data extraction workflow
4. Implement comprehensive test coverage matching existing unit tests