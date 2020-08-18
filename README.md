# PaPILO --- Parallel Presolve for Integer and Linear Optimization

PaPILO, a C++14 based software package, provides parallel presolve routines for (mixed integer) linear programming problems. The routines
are implemented using templates which allows switching to higher precision or rational arithmetic using the boost multiprecision package.

Additionally to the distribution here in github under the LGPLv3, PaPILO is also distributed as part of the
SCIP Optimization Suite which is available under https://scipopt.org/.

PaPILO can be used as a header-based library and also provides an executable.
Using the executable it is possible to presolve and postsolve MILP instances based on files.
Additionaly PaPILO can be linked to SCIP, SoPlex, and HiGHS (https://github.com/ERGO-Code/HiGHS) solvers and act as a frontend. In this setting PaPILO passes the presolved problem to those solvers and applies the postsolve step to the optimal solution.

When PaPILO is compiled as part of the SCIP Optimization Suite linking of SoPlex and SCIP solvers is performed automatically.

# Dependencies

External dependencies that need to be installed by the user are the Intel TBB runtime library and boost headers.
The executable additional requires some of the boost runtime libraries that are not required when papilo is used as
a library.
Under the folder external/ there are additional packages that are directly included within PaPILO and have a
liberal open-source license.

Intel TBB is also included and PaPILO tries to compile a static version of TBB if the runtime library is missing.
This fails on some systems currently and if any problems occur it is recommended to install an Intel TBB runtime library
from, e.g., the systems package manager (`libtbb2` on debian/ubuntu based distributions).

# Building

Building PaPILO works with the standard cmake workflow:
```
mkdir build
cd build
cmake ..
make
```

Solvers that are found in the system are automatically linked to the executable.
Additionally one can speciy the locations of solvers, e.g. with -DSCIP_DIR=<location of scip-config.cmake>, to allow
PaPILO to find them in non-standard locations.

# Usage of the binary

The PaPILO binary provides a list of all available functionality when the help flag `-h` or `--help` is specified.
The binary provides the three subcommands `solve`, `presolve`, and `postsolve`. If no solvers are linked the `solve` subcommand will fail and print that as a message.

Next we provide a small example of how the binary can be used to apply presolving and postsolving to problem based on files.

Assuming a problem instance is stored in the file `problem.mps` the following call will apply presolving with standard settings and write the reduced problem to `reduced.mps` and all information that is needed for postsolve to the binary archive `reduced.postsolve`.
```
bin/papilo presolve -f problem.mps -r reduced.mps -v reduced.postsolve
```

Now we can use the reduced problem `reduced.mps` to obtain a solution
using any solver or from any other source to the file `reduced.sol`.
The format of the solution should be one value per line given like this:
```
<variable name>        <value>
```
This is compatible with the solutions given on the MIPLIB 2017 homepage https://miplib.zib.de and with the solutions written by the SCIP solver.
Variable names that are not found in reduced problem are ignored.

The command for applying the postsolve step to the solution `reduced.sol` is then
```
bin/papilo postsolve -v reduced.postsolve -u reduced.sol -l problem.sol
```
Giving the parameter `-l problem.sol` is optional and will store the solution transformed to the original space under `problem.sol`.
The output of papilo contains some information about violations and objective value of the solution.

If PaPILO was linked to a suitable solver, then the above can also be achieved by using the `solve` subcommand like this:
```
bin/papilo solve -f problem.mps -l problem.sol
```
This will presolve the problem, pass the reduced problem to a solver, and subsequently transform back the optimal solution returned by the solver and write it to problem.sol.

# Using PaPILO as a library

As PaPILO provides a templated C++ interface the type used for numerical computations must be specified. During configuration time PaPILO scans the system and provides the fastest available numeric types for quadprecision and for exact rational arithmetic in the file
`papilo/misc/MultiPrecision.hpp`. Including this file will introduce the types
```
papilo::Quad
papilo::Float100
papilo::Float500
papilo::Float1000
papilo::Rational
```
The numeric type used by papilo will be reffered to as REAL in the following section. It can be any of above types as well as simply `double` for using standard double precision arithmetic.

To avoid confusion with types a short note on types like `papilo:Vec` and `papilo::String`.
Those types are aliases for types from the standard library, `std::vector` and `std::string`, that possibly use an adjusted allocator. If nothing is altered regarding the allocator the type `papilo::Vec` will be exactly the same as `std::vector`.
It can be changed by adding a partial specialization of `papilo::AllocatorTraits<T>` after including `papilo/misc/Alloc.hpp` but before including any other header of papilo.

The C++ interface for using PaPILO mainly evolves around the classes
`papilo::Presolve<REAL>`, which controls the presolving routines, `papilo::Problem<REAL>` which holds the problem instance, and `papilo::Postsolve<REAL>`, which can transform solutions in the reduced space into solutions for the original problem space. The includes for those classes
under `papilo/core/{Problem,Postsolve,Presolve}.hpp`.

## Creating an instance of `papilo::Problem<REAL>`
The PaPILO binary uses the MPS parsing routine to construct an instance of `papilo::Problem<REAL>` with the call `papilo::MpsParser<REAL>::loadProblem("problem.mps")`.

For feeding a problem to PaPILO programmatically, there is the class
`papilo::ProblemBuilder<REAL>`. This class allows to efficiently build an `papilo::Problem<REAL>` incrementally.
The problem definition that PaPILO supports does not use a row sense, but uses left and right hand sides $l$ and $u$ of rows defined as
$\text{l} \leq a^\top x \leq \text{u}$. For defining a row that is an equation with right hand side $b$ one has to set $l = u = b$. For inequalities either $l$ or $u$ are set to infinity.
One thing where PaPILO differs from many solvers is how infinite values are encoded. Whether for column bounds or left/right hand sides of rows PaPILO encodes the infinite value as a logical condition.
This ensures that regardless of numerical type used for `REAL`, that infinite values are always treated the same.

The member functions for initializing the rows and columns are
```
setNumCols( int ncols )
setNumRows( int nrows )
```
After calling those functions the problem will have no nonzero entries but the given number of columns and rows. The left and right hand side values of rows the rows are set to $0$ as well as the lower and upper bounds of the columns. Additionally all columns are initialized as continuous columns.

Next the following member functions can be used to alter the bound information about rows and columns as well as the objective coefficients and integrality flags of the columns.

*   alter objective coefficient for columns
    ```
    setObj( int col, REAL val )
    ```
*   alter bound information for columns
    ```
    setColLb( int col, REAL lb )
    setColLbInf( int col, bool isInfinite )
    setColUb( int col, REAL ub )
    setColUbInf( int col, bool isInfinite )
    ```
*   mark column to be restricted to integral values or not
    ```
    setColIntegral( int col, bool isIntegral )
    ```
*   alter left and right hand sides of rows
    ```
    setRowLhsInf( int row, bool isInfinite )
    setRowRhsInf( int row, bool isInfinite )
    setRowLhs( int row, REAL lhsval )
    setRowRhs( int row, REAL rhsval )
    ```
*   set names of rows, columns, and the problem
    ```
    setRowName( int row, Str&& name )
    setColName( int col, Str&& name )
    setProblemName( Str&& name )
    ```

Adding nonzeros to the problem can be done with individual nonzero values, row-based, or column-based using the following functions:
```
   /// add the nonzero entries for the given column
   addColEntries( int col, int len, const int* rows, const R* vals )
   /// add a nonzero entry for the given row and column
   addEntry( int row, int col, const REAL& val )
   /// add the nonzero entries for the given row
   addRowEntries( int row, int len, const int* cols, const R* vals )
```
All those functions can be called multiple times, but a nonzero entry for a particular column in a particular row should only be added once.
For maximum efficiency the member function
`papilo::ProblemBuilder<REAL>::reserve(int nnz, int nrows, int ncols)` should be used to reserve all required memory before adding nonzeros.

Finally calling `papilo::ProblemBuilder<REAL>::build()` will return an instance of `papilo::Problem<REAL>` with the information that was given to the builder. The builder can be reused afterwards.

## Presolving an instance of `papilo::Problem<REAL>`

For this section we assume a problem instance is stored in a variable `problem` of type `papilo::Problem<REAL>`.

In order to presolve a problem instance we need to setup an instance of `papilo::Presolve<REAL>` and then call `papilo::Presolve<REAL>::apply(problem)`.

The same instance of `papilo::Presolve<REAL>` can be used for presolving multiple problem instances.
It stores the basic configuration of the presolving routine and constructing it with the default presolvers and settings to presolve the problem is straight forward:
```
papilo::Presolve<REAL> presolve;
presolve.addDefaultPresolvers();
papilo::PresolveResult<REAL> result = presolve.apply(problem);
```

After the above call `result.status` will contain an enum class of the following type:
```
/// result codes of a presolving routine
enum class PresolveStatus : int
{
   /// problem was not changed
   kUnchanged = 0,

   /// problem was reduced
   kReduced = 1,

   /// problem was detected to be unbounded or infeasible
   kUnbndOrInfeas = 2,

   /// problem was detected to be unbounded
   kUnbounded = 3,

   /// problem was detected to be infeasible
   kInfeasible = 4,
};
```

And `result.postsolve` contains an instance of the class `papilo::Postsolve<REAL>`.

## Postsolve of a solution in the reduced problem space

First we construct a `papilo::Solution<REAL>` from a `papilo::Vec<REAL>` of reduced solution values and an empty instance of `papilo::Solution<REAL>` to hold the original space solution.
The interface here is not the simplest for the current functionality. It is like this to support
postsolve of dual solutions in the future. The class `papilo::Solution<REAL>` cannot only hold primal solutions but also dual solution values, even though the postsolve routine does not yet support this step.

Obtaining the original space solution
```
papilo::Vec<REAL> reducedsolvals;
...
// set up the values of the reduced solution in the reduced index space
...
// create reduced solution and original solution
papilo::Solution<REAL> reducedsol(std::move(reducedsolvals));
papilo::Solution<REAL> origsol;

// transform the reduced solution into the original problem space
PostsolveStatus status = result.postsolve.undo(reducedsol, origsol);
```

The value of `status` is `PostsolveStatus::kOk` if everything worked or `PostsolveStatus::kFail` otherwise.
If everything worked then `origsol.primal` contains the primal solution values in the original problem space.

# Algorithmic and implementation details

The release report of the SCIP Optimization Suite 7.0 contains a section about PaPILO. The report is available under http://www.optimization-online.org/DB_HTML/2020/03/7705.html.

Most of the presolve methods implemented in PaPILO are described in the paper "Presolve Reductions in Mixed Integer Programming" by Achterberg et al.
which is available under https://opus4.kobv.de/opus4-zib/files/6037/Presolve.pdf.

Some details on how PaPILO works internally are presented in a talk given during the SCIP workshop 2020 which has been recorded
and is available under https://www.youtube.com/watch?v=JKAyyWXGeQM.

# Licensing

To avoid confusion about licensing a short note on the LGPLv3.
This note is just an explanation and legally only the license text itself is of relevance.

When PaPILO is used as a header-based library then only section 3 of LGPLv3 is relevant, and not section 4. Therefore PaPILO in that setting could be used in a software that is distributed under the terms of a different license when the conditions of section 3 are met, which are

    a) Give prominent notice with each copy of the object code (refers to binary distributions of your software) that the Library (refers to PaPILO) is used in it and that the Library and its use are covered by this License (refers to the LGPLv3).
    b) Accompany the object code with a copy of the GNU GPL and this license document.

Modifications of PaPILO itself, however, must be distributed under LGPLv3.

For other licensing options we refer to https://scipopt.org/, where PaPILO can be obtained as part of the SCIP Optimization Suite.

# Contributors

Leona Gottwald ([@lgottwald](https://github.com/lgottwald)) --- main author

[Ivet Galabova](https://sites.google.com/view/ivetgalabova) ([@galabovaa](https://github.com/galabovaa)) --- groundwork for dual postsolve, KKT checker

[Kathrin Halbig](https://en.www.math.fau.de/edom/team/katrin-halbig/) --- presolver for GCD based reductions on inequalities(simplifyineq)

Anass Meskini --- general development and contributions to substitution presolver in terms of internship

[Daniel Rehfeldt](https://www.zib.de/members/rehfeldt) ([@dRehfeldt](https://github.com/dRehfeldt)) --- core data structures and MPS parsing
