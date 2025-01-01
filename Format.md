## `.sol` file format

In order to postsolve PaPILO requires the solution of the reduced file. To parse the solution file PaPILO uses the `.sol` file format (as supported by solvers like [SCIP](https://www.scipopt.org/doc/html/reader__sol_8h.php), or [Gurobi](https://www.gurobi.com/documentation/current/refman/sol_format.html), and used for example in [MIPLIB](https://miplib.zib.de/)). 
For the use with PaPILO it:
- is a tab-separated file, (basically a TSV with the extension `.sol`), containing up to three columns
- makes use of the same variable and constraint names as an initially supplied `.mps` file
- is not bound to a specific number type (`obj(0)` is treated and understood in the same way as `obj(0.00000000000000)`)
- can contain primal or dual solution values, but not in the same file
- if variable does not appear in the solution file it is supposed to be zero

Each row can either be an "auxiliary" information (e.g., the objective value) or give the solution for a specific variable/constraint.

### Row types
The following types of rows are allowed:
- starting with `#`: indicates a comment line and will be ignored (in-line comments are not allowed)
- containing `=obj=` in the first column, indicates that the problems objective value is given in the second column
- containing the name of a variable or constraint in the first column, indicating that a respective solution value is given in the second column, with the third column optionally - for variables - containing `obj(...)` with `...` being any integer or floating-point number (as documented by [SCIP](https://www.scipopt.org/doc/html/reader__sol_8h.php), the latter contains the  _"objective coefficient of variable"_

### Examples
The following gives a few short examples to better visualize the file format.

#### Primal solution
This may contain a primal solution that you may pass to PaPILO's `postsolve` command:
```txt
x    0.0
y    0.12508
z    -8815
```
Note that you should provide a proper primal result for all variables contained in the reduced `.mps` that PaPILO generated. Giving an objective value is not mandatory, PaPILO will derive that automatically, but you may want to check out the `--reference-objective` (`-o`) flag to pass an objective value obtained from your solver as reference for validation (see `--help ` for more information).

> Output from PaPILO: After running `postsolve` (or if you are using `solve` in the first place), you will get a similar `.sol` file from PaPILO, with the addition that it will always contain `=obj=` before any variable values, and that it applies some column alignment to make the file more readable, which extends the "tab separation" into an undefined amount of spaces used to align columns, where reasonable.

#### Dual solution
> Important: To allow PaPILO to consider a dual solution in `postsolve`, you need to properly configure it during the `presolve` step (otherwise it may apply presolve-algorithms that are - currently - not "reversable" during `postsolve`). Consult [lp_presolvers_with_basis.set](https://github.com/scipopt/papilo/blob/main/settings/lp_presolvers_with_basis.set) (or `...without_basis.set` and pass a proper set of parameters using the `--parameter-settings` (`-p`) flag of `presolve`. You should see a information log containing `... with dual-postsolve activated` when running PaPILO. If it states `... with dual-postsolve de-activated` your configuration is not properly set.

Dual solutions can be passed in a similar fashion, but are split up into results related to variables (sometimes also called "reduced costs") and to constraints (sometimes also called "shadow prices"). For a full handling of all available dual results, you therefore need to pass two files, with the same format, using the parameters `--dual-reduced-solution` (constraints) and `--costs-reduced-solution` (variables).

For example, dual results linked to constraints of the original problem (the values of the associated dual variables) may then look as `.sol` file like:
```txt
constr_1    -0.0
constr_sum_xy    10.0
```

> Note: This is based on the convention that "reduced costs" are the values of the dual variables associated with the bounds of a variable (primal variables themselves do not have a "dual solution"). If a variable has an upper and lower bound, this refers to (if existing) the bound-constraint that is actually binding. If none is binding, the reduced cost is zero.

Note that PaPILO will include the leading `=obj=` row in both resulting dual solution `.sol` files ("dual" and "cost").

# `.bas` file format

Similar to the solution format, PaPILO uses the basis format (.bas) used also in [SoPlex](https://github.com/scipopt/soplex).
For the use with PaPILO it:
* is a tab-separated file, (basically a TSV with the extension `.bas`), each line containing up to three entries
* the first entry identifies the status, the following entries may describe a column-row pair
  * `LL` column is at the lower bound (no row should be present)
  * `UL` column is at the upper bound (no row should be present)
  * `XL` column is in the basis and row is at the left hand side
  * `XU` column is in the basis and row is at the right hand side
- makes use of the same column and row names as the previously supplied `.mps` file
- if a column does not appear in the basis file, it is considered non-basic at its lower bound if existent, or at value zero if free

### Example

```txt
NAME  soplex.bas
 UL x6
 XL x8 c1
 UL x13
```
Only column x8 is basic determined by the left hand side of row c1, the columns x6 and x13 are at their upper bounds, and all remaining columns are at their lower bounds or zero.
