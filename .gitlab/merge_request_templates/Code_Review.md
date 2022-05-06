### Description

[Add a description of the introduced changes here.]

### Code review

* [ ] The code change is correct.
* [ ] The naming and place of new methods is clear and consistent **or** no new methods have been added.
* [ ] The code is sufficiently documented.
* [ ] The CHANGELOG is up to date (including API changes if present in this MR).
* [ ] The new code is sufficiently covered by tests (perhaps, new coverage settings or new unit tests have been added).

### Additional Testing (depending on the changes)

* [ ] PaPILO with Boost version 1.78 passes `jenkin rerun boost78`
* [ ] Supported integrated solver (besides SCIP and SoPlex) pass `jenkins rerun build_with_{gurobi,glop}`). _highs currently disabled_
* [ ] PaPILO without TBB is build `jenkins rerun build_no_TBB`

### Does this merge request introduce an API change? :warning:

* [ ] No, **or** as far as possible, the code ensures backwards compatibility.
* [ ] No, **or** the `PAPILO_MINOR_VERSION` or `PAPILO_MAJOR_VERSION` is updated **AND** `SCIP` and `SoPlex` can handle these changes.
