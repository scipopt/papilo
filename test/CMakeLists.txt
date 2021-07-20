find_package(Boost REQUIRED OPTIONAL_COMPONENTS serialization program_options iostreams)

if (TARGET papilolib)
    set(PAPILOLIB_TESTS
            "papilolib"
            "papilolib-no-rows"
            )
    set(PAPILOLIB_TEST_FILE PapiloLib.cpp)
    set(PAPILOLIB_TARGET papilolib)
else ()
    set(PAPILOLIB_TESTS "")
    set(PAPILOLIB_TEST_FILE "")
    set(PAPILOLIB_TARGET "")
endif ()

if (Boost_IOSTREAMS_FOUND AND Boost_SERIALIZATION_FOUND AND Boost_PROGRAM_OPTIONS_FOUND)
    configure_file(resources/dual_fix_neg_inf.postsolve resources/dual_fix_neg_inf.postsolve COPYONLY)
    configure_file(resources/dual_fix_pos_inf.postsolve resources/dual_fix_pos_inf.postsolve COPYONLY)
    configure_file(instances/dual_fix_neg_inf.mps resources/dual_fix_neg_inf.mps COPYONLY)
    set(BOOST_REQUIRED_TESTS
            "finding-the-right-value-in-postsolve-for-a-column-fixed-pos-inf"
            "finding-the-right-value-in-postsolve-for-a-column-fixed-neg-inf"
            "mps-parser-loading-simple-problem"
            )
    set(BOOST_REQUIRED_TEST_FILES
            papilo/core/PostsolveTest.cpp
            papilo/io/MpsParserTest.cpp
            )
else ()
    set(BOOST_REQUIRED_TESTS "")
    set(BOOST_REQUIRED_TEST_FILES)
endif ()

add_executable(unit_test TestMain.cpp
        papilo/InstancesTest.cpp

        papilo/core/MatrixBufferTest.cpp
        papilo/core/SparseStorageTest.cpp
        papilo/core/PresolveTest.cpp
        papilo/core/ProblemUpdateTest.cpp
        papilo/misc/VectorUtilsTest.cpp

        papilo/presolve/CoefficientStrengtheningTest.cpp
        papilo/presolve/ConstraintPropagationTest.cpp
        papilo/presolve/DualFixTest.cpp
        papilo/presolve/DominatedColsTest.cpp
        papilo/presolve/DualInferTest.cpp
        papilo/presolve/FixContinuousTest.cpp
        papilo/presolve/FreeVarSubstitutionTest.cpp
        papilo/presolve/ImplIntDetectionTest.cpp
        papilo/presolve/ParallelRowDetectionTest.cpp
        papilo/presolve/ParallelColDetectionTest.cpp
        papilo/presolve/ProbingTest.cpp
        papilo/presolve/SingletonColsTest.cpp
        papilo/presolve/SimpleProbingTest.cpp
        papilo/presolve/SimpleSubstitutionTest.cpp
        papilo/presolve/SimplifyInequalitiesTest.cpp
        papilo/presolve/SparsifyTest.cpp

        ${BOOST_REQUIRED_TEST_FILES}
        ${PAPILOLIB_TEST_FILE})

target_link_libraries(unit_test papilo ${PAPILOLIB_TARGET} ${Boost_LIBRARIES})

set(unit_tests
        "accurate-numerical-statistics"

        "matrix-buffer"
        "vector-comparisons"
        "matrix-comparisons"

        "replacing-variables-is-postponed-by-flag"
        "happy-path-replace-variable"
        "happy-path-substitute-matrix-coefficient-into-objective"
        "happy-path-aggregate-free-column"

        #ProblemUpdate
        "happy-path-presolve-singleton-row"

        "problem-comparisons"

        "happy-path-coefficient-strengthening"

        "happy-path-constraint-propagation"

        #DomCol
        "domcol-happy-path"
        "domcol-parallel-columns"
        "domcol-multiple-parallel-cols-generate_redundant-reductions"
        "domcol-multiple-column"

        # Fix Continuous
        "happy-path-presolve-fix-continuous"
        "happy-path-no-presolve-fix-continuous"

        #FreeVarSubstitution
        "happy-path-test-free-variable-detection"

        "happy-path-implied-integer-detection"

        #Parallel row Detection
        "parallel-row-unchanged"
        "parallel-row-two-equations-infeasible-second-row-dominant"
        "parallel-row-two-equations-infeasible-first-row-dominant"
        "parallel-row-two-equations-feasible-second-row-dominant"
        "parallel-row-two-equations-feasible-first-row-dominant"
        "parallel-row-two-inequalities-redundant-row-second-row-dominant"
        "parallel-row-two-inequalities-redundant-row-first-row-dominant"
        "parallel-row-two-inequalities-tighten-lower-bound-second-row-dominant"
        "parallel-row-two-inequalities-tighten-lower-bound-first-row-dominant"
        "parallel-row-two-inequalities-tighten-upper-bound-second-row-dominant"
        "parallel-row-two-inequalities-tighten-upper-bound-first-row-dominant"
        "parallel-row-two-inequalities-infeasible-first-row-dominant"
        "parallel-row-two-inequalities-infeasible-second-row-dominant"
        "parallel-row-two-inequalities-tighten-upper-bound-first-row-neg"
        "parallel-row-overwrite-inf-first-row-rhs-inf"
        "parallel-row-overwrite-inf-first-row-lhs-inf"
        "parallel-row-overwrite-inf-first-row-lhs-inf-neg-factor"
        "parallel-row-mixed-infeasible-first-row-equation"
        "parallel-row-mixed-second-row-equation"
        "parallel-row-mixed-infeasible-second-row-equation"
        "parallel-row-multiple-parallel-rows"

        #parallel_col_detection
        "parallel_col_detection_2_integer_columns"
        "parallel_col_detection_2_continuous_columns"
        "parallel_col_detection_int_cont_merge_possible"
        "parallel_col_detection_cont_int_merge_possible"
        "parallel_col_detection_cont_int_merge_failed"
        "parallel_col_detection_int_cont_merge_failed"
        "parallel_col_detection_int_merge_failed_hole"
        "parallel_col_detection_obj_not_parallel"
        "parallel_col_detection_multiple_parallel_columns"

        #Probing
        "happy-path-probing"
        "failed-path-probing-on-not-binary-variables"

        #singleton Column
        "happy-path-singleton-column"
        "happy-path-singleton-column-equation"
        "happy-path-singleton-column-implied-bounds-negative-coeff-pos-bounds"
        "happy-path-singleton-column-implied-bounds-negative-coeff-neg-bounds"
        "happy-path-singleton-column-implied-bounds-positive-coeff-pos-bounds"
        "happy-path-singleton-column-implied-bounds-positive-coeff-neg-bounds"

        #Simple Probing
        "happy-path-simple-probing"
        "happy-path-simple-probing-only-negative-coeff"
        "happy-path-simple-probing-only-binary-negative-coefficient"
        "happy-path-simple-probing-only-binary-positive-coefficient"


        #Simple Substitution
        "happy-path-simple-substitution-for-2-int"
        "happy-path-simple-substitution-for-2-continuous"
        "happy-path-simple-substitution-for-continuous-and-integer"
        "happy-path-simple-substitution-for-int-continuous-coeff"
        "example_10_1_in_constraint_integer_programming"
        "should_return_infeasible_if_gcd_of_coeff_is_in_rhs"
        "should_return_feasible_if_gcd_of_coeff_is_in_rhs"

        #DualFix
        "happy-path-dual-fix"
        "happy_path_dual_substitution_for_equations"
        "happy_path_dual_substitution"
        "failed-path-simple-substitution-for-2-int"
        "happy_path_dual_substitution_unbounded_variables"
        "trivial-column-presolve-does-dual-presolve-already"
        "happy_path_dual_substitution_rounding"
        "happy_path_dual_fix_on_infinity"

        "happy-path-simplify-inequalities-only-greatest-divisor"
        "simplify_inequ_doesnt_lock_more_rows"
        "simplify_inequ_doesnt_apply_lb_and_ub_on_one_row"

        #Sparsify
        "happy-path-sparsify"
        "happy-path-sparsify-two-equalities"
        "failed-path-sparsify"

        "example-1-from-3.5-Presolve-Reductions-in-MIP" # not working
        "example-2-from-3.5-Presolve-Reductions-in-MIP" # not working
        "example-3-from-3.6-Presolve-Reductions-in-MIP" # trivial presolve fixes this case
        "example-4-from-4.4-Presolve-Reductions-in-MIP" # not working
        "example-5-from-4.4-Presolve-Reductions-in-MIP"
        "example-10-from-7.5-Presolve-Reductions-in-MIP"


        "integration-test-for-flugpl"
        ${PAPILOLIB_TESTS}
        ${BOOST_REQUIRED_TESTS}
        )

# add a test to build the SCIP binary that all further tests depend on
add_test(NAME unit-test-build
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target unit_test
        )

set_tests_properties(unit-test-build
        PROPERTIES
        RESOURCE_LOCK unittestbin)

foreach (test ${unit_tests})
    add_test(NAME unit-test-${test} COMMAND unit_test ${test})

    set_tests_properties(unit-test-${test}
            PROPERTIES
            DEPENDS unit-test-build)
endforeach ()

if (TARGET papilo-executable AND PAPILO_HAVE_SCIP)


    # MIP instances
    set(instances_MIP
            "instances/bell5.mps\;8966406.49152"
            "instances/blend2.mps\;7.598985"
            "instances/dcmulti.mps\;188182"
            "instances/egout.mps\;568.1007"
            "instances/enigma.mps\;0"
            "instances/flugpl.mps\;1201500"
            "instances/gt2.mps\;21166"
            "instances/lseu.mps\;1120"
            "instances/misc03.mps\;3360"
            "instances/p0548.mps\;8691"
            "instances/rgn.mps\;82.19999924"
            "instances/dual_fix_pos_inf.mps\;-7.5"
            "instances/dual_fix_neg_inf.mps\;-10"
            )

    set(infeas
            "instances/presolved_ns2080781.mps"
            )

    set(numtypes
            "d"
            "q"
            "r")

    # add a test to build the SCIP binary that all further tests depend on
    add_test(NAME papilo-build
            COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target papilo-executable
            )
    # avoid that several build jobs try to concurrently build the binaries
    set_tests_properties(papilo-build
            PROPERTIES
            RESOURCE_LOCK papilobin)

    # macro to split an instance into its relevant information
    # - path
    # - optval
    # - basename
    macro(split_instance instance)
        list(GET instance 0 path)
        list(GET instance 1 optval)
        get_filename_component(basename ${path} NAME)
    endmacro(split_instance)

    ##solving these instances returns the defined optimal value for every numtype
    macro(add_instancetests instances numtypes)
        foreach (instance ${${instances}})
            split_instance(instance)
            file(TO_NATIVE_PATH "${PROJECT_SOURCE_DIR}/test/${path}" instance_file)
            foreach (numtype ${${numtypes}})
                add_test(NAME ${numtype}-solve-${basename}
                        COMMAND $<TARGET_FILE:papilo-executable> solve -a ${numtype} -f ${instance_file} -o ${optval}
                        )
                set_tests_properties(${numtype}-solve-${basename}
                        PROPERTIES
                        PASS_REGULAR_EXPRESSION "validation: SUCCESS"
                        DEPENDS papilo-build
                        )
            endforeach (numtype)
        endforeach (instance)
    endmacro(add_instancetests)

    ##solving returns infeasible for these instances and every potential numtype
    macro(add_infeasible_instancetests instances numtypes)
        foreach (instance ${${instances}})
            file(TO_NATIVE_PATH "${PROJECT_SOURCE_DIR}/test/${instance}" instance_file)
            foreach (numtype ${${numtypes}})
                add_test(NAME ${numtype}-infeas-${instance}
                        COMMAND $<TARGET_FILE:papilo-executable> solve -a ${numtype} -f ${instance_file}
                        )
                set_tests_properties(${numtype}-infeas-${instance}
                        PROPERTIES
                        PASS_REGULAR_EXPRESSION "presolve detected infeasible problem"
                        DEPENDS papilo-build
                        )
            endforeach (numtype)
        endforeach (instance)
    endmacro(add_infeasible_instancetests)


    add_instancetests(instances_MIP numtypes)
    add_infeasible_instancetests(infeas numtypes)
endif ()



set(filename_for_solutions "solutions.test")
# the following code block checks if a solution is contained in a reduced problem.
# The input for solutions.test file can be done via a file with the following format:
# PATH_TO_MPS_FILE_1,PATH_TO_SOL_FILE_1
# PATH_TO_MPS_FILE_2,PATH_TO_SOL_FILE_2
# ...
if (TARGET papilo-executable AND EXISTS "${PROJECT_SOURCE_DIR}/${filename_for_solutions}")

    FILE(READ "${PROJECT_SOURCE_DIR}/${filename_for_solutions}" content)
    string(REPLACE "\n" ";" instances ${content})

    macro(split_solution solution)
        string(REPLACE "," ";" sol ${solution})
        list(GET sol 0 path)
        list(GET sol 1 path_solution)
        get_filename_component(basename ${path} NAME)
    endmacro(split_solution)

    #solution is still contained after presolving
    macro(validate_solution instance settings)
        split_solution(${instance})
        file(TO_NATIVE_PATH "${path}" instance_file)
        file(TO_NATIVE_PATH "${path_solution}" solution_file)
        file(TO_NATIVE_PATH "${PROJECT_SOURCE_DIR}/settings/${settings}" settings_file)
        #TODO find better testname -> solution file can be occur multiple times
        add_test(NAME validate-solution-${settings}-${solution_file}
                COMMAND $<TARGET_FILE:papilo-executable> presolve -f ${instance_file} -b ${solution_file} -p ${settings_file}
                )
        set_tests_properties(validate-solution-${settings}-${solution_file}
                PROPERTIES
                PASS_REGULAR_EXPRESSION "validation: SUCCESS\n"
                DEPENDS papilo-build
                )
    endmacro(validate_solution)

    foreach (instance ${instances})
        validate_solution(${instance}, "no_strong_dual_reductions.set")
#        validate_solution(${instance}, "scip_no_strong_dual_reductions.set")
    endforeach (instance)
endif ()
