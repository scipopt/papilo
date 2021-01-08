if(TARGET papilolib)
   set(PAPILOLIB_TESTS
      "papilolib"
      "papilolib-no-rows"
      )
   set(PAPILOLIB_TEST_FILE PapiloLib.cpp)
   set(PAPILOLIB_TARGET papilolib)
else()
   set(PAPILOLIB_TESTS "")
   set(PAPILOLIB_TEST_FILE "")
   set(PAPILOLIB_TARGET "")
endif()

add_executable(unit_test TestMain.cpp
        papilo/InstancesTest.cpp
        papilo/DebugTest.cpp

        papilo/core/MatrixBufferTest.cpp
        papilo/core/SparseStorageTest.cpp
        papilo/core/ProblemUpdateTest.cpp
        papilo/misc/VectorUtilsTest.cpp

        papilo/presolve/CoefficientStrengtheningTest.cpp
        papilo/presolve/ConstraintPropagationTest.cpp
        papilo/presolve/DualFixTest.cpp
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

        ${PAPILOLIB_TEST_FILE})

target_link_libraries(unit_test papilo ${PAPILOLIB_TARGET} ${Boost_LIBRARIES})


set(unit_tests
        "accurate-numerical-statistics"

        "matrix-buffer"
        "vector-comparisons"
        "matrix-comparisons"

        #ProblemUpdate
        "happy-path-presolve-singleton-row"

        "problem-comparisons"

        "happy-path-coefficient-strengthening"

        "happy-path-constraint-propagation"

        # Fix Continuous
        "happy-path-presolve-fix-continuous"
        "happy-path-no-presolve-fix-continuous"

        #FreeVarSubstitution
        "test-free-variable-detection"

        "happy-path-implied-integer-detection"

        "happy-path-parallel-row-detection"

        "happy-path-parallel-column-detection"

        #Probing
        "happy-path-probing"
        "failed-path-probing-on-not-binary-variables"

        #singleton Column
        "happy-path-singleton-column"
        "happy-path-singleton-column-equation"

        #Simple Probing
        "happy-path-simple-probing"

        #Simple Substitution
        "happy-path-simple-substitution-for-2-int"
        "happy-path-simple-substitution-for-2-continuous"
        "happy-path-simple-substitution-for-continuous-and-integer"

        #DualFix
        "happy-path-dual-fix"
        "happy_path_dual_substitution_for_equations"
        "happy_path_dual_substitution"
        "failed-path-simple-substitution-for-2-int"
        "happy-path-simplify-inequalities-only-greatest-divisor"
        "happy_path_dual_substitution_unbounded_variables"
        "trivial-column-presolve-does-dual-presolve-already"
        "happy_path_dual_substitution_rounding"

        "example-1-from-3.5-Presolve-Reductions-in-MIP" # not working
        "example-2-from-3.5-Presolve-Reductions-in-MIP" # not working
        "example-3-from-3.6-Presolve-Reductions-in-MIP" # trivial presolve fixes this case
        "example-4-from-4.4-Presolve-Reductions-in-MIP" # not working
        "example-5-from-4.4-Presolve-Reductions-in-MIP"
        # "example-6-from-4.5-Presolve-Reductions-in-MIP"
        # "example-7-from-5.5-Presolve-Reductions-in-MIP"
        "example-8-from-6.3-Presolve-Reductions-in-MIP" # not working because one Coefficients are excluded
        # "example-9-from-7.1-Presolve-Reductions-in-MIP"
        "example-10-from-7.5-Presolve-Reductions-in-MIP"


        ${PAPILOLIB_TESTS}
        )

# add a test to build the SCIP binary that all further tests depend on
add_test(NAME unit-test-build
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target unit_test
        )

set_tests_properties(unit-test-build
                     PROPERTIES
                     RESOURCE_LOCK unittestbin)

foreach(test ${unit_tests})
   add_test(NAME unit-test-${test} COMMAND unit_test ${test} )

   set_tests_properties(unit-test-${test}
                        PROPERTIES
                        DEPENDS unit-test-build)
endforeach()

if(TARGET papilo-executable AND PAPILO_HAVE_SCIP)

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
       )

   set(numtypes
      "d"
      "q"
      "r")

   #
   # add a test to build the SCIP binary that all further tests depend on
   #
   add_test(NAME papilo-build
           COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target papilo-executable
           )
   # avoid that several build jobs try to concurrently build the binaries
   set_tests_properties(papilo-build
                        PROPERTIES
                        RESOURCE_LOCK papilobin)

   #
   # macro to split an instance into its relevant information
   # - path
   # - optval
   # - basename
   #
   macro(split_instance instance)
       list(GET instance 0 path)
       list(GET instance 1 optval)
       get_filename_component(basename ${path} NAME)
   endmacro(split_instance)

   macro(add_instancetests instances numtypes)
      #
      # loop over the instances
      #
      foreach(instance ${${instances}})
         split_instance(instance)

         foreach(numtype ${${numtypes}})
            #
            # treat the instance as a tuple (list) of two values
            #
            file(TO_NATIVE_PATH ${PROJECT_SOURCE_DIR}/test/${path} instance_file)
            add_test(NAME ${numtype}-solve-${basename}
                     COMMAND $<TARGET_FILE:papilo-executable> solve -a ${numtype} -f ${instance_file} -o ${optval}
                     )
            set_tests_properties(${numtype}-solve-${basename}
                                 PROPERTIES
                                    PASS_REGULAR_EXPRESSION "validation: SUCCESS"
                                    DEPENDS papilo-build
                                    )
         endforeach(numtype)
      endforeach(instance)
   endmacro(add_instancetests)

   add_instancetests(instances_MIP numtypes)
endif()
