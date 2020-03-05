#add_executable(sparse_storage SparseStorage.cpp)
#add_executable(common Common.cpp)
#add_executable(coef_stren CoefStren.cpp)
#add_executable(aggregation Aggregation.cpp)
#todo add if fixed
#add_executable(parallel_rows ParallelRows.cpp)

#target_link_libraries(sparse_storage ${TBB_IMPORTED_TARGETS})
#target_link_libraries(common ${TBB_IMPORTED_TARGETS})
#target_link_libraries(coef_stren ${TBB_IMPORTED_TARGETS})
#target_link_libraries(aggregation ${TBB_IMPORTED_TARGETS})
#target_link_libraries(parallel_rows ${TBB_IMPORTED_TARGETS})

#add_test(NAME SparseStorage COMMAND sparse_storage )
#add_test(NAME Common COMMAND common )
#add_test(NAME CoefStren COMMAND coef_stren)
#add_test(NAME Aggregation COMMAND aggregation)
#add_test(NAME ParallelRows COMMAND parallel_rows)

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

add_executable(unit_test TestMain.cpp MatrixBuffer.cpp VectorUtils.cpp ${PAPILOLIB_TEST_FILE})
target_link_libraries(unit_test papilo ${PAPILOLIB_TARGET})

set(unit_tests
 "matrix-buffer"
 "vector-comparisons"
 ${PAPILOLIB_TESTS}
  )

#
# add a test to build the SCIP binary that all further tests depend on
#
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
   #
   # MIP instances
   #
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

