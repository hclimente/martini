set(R_LIBRARIES "~/anaconda/lib")

file(GLOB_RECURSE SOURCES ${R_LIBRARIES}/R/library/Matrix/include/*h
		${R_LIBRARIES}/R/library/Rcpp/include/*h
		${R_LIBRARIES}/R/library/RcppEigen/include/*h
		${R_LIBRARIES}/R/include/*h)

link_directories(${R_LIBRARIES}/R/library/Rcpp/libs/
		${R_LIBRARIES}/R/lib/
		${R_LIBRARIES}
		${PROJECT_SOURCE_DIR}/src/scones2
		${PROJECT_SOURCE_DIR}/src/scones2/External/maxflow
		${PROJECT_SOURCE_DIR}/src/scones2/External/Cephes
		${PROJECT_SOURCE_DIR}/src/scones2/External/eigen)

include_directories(${R_LIBRARIES}/R/library/Matrix/include/
		${R_LIBRARIES}/R/library/Rcpp/include
		${R_LIBRARIES}/R/library/RcppEigen/include
		${R_LIBRARIES}/R/include
		${PROJECT_SOURCE_DIR}/src/scones2/src/
		${PROJECT_SOURCE_DIR}/src/scones2/External
		${PROJECT_SOURCE_DIR}/src/scones2/External/maxflow
		${PROJECT_SOURCE_DIR}/src/scones2/External/Cephes
		${PROJECT_SOURCE_DIR}/src/scones2/External/eigen)

# add_executable(DBtest_associations_gridsearch DBtest_associations_gridsearch.cpp ${SOURCES})
# target_link_libraries(DBtest_associations_gridsearch CEasyGWAS R gfortran)

add_executable(DBreadNetwork DBreadNetwork.cpp ${SOURCES})
target_link_libraries(DBreadNetwork CEasyGWAS R gfortran)