file(GLOB_RECURSE SOURCES /bioinfo/users/hcliment/anaconda3/lib/R/library/Matrix/include/*h
        /bioinfo/users/hcliment/anaconda3/lib/R/library/Rcpp/include/*h
        /bioinfo/users/hcliment/anaconda3/lib/R/library/RcppEigen/include/*h
        /bioinfo/users/hcliment/anaconda3/lib/R/include/*h)

link_directories(/bioinfo/users/hcliment/anaconda3/lib/R/library/Rcpp/libs/
        /bioinfo/users/hcliment/anaconda3/lib/R/lib/
        /bioinfo/users/hcliment/anaconda3/lib )

include_directories(/bioinfo/users/hcliment/anaconda3/lib/R/library/Matrix/include/
        /bioinfo/users/hcliment/anaconda3/lib/R/library/Rcpp/include
        /bioinfo/users/hcliment/anaconda3/lib/R/library/RcppEigen/include
        /bioinfo/users/hcliment/anaconda3/lib/R/include)

add_executable(DBtest_associations_gridsearch DBtest_associations_gridsearch.cpp ${SOURCES})
target_link_libraries(DBtest_associations_gridsearch CEasyGWAS R gfortran)

add_executable(DBreadNetwork DBreadNetwork.cpp ${SOURCES})
target_link_libraries(DBreadNetwork CEasyGWAS R gfortran)