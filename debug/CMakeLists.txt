file(GLOB_RECURSE SOURCES /Users/hclimente/anaconda/lib/R/library/Matrix/include/*h
        /Users/hclimente/anaconda/lib/R/library/Rcpp/include/*h
        /Users/hclimente/anaconda/lib/R/library/RcppEigen/include/*h
        /Users/hclimente/anaconda/lib/R/include/*h)

link_directories(/Users/hclimente/anaconda/lib/R/library/Rcpp/libs/
        /Users/hclimente/anaconda/lib/R/lib/
        /Users/hclimente/anaconda/lib )

include_directories(/Users/hclimente/anaconda/lib/R/library/Matrix/include/
        /Users/hclimente/anaconda/lib/R/library/Rcpp/include
        /Users/hclimente/anaconda/lib/R/library/RcppEigen/include
        /Users/hclimente/anaconda/lib/R/include)

add_executable(DBtest_associations_gridsearch DBtest_associations_gridsearch.cpp ${SOURCES})
target_link_libraries(DBtest_associations_gridsearch CEasyGWAS R gfortran)

add_executable(DBreadNetwork DBreadNetwork.cpp ${SOURCES})
target_link_libraries(DBreadNetwork CEasyGWAS R gfortran)