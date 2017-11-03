AC_INIT(martini, version-0.1)
AC_LANG(C++)
AC_PROG_CXX

dnl check that gin is correctly installed
AC_CHECK_HEADERS(gin/feature_selection/shake.h)
AC_SEARCH_LIBS([libgin_present], [gin, gin-devel], [], [AC_MSG_ERROR([gin not found])])

AC_ARG_VAR([GIN_PATH], "Path to the gin dynamic library.")

dnl fill rpath argument when necessary
RPATH_FLAG=''
unamestr=`uname`
if [[ "$unamestr" == 'Darwin' ]]; then
   RPATH_FLAG="-rpath $GIN_PATH"
fi

AC_SUBST([RPATH_FLAG], $RPATH_FLAG)


AC_OUTPUT(src/Makevars)