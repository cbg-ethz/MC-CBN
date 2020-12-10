dnl Initialise Autoconf
## require at least autoconf 2.69
AC_PREREQ([2.69])


AC_INIT(
	[mccbn],
	m4_esyscmd_s([awk '/^Version:/ {print $2}' DESCRIPTION]),
	[],
	[mccbn],
	[https://github.com/cbg-ethz/MC-CBN])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_SRCDIR([src/mcem.cpp])

dnl Setup RBin
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}";
then
  AC_MSG_ERROR([no R_HOME found])
else
  AC_MSG_NOTICE([R home found in ${R_HOME}])

  R_BIN="${R_HOME}/bin/R"
  CXX11=`"${R_BIN}" CMD config CXX11`
  if test -z "$CXX11"; then
    AC_MSG_ERROR([No C++11 compiler is available])
  fi
  CXX11STD=`"${R_BIN}" CMD config CXX11STD`
  CXX="${CXX11} ${CXX11STD}"
  CXXFLAGS=`"${R_BIN}" CMD config CXX11FLAGS`
  CPPFLAGS=`"${R_BIN}" CMD config CPPFLAGS`
  AC_MSG_NOTICE([setting CXX to ${CXX}])
fi

# We are using C++
AC_LANG(C++)
AC_REQUIRE_CPP

# Check the C++ compiler using the CXX value set
AC_PROG_CXX
AC_PROG_SED
AC_PROG_RANLIB

AC_ARG_ENABLE(debug,
AS_HELP_STRING([--enable-debug],
               [enable debugging, default: no]),
[case "${enableval}" in
             yes) debug=yes ;;
             no)  debug=no ;;
             *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
esac],
[debug=no])
if (test "x$debug" = "xyes"); then
       AC_DEFINE(DEBUG, 1, [Define if debug code should be enabled.])
       GLOBAL_CFLAGS+=" -g"
fi

AC_SEARCH_LIBS([cblas_dgemm], [openblas],
[
	AC_MSG_NOTICE([using openblas])
	mccbn_LDLIBS="${mccbn_LDLIBS} -lopenblas"
],
[
	AC_MSG_NOTICE([using rblas])
	RBLAS=`"${R_BIN}" CMD config BLAS_LIBS`
	mccbn_LDLIBS="${mccbn_LDLIBS} $RBLAS"
])

RLAPACK=`"${R_BIN}" CMD config LAPACK_LIBS`
RFORT=`"${R_BIN}" CMD config FLIBS`
mccbn_LDLIBS="${mccbn_LDLIBS} $RLAPACK $RFORT"

# Other dependencies
AC_OPENMP
AX_BOOST_BASE([1.59.0], [], [AC_MSG_ERROR([did not find the boost headers])])
AC_ARG_WITH([mklcxxflags], AS_HELP_STRING([--with-mkl], [Compiler flags for Intel MKL (optional)]))
AC_ARG_WITH([mklldflags], AS_HELP_STRING([--with-mklldflags], [Linker flags for Intel MKL (optional)]))

# Add MKL FLAGS
AS_IF([test "x$with_mklcxxflags" != "xno" && test "x$with_mklldflags" != "xno"], [
	AX_APPEND_COMPILE_FLAGS([$with_mklcxxflags], [CXXFLAGS])
	CXXFLAGS=$( echo ${CXXFLAGS} | $SED -e 's/^ *//' -e 's/ *$//' )

  AX_APPEND_LINK_FLAGS([$with_mklldflags], [LDFLAGS])
  LDFLAGS=$( echo ${LDFLAGS} | $SED -e 's/^ *//' -e 's/ *$//' )

  AC_DEFINE(MKL_ENABLED, [1], "Use Intel MKL")
  AC_SUBST(MKL_CXXFLAGS, [$with_mklcxxflags])
  AC_SUBST(MKL_LDFLAGS, [$with_mklldflags])
])

# Write the flags into the src/Makevars file
AC_SUBST([PKG_CPPFLAGS], ["${PKG_CPPFLAGS}"])
AC_SUBST([PKG_LIBS], ["${LIBS} ${PKG_LIBS}"])
AC_SUBST([mccbn_LDLIBS])
AC_SUBST([mccbn_CPPFLAGS])
AC_SUBST([OPENMP_CXXFLAGS])

AC_CONFIG_FILES([src/Makevars Makefile])
AM_INIT_AUTOMAKE([foreign])
AC_OUTPUT

AC_MSG_RESULT([
    $PACKAGE_NAME

	CXX:            ${CXX}
	CXXFLAGS:       ${CXXFLAGS} ${OPENMP_CXXFLAGS}
	CPPFLAGS:       ${CPPFLAGS} ${mccbn_CPPFLAGS} ${PKG_CPPFLAGS}
	Boost_CPPFLAGS: ${BOOST_CPPFLAGS}
  MKL_CXXFLAGS:   ${MKL_CXXFLAGS}

	LIBS:           ${PKG_LIBS} ${mccbn_LDLIBS} ${BOOST_LDFLAGS}
	MKL_LDFLAGS:    ${MKL_LDFLAGS}
])
