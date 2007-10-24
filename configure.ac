AC_INIT(qdp-lapack, 1.0.0, kostas@jlab.org)
AC_CONFIG_AUX_DIR(config)

dnl check for compiler
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB
     
dnl initialise automake
AM_INIT_AUTOMAKE

dnl put defines in include config file
AC_CONFIG_HEADERS([include/qdp-lapack_config_internal.h])
dnl dnl **************************************************************
dnl dnl ***                 Get User Options                       ***
dnl dnl **************************************************************

dnl What type of wrappings to use
dnl NOTE: there is uncertainty here. Different OS might need different
dnl bindings lapack functions. Dunno what they are (the known set of unknowns).
dnl This might well need more expanding
AC_ARG_ENABLE(lapack,
  AC_HELP_STRING([--enable-lapack=<type>],
    [Use the type of lapack bindings specified by <type>. Currently, this can be lapack]),
  [LAPACK_BINDING="${enableval}"]
)

dnl with QDP -- QDP will override everything else
dnl since a parallel QDP will know its architecture
dnl precition and have all the flags for QMP set already
dnl in its CXXFLAGS etc etc
dnl 
dnl This flag and all the other flags are complimentary
dnl to each other therefore. Also --with-qdp is also 
dnl passed down the command line when compiling Chroma,
dnl so a sub-configure will also work
dnl we can also put the QDP++ packers in here too...
AC_ARG_WITH(qdp,
  AC_HELP_STRING(
     [--with-qdp=DIR],
     [Build package with QDP++, where QDP++ is installed in DIR.]
  ),
  [QDPXX_HOME="${with_qdp}"
   QDPXX_GIVEN="yes"
  ],
  [QDPXX_GIVEN="no"]
)

dnl dnl ******************** DONE WITH USER INPUT *********************
dnl dnl ***        Now Configuration                               ****
dnl dnl ***************************************************************
dnl 
if test "X${QDPXX_HOME}X" = "XX" ; then
  AC_PATH_PROG(QDPXX_CONFIG, [qdp++-config], [])
else
  AC_PATH_PROG(QDPXX_CONFIG, [qdp++-config], [], [${QDPXX_HOME}/bin:${PATH}])
fi

if test "X${QDPXX_CONFIG}X" = "XX" ; then
  AC_MSG_ERROR([QDP++ configuration program qdp++-config not found.])
fi

AC_MSG_NOTICE([Found QDP++ configuration program ${QDPXX_CONFIG}])
AC_SUBST(CXX, "`${QDPXX_CONFIG} --cxx`")
AC_MSG_NOTICE([QDP++ compiler: ${CXX}])
AC_SUBST(QDPXX_CXXFLAGS, "`${QDPXX_CONFIG} --cxxflags`")
AC_MSG_NOTICE([QDP++ compile flags: ${QDPXX_CXXFLAGS}])
AC_SUBST(QDPXX_LDFLAGS,  "`${QDPXX_CONFIG} --ldflags`")
AC_MSG_NOTICE([QDP++ linking flags: ${QDPXX_LDFLAGS}])
AC_SUBST(QDPXX_LIBS,     "`${QDPXX_CONFIG} --libs`")
AC_MSG_NOTICE([QDP++ libraries flags: ${QDPXX_LIBS}])
AC_SUBST(QDPXX_ND,       "`${QDPXX_CONFIG} --Nd`")
if test "X${QDPXX_ND}X" = "XX" ; then
  AC_MSG_ERROR([Cannot determine QDP++ spacetime dimension])
else
  AC_MSG_NOTICE([QDP++ spacetime dimension: ${QDPXX_ND}])
fi
AC_SUBST(QDPXX_NC,       "`${QDPXX_CONFIG} --Nc`")
AC_MSG_NOTICE([QDP++ number of colors: ${QDPXX_NC}])
AC_SUBST(QDPXX_NS,       "`${QDPXX_CONFIG} --Ns`")
if test "X${QDPXX_NS}X" = "XX" ; then
  AC_MSG_ERROR([Cannot determine number of spin components in QDP++])
else
  AC_MSG_NOTICE([QDP++ number of spin components: ${QDPXX_NS}])
fi

AC_SUBST(QDPXX_PARALLEL_ARCH, "`${QDPXX_CONFIG} --parallel-arch`")
AC_MSG_NOTICE([QDP++ parallel arch: ${QDPXX_PARALLEL_ARCH}])

AC_SUBST(QDPXX_PRECISION, "`${QDPXX_CONFIG} --precision`")
AC_MSG_NOTICE([QDP++ precision: ${QDPXX_PRECISION}])

dnl Try to compile a QDP++ program to check the --with options
AC_MSG_CHECKING([if we can compile/link a simple QDP++ program])
PAC_QDPXX_LINK_CXX_FUNC(
  ${QDPXX_CXXFLAGS},
  ${QDPXX_LDFLAGS},
  ${QDPXX_LIBS},
  ,
  ,
  [qdpxx_link_ok=yes],
  [qdpxx_link_ok=no]
)
if test "X${qdpxx_link_ok}X" = "XyesX"; 
then 
  AC_MSG_RESULT(yes)
else 
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([Cannot compile/link a program with QDP++. Use --with-qdp++=<dir> to select a working version.])
fi 


dnl ************************************************************************
dnl **** Lapack bindings                                                ****
dnl ************************************************************************

case "${LAPACK_BINDING}" in
lapack)
	AC_MSG_NOTICE([Configuring to use Lapack Fortran bindings])
	AC_DEFINE(FORTRAN_LAPACK_BINDING,[],[Use Lapack Fortran bindings])
        LAPACK_STUBS="no"
	;;
*)
	AC_MSG_NOTICE([Configuring to use dummy stubs. If you want a real binding, use the --enable-lapack option])
	AC_DEFINE(DUMMY_LAPACK_BINDING,[],[Use dummy bindings])
        LAPACK_STUBS="yes"
	;;
esac

dnl There is one of these for each case above in the LAPACK_BINDING
AM_CONDITIONAL(BUILD_DUMMY_LAPACK, [ test "X${LAPACK_STUBS}X" = "XyesX" ])
AM_CONDITIONAL(BUILD_FORTRAN_LAPACK, [ test "X${LAPACK_BINDING}X" = "XlapackX" ])

dnl dnl produce output
AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(include/Makefile)
AC_CONFIG_FILES(lib/Makefile)
AC_OUTPUT
