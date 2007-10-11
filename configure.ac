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

dnl use c-lapack wrappings
AC_ARG_WITH(clapack,
  AC_HELP_STRING(
    [--with-clapack=DIR],
    [Use the BAGEL Clover Apply Library installed in DIR]
  ),
    [bagel_clover_enabled="yes"
     BAGEL_CLOVER_HOME="${with_bagel_clover}"]
)

AC_ARG_ENABLE(parallel-arch,
  AC_HELP_STRING([--enable-parallel-arch=<arch>],
    [Build Assembler for parallel architecture <arch>]),
  [ PARALLEL_ARCH="${enableval}"]
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
     [ Build Assembler with QDP++, where QDP++ is installed in DIR. ]
  ),
  [ QDPXX_HOME="${with_qdp}"
    QDPXX_GIVEN="yes"
  ],
  [ QDPXX_GIVEN="no" ]
)

dnl dnl ******************** DONE WITH USER INPUT *********************
dnl dnl ***        Now Configuration                               ****
dnl dnl ***************************************************************
dnl 
dnl     ***************************************************************
dnl     ***       There are two distinct branches depending         ***
dnl     ***      on whether we are compiling with QDP++ or not      ***
dnl     ***************************************************************

dnl 
if test "X${QDPXX_GIVEN}X" = "XyesX"; then 
     dnl Find the QDP++ config program 
     if test "X${QDPXX_HOME}X" = "XX" ; 
     then
	AC_PATH_PROG(QDPXX_CONFIG, [qdp++-config], [])
     else
	AC_PATH_PROG(QDPXX_CONFIG, [qdp++-config], [], [${QDPXX_HOME}/bin:${PATH}])
     fi

     if test "X${QDPXX_CONFIG}X" = "XX" ; then
	AC_MSG_ERROR([QDP++ configuration program qdp++-config not found.])
     fi

     AC_MSG_NOTICE([Found QDP++ configuration program ${QDPXX_CONFIG}])
     AC_MSG_CHECKING([if QDP++ uses checkerboarded layout])
     layout=`${QDPXX_CONFIG} --layout`
     if test "X${layout}X" == "Xcb2X"; then 
	AC_MSG_RESULT([yes])
     else
        AC_MSG_RESULT([no])
	AC_MSG_ERROR([This Dslash Operator only works for checkerboarded QDP++. Sorry])
     fi

     AC_SUBST(CC, ${CC})
     AC_MSG_NOTICE([QDP++ C compiler: ${CC}])

     AC_SUBST(CXX,"`${QDPXX_CONFIG} --cxx`")
     AC_MSG_NOTICE([QDP++ C++ compiler: ${CXX}])

     AC_SUBST(SSEXX_AUX_CFLAGS, "`${QDPXX_CONFIG} --cxxflags`")
     AC_MSG_NOTICE([QDP++ compile flags: ${SSEXX_AUX_CFLAGS}])

     AC_SUBST(SSEXX_AUX_LDFLAGS,  "`${QDPXX_CONFIG} --ldflags`")
     AC_MSG_NOTICE([QDP++ linking flags: ${SSEXX_AUX_LDFLAGS}])

     AC_SUBST(SSEXX_AUX_LIBS,     "`${QDPXX_CONFIG} --libs`")
     AC_MSG_NOTICE([QDP++ libraries flags: ${SSEXX_AUX_LIBS}])

     QDPXX_ND="`${QDPXX_CONFIG} --Nd`"
     AC_MSG_NOTICE([QDP++ spacetime dimension: ${QDPXX_ND}])
     
     QDPXX_NC="`${QDPXX_CONFIG} --Nc`"
     AC_MSG_NOTICE([QDP++ number of colors: ${QDPXX_NC}])

     QDPXX_NS="`${QDPXX_CONFIG} --Ns`"
     AC_MSG_NOTICE([QDP++ number of spin components: ${QDPXX_NS}])

     PARALLEL_ARCH="`${QDPXX_CONFIG} --parallel-arch`"
     AC_MSG_NOTICE([QDP++ has parallel arch: ${PARALLEL_ARCH} ] )

     PRECISION="`${QDPXX_CONFIG} --precision`"
     AC_MSG_NOTICE([QDP++ has precision: ${PRECISION} ] )

     dnl This function is defined in acinclude.m4
     dnl it tries to try to compile the program 
     PAC_QDPXX_LINK_CXX_FUNC(
	${SSEXX_AUX_CFLAGS},
	${SSEXX_AUX_LDFLAGS},
	${SSEXX_AUX_LIBS},
	,
	,
	[qdpxx_link_ok=yes],
	[qdpxx_link_ok=no]
     )

     AC_MSG_CHECKING([if we can compile/link a simple QDP++ program])
     if test "X${qdpxx_link_ok}X" = "XyesX" ; then 
	AC_MSG_RESULT(yes)
     else 
	AC_MSG_RESULT(no)
	AC_MSG_ERROR([Cannot compile/link a program with QDP++.
    Use --with-qdp++=<dir> to select a working version.])
    fi 

    dnl *********************************************************
    dnl ***              end of section with QDP++            ***
    dnl *********************************************************
else 

dnl dnl produce output
AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(include/Makefile)
AC_CONFIG_FILES(lib/Makefile)
AC_OUTPUT
