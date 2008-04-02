# $Id: Makefile.am,v 1.6 2008-04-02 02:45:53 kostas Exp $
#
# Path for all files used in building the library

#
# Any other compiler flags
# 
#AM_CFLAGS = -I@top_srcdir@/include -I@top_builddir@/include @QDPXX_CFLAGS@ @CFLAGS@
AM_CFLAGS = -I@top_srcdir@/include -I@top_builddir@/include @CFLAGS@ -DUSE_QMP
AM_CXXFLAGS = -I@top_srcdir@/include -I@top_builddir@/include @QDPXX_CXXFLAGS@ @CXXFLAGS@
AM_CC  = $(CC)
AM_CXX = $(CXX)

lib_LIBRARIES = libqdp-lapack.a

#
# An AM_CONDITIONAL if statement to decide which fermion type
#
libqdp_lapack_a_SOURCES =

#if USE_QMP2
#AM_CFLAGS += -DUSE_QMP
#endif

if BUILD_DUMMY_LAPACK
libqdp_lapack_a_SOURCES += lapack_stub.cc 
endif

if BUILD_FORTRAN_LAPACK
libqdp_lapack_a_SOURCES += fortran_lapack.cc Complex.c  IncrEigpcg.c  numerical.c  restart_X.c eigpcg.c   ortho.c
endif

