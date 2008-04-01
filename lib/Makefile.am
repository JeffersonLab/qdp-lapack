# $Id: Makefile.am,v 1.5 2008-04-01 19:59:46 kostas Exp $
#
# Path for all files used in building the library

#
# Any other compiler flags
# 
#AM_CFLAGS = -I@top_srcdir@/include -I@top_builddir@/include @QDPXX_CFLAGS@ @CFLAGS@
AM_CFLAGS = -I@top_srcdir@/include -I@top_builddir@/include @CFLAGS@
AM_CXXFLAGS = -I@top_srcdir@/include -I@top_builddir@/include @QDPXX_CXXFLAGS@ @CXXFLAGS@
AM_CC  = $(CC)
AM_CXX = $(CXX)

lib_LIBRARIES = libqdp-lapack.a

#
# An AM_CONDITIONAL if statement to decide which fermion type
#
libqdp_lapack_a_SOURCES =

if BUILD_DUMMY_LAPACK
libqdp_lapack_a_SOURCES += lapack_stub.cc 
endif

if BUILD_FORTRAN_LAPACK
libqdp_lapack_a_SOURCES += fortran_lapack.cc Complex.c  IncrEigpcg.c         numerical.c  restart_X.c eigpcg.c   ortho.c
endif

