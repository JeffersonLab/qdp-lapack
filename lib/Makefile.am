# $Id: Makefile.am,v 1.1 2007-10-11 20:21:11 edwards Exp $
#
# Path for all files used in building the library

#
# Any other compiler flags
# 
AM_CFLAGS = -I@top_srcdir@/include -I@top_builddir@/include @QDPXX_CFLAGS@ @CFLAGS@ 
AM_CXXFLAGS = -I@top_srcdir@/include -I@top_builddir@/include @QDPXX_CXXFLAGS@ @CXXFLAGS@
AM_CC  = $(CC)
AM_CXX = $(CXX)

lib_LIBRARIES = libqdp-lapack.a

#
# An AM_CONDITIONAL if statement to decide which fermion type
#
libqdp_lapack_a_SOURCES =

if BUILD_STUBS
libqdp_lapack_a_SOURCES += lapack_stub.cc
endif

if BUILD_CLAPACK
libqdp_lapack_a_SOURCES += clapack.cc
endif

