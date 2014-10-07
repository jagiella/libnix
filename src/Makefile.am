AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = abc gp
lib_LTLIBRARIES = libnix.la


abc_SOURCES = abc.cpp
abc_LDADD = -lgsl -lgslcblas -lgomp libnix.la

gp_SOURCES = gpopt.cpp
gp_LDADD = -lgsl -lgslcblas -lgomp libnix.la

libnix_la_SOURCES = gp.cpp model.cpp
libnix_la_LDFLAGS = -version-info 0:0:0