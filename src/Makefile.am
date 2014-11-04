AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = nix-gp nix-tumor1d
lib_LTLIBRARIES = libnix.la

nix_gp_SOURCES = gpopt.cpp
nix_gp_LDADD = -lgsl -lgslcblas libnix.la
nix_gp_LDFLAGS = -fopenmp

nix_tumor1d_SOURCES = tumor1d.cpp
nix_tumor1d_LDADD = libnix.la


libnix_la_SOURCES = gp.cpp model.cpp hash.cpp statistics.cpp optimization.cpp Mathematix.cpp testfunctions.cpp
libnix_la_LDFLAGS = -version-info 0:0:0