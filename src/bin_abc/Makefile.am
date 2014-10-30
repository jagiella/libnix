AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = nix-abc

nix_abc_SOURCES = bin_abc.cpp
nix_abc_LDADD = -lgsl -lgslcblas  ../libnix.la
nix_abc_LDFLAGS = -fopenmp
