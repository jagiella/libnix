AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = nix-abc nix-gp
lib_LTLIBRARIES = libnix.la


nix_abc_SOURCES = abc.cpp
nix_abc_LDADD = -lgsl -lgslcblas -lgomp libnix.la

nix_gp_SOURCES = gpopt.cpp
nix_gp_LDADD = -lgsl -lgslcblas -lgomp libnix.la

libnix_la_SOURCES = gp.cpp model.cpp
libnix_la_LDFLAGS = -version-info 0:0:0