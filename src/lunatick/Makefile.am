AUTOMAKE_OPTIONS = subdir-objects

#bin_PROGRAMS = lunatick
EXTRA_PROGRAMS = lunatick

lunatick_SOURCES = main.cpp
lunatick_LDADD = -lncurses ../libnix.la
lunatick_LDFLAGS = -framework SDL2 -framework SDL2_image -framework SDL2_ttf
#lunatick_CXXFLAGS = -msse --version