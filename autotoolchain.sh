echo "test" > NEWS
echo "test" > README
echo "test" > AUTHORS
echo "test" > ChangeLog

autoscan;
sed -e '7a\'$'\n''AM_INIT_AUTOMAKE' configure.scan > configure.scan1
sed -e '8a\'$'\n''CC=${CC-/opt/local/bin/g++-mp-4.9}' configure.scan1 > configure.scan2
sed -e '9a\'$'\n''CPPFLAGS="$CPPFLAGS -I/opt/local/lib/gcc49/gcc/x86_64-apple-darwin14/4.9.1/include -I/opt/local/include"' configure.scan2 > configure.scan1
sed -e '9a\'$'\n''LDFLAGS="$LDFLAGS -L/opt/local/lib -L/opt/local/lib/gcc49"' configure.scan1 > configure.ac
rm configure.scan?
aclocal
autoheader
autoconf

automake --add-missing 
