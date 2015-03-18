AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = offlattice #offlattice-gui

offlattice_SOURCES = main.cpp Agent.cpp Kinetics.cpp Covariance.cpp
offlattice_LDADD = 
offlattice_LDFLAGS =
offlattice_CXXFLAGS = -I.. -DNOGUI -g

offlattice_gui_SOURCES = main.cpp Agent.cpp Kinetics.cpp Color.cpp Covariance.cpp Window.cpp moc_Window.cpp moc_WindowGL.cpp
offlattice_gui_LDADD = 
offlattice_gui_LDFLAGS = -F/Users/jagiella/Qt/5.4/clang_64/lib -framework GLUT -framework QtCore -framework OpenGL -framework QtOpenGL -framework QtWidgets -framework QtGui 
offlattice_gui_CXXFLAGS = -I.. -I/Users/jagiella/Qt/5.4/clang_64/include/
