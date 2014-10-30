AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = nix-tumor3d nix-tumor2d nix-compare3d nix-compare2d

nix_tumor3d_SOURCES = main.cpp Montecarlo.cpp Agent.cpp Angiogenesis.cpp CellProcesses.cpp DiffusionReactionEquation.cpp finiteDifferences.cpp Interpolation.cpp Mathematix.cpp Molecule.cpp SparseMatrix.cpp Substrate.cpp EpsIO.cpp PovrayIO.cpp Vessel_Graph.cpp Vessel_Unit.cpp VesselNetwork.cpp VoronoiDiagramExtended.cpp
nix_tumor3d_LDADD = -lm 
nix_tumor3d_CXXFLAGS = -DDIMENSIONS=3

nix_tumor2d_SOURCES = main.cpp Montecarlo.cpp Agent.cpp Angiogenesis.cpp CellProcesses.cpp DiffusionReactionEquation.cpp finiteDifferences.cpp Interpolation.cpp Mathematix.cpp Molecule.cpp SparseMatrix.cpp Substrate.cpp EpsIO.cpp PovrayIO.cpp Vessel_Graph.cpp Vessel_Unit.cpp VesselNetwork.cpp VoronoiDiagramExtended.cpp
nix_tumor2d_LDADD = -lm 
nix_tumor2d_CXXFLAGS = -DDIMENSIONS=2

nix_compare3d_SOURCES = compare.cpp Montecarlo.cpp Agent.cpp Angiogenesis.cpp CellProcesses.cpp DiffusionReactionEquation.cpp finiteDifferences.cpp Interpolation.cpp Mathematix.cpp Molecule.cpp SparseMatrix.cpp Substrate.cpp EpsIO.cpp PovrayIO.cpp Vessel_Graph.cpp Vessel_Unit.cpp VesselNetwork.cpp VoronoiDiagramExtended.cpp
nix_compare3d_LDADD = -lm ../libnix.la
nix_compare3d_CXXFLAGS = -DDIMENSIONS=3

nix_compare2d_SOURCES = compare.cpp Montecarlo.cpp Agent.cpp Angiogenesis.cpp CellProcesses.cpp DiffusionReactionEquation.cpp finiteDifferences.cpp Interpolation.cpp Mathematix.cpp Molecule.cpp SparseMatrix.cpp Substrate.cpp EpsIO.cpp PovrayIO.cpp Vessel_Graph.cpp Vessel_Unit.cpp VesselNetwork.cpp VoronoiDiagramExtended.cpp
nix_compare2d_LDADD = -lm ../libnix.la
nix_compare2d_CXXFLAGS = -DDIMENSIONS=2
