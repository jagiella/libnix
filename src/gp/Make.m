mex  src/gplikelihood.cpp           -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas -DMATLAB
mex  src/gp.cpp                     -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas -DMATLAB
mex  src/gpevaluation.cpp           -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas -DMATLAB
mex  src/gptraining.cpp src/cov.cpp -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas -DMATLAB
mex  src/covMatern5.cpp src/cov.cpp -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas -DMATLAB
mex  src/matrix.cpp                 -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas -DMATLAB

movefile( '*.mexmaci64', 'matlab')
