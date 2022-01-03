CXXFLAGS+=$(mpic++ --showme:compile)
LIBS+=$(mpic++ --showme:link)
	
CC = mpicc
CXX = mpic++
LIB = 	-lmpi
MPI= -oversubscribe

