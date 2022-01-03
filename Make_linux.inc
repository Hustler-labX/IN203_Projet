CC = gcc
CXX = g++-10
LIB = -fopenmp
CXXFLAGS+=$(mpic++ --showme:compile)
LIBS+=$(mpic++ --showme:link)
	
