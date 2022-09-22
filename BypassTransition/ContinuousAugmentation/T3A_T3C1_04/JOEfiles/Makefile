default: joe_1 joe_2

joe_1: joe.cpp
	mpicxx $(JOE_COMPILE_PATH) joe.cpp $(JOE_COMPILE_LIBS) -o joe_direct

joe_2: joe.cpp
	mpicxx -DADJOINT_MODE $(JOE_COMPILE_PATH) joe.cpp $(JOE_COMPILE_LIBS) -o joe_adjoint
