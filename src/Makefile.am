#SUBDIRS = cubature
AM_CXXFLAGS = -std=c++11 -fopenmp
#AM_LIBADD= -gsl -gslcblas -lopenblas -lgyoto -lgyoto-stdplug -lxerces-c


bin_PROGRAMS=growthsampling testlimits testlimits_old testgrowthsampling testLookupTable testLookupTableMC testLookupTableRead 

if USE_GYOTO
bin_PROGRAMS+=rungyoto testconstantcore
testconstantcore_SOURCES=testconstantcore.cpp grid.cpp helperFunctions.cpp misc.cpp rpoly_ak1.cpp sampling.cpp cubature/hcubature.c
testconstantcore_LDADD= -lgyoto -lgyoto-stdplug -lxerces-c
rungyoto_SOURCES= rungyoto.cpp grid.cpp helperFunctions.cpp misc.cpp rpoly_ak1.cpp sampling.cpp cubature/hcubature.c
rungyoto_LDADD= -lgyoto -lgyoto-stdplug -lxerces-c
endif

growthsampling_SOURCES=helperFunctions.cpp misc.cpp growthsampling.cpp rpoly_ak1.cpp sampling.cpp cubature/hcubature.c haloFunctions.cpp
testgrowthsampling_SOURCES=helperFunctions.cpp misc.cpp testgrowthsampling.cpp rpoly_ak1.cpp sampling.cpp cubature/hcubature.c haloFunctions.cpp
testlimits_SOURCES=helperFunctions.cpp misc.cpp rpoly_ak1.cpp sampling.cpp cubature/hcubature.c testlimits.cpp
#testlimits_LDADD=testlimits.o
testlimits_old_SOURCES= helperFunctions.cpp misc.cpp testlimits_old.cpp rpoly_ak1.cpp sampling.cpp cubature/hcubature.c
testLookupTable_SOURCES = testLookupTable.cpp helperFunctions.cpp rpoly_ak1.cpp sampling.cpp cubature/hcubature.c
testLookupTableMC_SOURCES = testLookupTableMC.cpp helperFunctions.cpp rpoly_ak1.cpp sampling.cpp cubature/hcubature.c
testLookupTableRead_SOURCES = testLookupTableRead.cpp helperFunctions.cpp rpoly_ak1.cpp sampling.cpp cubature/hcubature.c haloFunctions.cpp misc.cpp lookupTable.cpp GrowBlackHole.cpp

