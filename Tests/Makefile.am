.PHONY:

NPARTICLES=10

# TESTS
############################
./Data/limits.dat: ../src/testlimits
	OMP_NUM_THREADS=8 ./$< $(NPARTICLES) 1 0 200 $@

./Data/limits_old.dat: ../src/testlimits_old
	OMP_NUM_THREADS=4 ./$< $(NPARTICLES) runge_kutta_fehlberg78 1 0 200 30 1 25 $@ 1

./Data/limits.npz: ./Data/limits.dat ../PostProcessing/processData.py
	python ../PostProcessing/processData.py Data/limits.dat ./Data/limits

./Plots/Density.dat: ./Data/limits.npz ../PostProcessing/calculateDensity.py
	python ../PostProcessing/calculateDensity.py ./Data/limits $@ 1 200

./Plots/Density.pdf: ./Plots/Density.dat ../PostProcessing/plotDensity.py
	python ../PostProcessing/plotDensity.py $< $@ log

./Plots/Comparison.pdf: ./Data/limits_old.dat ./Data/limits.dat ./comparison.py ../Paris/Analysis/GR.py
	python ./comparison.py ./Data/limits.dat ./Data/limits_old.dat $@ log

./Plots/ComparisonNormal.pdf: ./Data/limits_old.dat ./Data/limits.dat ./comparison.py ../Paris/Analysis/GR.py
	python ./comparison.py ./Data/limits.dat ./Data/limits_old.dat $@ normal

./Plots/ComparisonProduction.pdf: ./Data/limits_old.dat ./Data/limits.dat ./comparison.py ../Paris/Analysis/GR.py
	python ./comparison.py ./Data/limits.dat ./Data/limits_old.dat $@ production

# Test COM Energy computation
Data/limitsCOM.dat: calculateCOM.py Data/limits.npz
	python $< Data/limits $@

Plots/limitsCOM.pdf: plotCOM.py Data/limitsCOM.dat
	python $< Data/limitsCOM.dat $@ normal

testsdata: ./Data/limits.dat ./Data/limits_old.dat ./Plots/Density.dat ./Data/limits.npz Data/limitsCOM.dat
testsplot: ./Plots/Density.pdf ./Plots/Comparison.pdf ./Plots/ComparisonNormal.pdf ./Plots/ComparisonProduction.pdf Plots/limitsCOM.pdf
############################

# TESTS (spin = 1)
############################
./Data/limitsspin.dat: ../src/testlimits
	OMP_NUM_THREADS=4 ./$< $(NPARTICLES) 1 0.998 160 $@

./Data/limits_oldspin.dat: ../src/testlimits_old
	OMP_NUM_THREADS=4 ./$< $(NPARTICLES) runge_kutta_fehlberg78 1 0.998 200 30 1 25 $@ 1

./Data/limitsspin.npz: ./Data/limitsspin.dat ../PostProcessing/processData.py
	python ../PostProcessing/processData.py Data/limitsspin.dat ./Data/limitsspin

./Plots/Densityspin.dat: ./Data/limitsspin.npz ../PostProcessing/calculateDensity.py
	python ../PostProcessing/calculateDensity.py ./Data/limitsspin $@ 1 200

./Plots/Densityspin.pdf: ./Plots/Densityspin.dat ../PostProcessing/plotDensity.py
	python ../PostProcessing/plotDensity.py $< $@ log

./Plots/Comparisonspin.pdf: ./Data/limits_oldspin.dat ./Data/limitsspin.dat ./comparison.py ../Paris/Analysis/GR.py
	python ./comparison.py ./Data/limitsspin.dat ./Data/limits_oldspin.dat $@ log

./Plots/ComparisonNormalspin.pdf: ./Data/limits_oldspin.dat ./Data/limitsspin.dat ./comparison.py ../Paris/Analysis/GR.py
	python ./comparison.py ./Data/limitsspin.dat ./Data/limits_oldspin.dat $@ normal

./Plots/ComparisonProductionspin.pdf: ./Data/limits_oldspin.dat ./Data/limitsspin.dat ./comparison.py ../Paris/Analysis/GR.py
	python ./comparison.py ./Data/limitsspin.dat ./Data/limits_oldspin.dat $@ production

testsdataspin: ./Data/limitsspin.dat ./Data/limits_oldspin.dat ./Plots/Densityspin.dat ./Data/limitsspin.npz
testsplotspin: ./Plots/Densityspin.pdf ./Plots/Comparisonspin.pdf ./Plots/ComparisonNormalspin.pdf ./Plots/ComparisonProductionspin.pdf
############################


# COM Tests with kerr black hole
############################
Data/COMKerrGrowth.dat: ../src/growthsampling
	OMP_NUM_THREADS=8 $< $(NPARTICLES) 1 500 0.998 3 $@ 0 0 0 0.0006671 2

Data/COMKerrGrowth.npz: Data/COMKerrGrowth.dat ../PostProcessing/processData.py
	python ../PostProcessing/processData.py $< Data/COMKerrGrowth

Data/COMKerrGrowth.radial: calculateCOM.py Data/COMKerrGrowth.npz
	python $< Data/COMKerrGrowth $@

Data/COMKerrGrowth.density: ../PostProcessing/calculateDensity.py Data/COMKerrGrowth.npz
	python $< Data/COMKerrGrowth $@

Plots/COMKerrGrowth.pdf: plotCOM.py Data/COMKerrGrowth.radial
	python $< Data/COMKerrGrowth.radial $@ normal

Plots/COMKerrGrowthDensity.pdf: ../PostProcessing/plotDensity.py Data/COMKerrGrowth.density
	python $< Data/COMKerrGrowth.density $@ normal

Data/COMKerrGrowthMax.radial: calculateCOM.py Data/COMKerrGrowth.npz
	python $< Data/COMKerrGrowth $@ 1

Plots/COMKerrGrowthMax.pdf: plotCOM.py Data/COMKerrGrowthMax.radial
	python $< Data/COMKerrGrowthMax.radial $@ normal



COMdata: Data/COMKerrGrowth.dat Data/COMKerrGrowth.npz Data/COMKerrGrowth.radial Data/COMKerrGrowth.density
COMPlots: Plots/COMKerrGrowth.pdf Plots/COMKerrGrowthDensity.pdf
############################

# Data for testconstantcore:
############################
if USE_GYOTO
Data/limits.processed: processConstantCoreData.py ./Data/limits.npz
	python $< Data/limits $@

Data/grid.dat: ../src/testconstantcore Data/limits.processed
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 Data/grid 1 Data/limits.processed

Plots/comparisonPlot0.pdf: createComparisonPlot.py Data/grid.dat Data/limits.npz
	python $< Data/grid.dat Data/limits $@


Data/limits1.dat: ../src/testlimits
	OMP_NUM_THREADS=8 ./$< $(NPARTICLES) 1 0.998 200 $@

Data/limits1.npz: ./Data/limits1.dat ../PostProcessing/processData.py
	python ../PostProcessing/processData.py Data/limits1.dat Data/limits1

Data/limits1.processed: processConstantCoreData.py ./Data/limits1.npz
	python $< Data/limits1 $@

Data/grid1.dat: ../src/testconstantcore Data/limits1.processed
	OMP_NUM_THREADS=40 $< runge_kutta_fehlberg78 200 30 1 200 Data/grid1 1 Data/limits1.processed

Plots/comparisonPlot1.pdf: createComparisonPlot.py Data/grid1.dat Data/limits1.npz
	python $< Data/grid1.dat Data/limits1 $@

testconstantcoredata: Data/limits.processed Data/grid.dat Data/limits1.dat Data/limits1.npz Data/limits1.processed Data/grid1.dat
testconstantcoreplot: Plots/comparisonPlot0.pdf Plots/comparisonPlot1.pdf
testconstantcore: testconstantcoredata testconstantcoreplot
else
testconstantcoredata: 
testconstantcoreplot: 
testconstantcore: 
endif
############################

# Final test with N(E,L,Q) sampled via the growthsampling procedure
# - Outputs first N(E,L,Q) and E',L',Q' into distribution.dat
# - We process this to create f(E,L,Q)=f(E',L',Q') AND N'(E',L',Q') and output results into distribution.npz
# - We then input the results into Gyoto and output the density into degeneracygrid.dat
# - We create a plot both using f(E',L',Q')->rho and rho from Gyoto and plot them both into Plots/finalComparison.pdf
############################
if USE_GYOTO
Data/degeneracydistribution.dat: ../src/growthsampling
	OMP_NUM_THREADS=40 $< $(NPARTICLES) 1 500 0.998 15 $@ 0 0 0 0.0006671 2

Data/degeneracydistribution.npz: ../PostProcessing/processData.py Data/degeneracydistribution.dat
	python $< Data/degeneracydistribution.dat Data/degeneracydistribution

Data/degeneracydistributiongyoto.dat: convertDataToGyoto.py Data/degeneracydistribution.npz
	python $< Data/degeneracydistribution $@


Data/degeneracygrid.dat: ../src/rungyoto Data/degeneracydistributiongyoto.dat
	OMP_NUM_THREADS=1 $< runge_kutta_fehlberg78 200 30 1 15 $@ 1 Data/degeneracydistributiongyoto.dat

Plots/finalcomparisonbetweengyotoandanalytical.pdf: createFinalComparison.py Data/degeneracygrid.dat Data/degeneracydistribution.npz
	python $< Data/degeneracygrid.dat Data/degeneracydistribution $@

# Intermediate results
###
Data/degeneracyrho.dat: ../PostProcessing/calculateDensity.py Data/degeneracydistribution.npz
	python $< Data/degeneracydistribution $@

Plots/degeneracyrho.pdf: ../PostProcessing/plotDensity.py Data/degeneracyrho.dat
	python $< Data/degeneracyrho.dat $@ normal

Plots/degeneracylogrho.pdf: ../PostProcessing/plotDensity.py Data/degeneracyrho.dat
	python $< Data/degeneracyrho.dat $@ log
###


finaltestdata: Data/degeneracyrho.dat Data/degeneracydistribution.dat Data/degeneracydistribution.npz Data/degeneracydistributiongyoto.dat Data/degeneracygrid.dat
finaltestplot: Plots/degeneracyrho.pdf Plots/degeneracylogrho.pdf Plots/finalcomparisonbetweengyotoandanalytical.pdf
finaltest: finaltestdata finaltestplot testsdata
else
finaltestdata:
finaltestplot:
finaltest:
endif
############################


# Sampling test
############################
Data/testgrowthsampling.dat: ../src/testgrowthsampling
	OMP_NUM_THREADS=1 $< 20 1 500 0 200 $@ 0 0 0 0.0006671 2

# Records whether we were successful in calculating the constants of motions after growth
Data/testgrowthsampling.success: testadiabaticgrowth.py Data/testgrowthsampling.dat
	python $< Data/testgrowthsampling.dat $@

testadiabaticgrowthsolverdata: Data/testgrowthsampling.dat Data/testgrowthsampling.success
############################


# Test lookup table
############################
LOOKUPMASS=500
Data/lookup_0.998.dat: ../src/testLookupTableMC
	OMP_NUM_THREADS=20 $< 10000000 $(LOOKUPMASS) 0.998 200 $@

Data/lookuptabledata.dat: ../src/testLookupTableRead Data/lookup_0.998.dat
	OMP_NUM_THREADS=1 $< 1000 1 $(LOOKUPMASS) 0.998 150 $@ 0 0 0 0.0006671 2 Data/lookup_0.998.dat

# Process data
Data/lookuptabledata.npz: ../PostProcessing/processData.py Data/lookuptabledata.dat
	python $< Data/lookuptabledata.dat Data/lookuptabledata

# Create density plot
Data/lookuprho.dat: ../PostProcessing/calculateDensity.py Data/lookuptabledata.npz
	python $< Data/lookuptabledata $@

Plots/lookuprho.pdf: ../PostProcessing/plotDensity.py Data/lookuprho.dat
	python $< Data/lookuprho.dat $@ normal

Plots/lookuplogrho.pdf: ../PostProcessing/plotDensity.py Data/lookuprho.dat
	python $< Data/lookuprho.dat $@ log

lookupdata: Data/lookup_0.998.dat Data/lookuptabledata.dat Data/lookuptabledata.npz Data/lookuprho.dat
lookupplots: Plots/lookuprho.pdf Plots/lookuplogrho.pdf
lookup: lookupdata lookupplots
############################

# Create plot of FFIO
############################
Plots/FFIOVolume.pdf: plotFFIOVolume.py
	python $< $@
############################

data:  testconstantcoredata testsdata testadiabaticgrowthsolverdata finaltestdata testsdataspin #lookupdata
plots: testconstantcoreplot testsplot finaltestplot testsplotspin #lookupplots
all: data plots
clean:
	rm Data/*dat Plots/*pdf
