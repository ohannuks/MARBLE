.PHONY : testplot productionfigures

data: grid.dat smallgrid.dat biggrid.dat spin0.dat
datap: grid.npz smallgrid.npz biggrid.npz
plotdata: Plots/plotgrid.dat Plots/plotsmallgrid.dat Plots/plotspin0.dat
plots: Plots/plotgrid.pdf Plots/plotsmallgrid.pdf Plots/plotspin0.pdf

#TODO: Change these to be from e.g. 1 to 400, 500, etc
############################
grid.dat: ..//src/growthsampling
	OMP_NUM_THREADS=4 ./$^ 10000 1 500 1 150 $@ 0 0 0 0.0006671 2

smallgrid.dat: ..//src/growthsampling
	OMP_NUM_THREADS=4 ./$^ 4000 1 500 1 150 $@ 0 0 0 0.0006671 2

biggrid.dat: ..//src/growthsampling
	OMP_NUM_THREADS=4 ./$^ 80000 1 500 1 150 $@ 0 0 0 0.0006671 2

# This does not seem to work
#spin0.dat: ..//src/growthsampling
#	OMP_NUM_THREADS=4 ./$^ 40000 1 400 0 150 $@ 0 0 0 0.0006671 2
############################

# Process data
############################
grid.npz: grid.dat processData.py
	python processData.py $< grid

smallgrid.npz: smallgrid.dat processData.py
	python processData.py $< grid

biggrid.npz: biggrid.dat processData.py
	python processData.py $< biggrid

spin0.npz: spin0.dat processData.py
	python processData.py $< spin0
############################

# Plots
############################
Plots/plotgrid.dat: grid.dat ./calculateDensity.py
	python ./calculateDensity.py $< $@ 1 20

Plots/plotgrid.pdf: Plots/plotgrid.dat ./plotDensity.py
	python ./plotDensity.py $< $@

Plots/plotsmallgrid.dat: smallgrid.dat ./plotDensity.py
	python ./calculateDensity.py $< $@ 1 20

Plots/plotsmallgrid.pdf: Plots/plotsmallgrid.dat ./plotDensity.py
	python ./plotDensity.py $< $@

Plots/plotspin0.dat: spin0.dat ./plotDensity.py
	python ./calculateDensity.py $< $@ 1 20

Plots/plotspin0.pdf: Plots/plotspin0.dat ./plotDensity.py
	python ./plotDensity.py $< $@
############################

Plots/FFIOZAMODensityComparison.pdf: FFIOZAMOComparison.py ProductionData/M10_r2_spin0998.radial
	python $< ProductionData/M10_r2_spin0998.radial $@

# Calculate COM
############################
Data/COMKerrGrowthSpin1.radial: ../Tests/calculateCOM.py ProductionData/M10_r2_spin0998.npz
	python $< ProductionData/M10_r2_spin0998 $@

Data/COMKerrGrowthSpin0.radial: ../Tests/calculateCOM.py ProductionData/M10_r2_spin0000.npz
	python $< ProductionData/M10_r2_spin0000 $@

Plots/COMKerrGrowthSpin1.pdf: ../Tests/plotCOM.py Data/COMKerrGrowthSpin1.radial
	python $< Data/COMKerrGrowthSpin1.radial $@ normal

Plots/COMKerrGrowthSpin0.pdf: ../Tests/plotCOM.py Data/COMKerrGrowthSpin0.radial
	python $< Data/COMKerrGrowthSpin0.radial $@ normal

Data/COMKerrGrowthMaxSpin1.radial: ../Tests/calculateCOM.py ProductionData/M10_r2_spin0998.npz
	python $< ProductionData/M10_r2_spin0998 $@ 1

Plots/COMKerrGrowthMaxSpin1.pdf: ../Tests/plotCOM.py Data/COMKerrGrowthMaxSpin1.radial
	python $< Data/COMKerrGrowthMaxSpin1.radial $@ normal

Data/COMKerrGrowthMaxSpin0.radial: ../Tests/calculateCOM.py ProductionData/M10_r2_spin0000.npz
	python $< ProductionData/M10_r2_spin0000 $@ 1

Plots/COMKerrGrowthMaxSpin0.pdf: ../Tests/plotCOM.py Data/COMKerrGrowthMaxSpin0.radial
	python $< Data/COMKerrGrowthMaxSpin0.radial $@ normal

Plots/COMKerrGrowthMaxComparison.pdf: ../Tests/plotCOMComparison.py Data/COMKerrGrowthMaxSpin0.radial Data/COMKerrGrowthMaxSpin1.radial
	python $< Data/COMKerrGrowthMaxSpin0.radial Data/COMKerrGrowthMaxSpin1.radial $@ normal

Plots/COMKerrGrowthMaxComparisonLog.pdf: ../Tests/plotCOMComparison.py Data/COMKerrGrowthMaxSpin0.radial Data/COMKerrGrowthMaxSpin1.radial
	python $< Data/COMKerrGrowthMaxSpin0.radial Data/COMKerrGrowthMaxSpin1.radial $@ log


comdata: Data/COMKerrGrowthSpin1.radial Data/COMKerrGrowthSpin0.radial Data/COMKerrGrowthMaxSpin1.radial Data/COMKerrGrowthMaxSpin0.radial
complots: Plots/COMKerrGrowthSpin1.pdf Plots/COMKerrGrowthSpin0.pdf Plots/COMKerrGrowthMaxSpin1.pdf Plots/COMKerrGrowthMaxSpin0.pdf Plots/COMKerrGrowthMaxComparison.pdf Plots/COMKerrGrowthMaxComparisonLog.pdf 
com: comdata complots
############################

# Development runs
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
ProductionData/M9_r1_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 0.900000 500 0 50 $@ 0 0 0 0.0006671 1 
ProductionData/M9_r1_spin0000.npz: processData.py ProductionData/M9_r1_spin0000.dat 
	python $< ProductionData/M9_r1_spin0000.dat ProductionData/M9_r1_spin0000 
ProductionData/M9_r1_spin0000.radial: ProductionDataVariables.py ProductionData/M9_r1_spin0000.npz 
	python $< ProductionData/M9_r1_spin0000 ProductionData/M9_r1_spin0000.radial 
# Intermediate results 
ProductionData/M9_r1_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M9_r1_spin0000.npz 
	python $< ProductionData/M9_r1_spin0000 $@ 
Plots/M9_r1_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M9_r1_spin0000rho.dat 
	python $< ProductionData/M9_r1_spin0000rho.dat $@ normal 
Plots/M9_r1_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M9_r1_spin0000rho.dat 
	python $< ProductionData/M9_r1_spin0000rho.dat $@ log 
ProductionData/M9_r1_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M9_r1_spin0000.npz 
	python $< ProductionData/M9_r1_spin0000 $@ 
ProductionData/M9_r1_spin0000grid.dat: ../src/rungyoto ProductionData/M9_r1_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M9_r1_spin0000Gyoto.dat 

ProductionData/M9_r1_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 0.900000 500 0.998 50 $@ 0 0 0 0.0006671 1 
ProductionData/M9_r1_spin0998.npz: processData.py ProductionData/M9_r1_spin0998.dat 
	python $< ProductionData/M9_r1_spin0998.dat ProductionData/M9_r1_spin0998 
ProductionData/M9_r1_spin0998.radial: ProductionDataVariables.py ProductionData/M9_r1_spin0998.npz 
	python $< ProductionData/M9_r1_spin0998 ProductionData/M9_r1_spin0998.radial 
# Intermediate results 
ProductionData/M9_r1_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M9_r1_spin0998.npz 
	python $< ProductionData/M9_r1_spin0998 $@ 
Plots/M9_r1_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M9_r1_spin0998rho.dat 
	python $< ProductionData/M9_r1_spin0998rho.dat $@ normal 
Plots/M9_r1_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M9_r1_spin0998rho.dat 
	python $< ProductionData/M9_r1_spin0998rho.dat $@ log 
ProductionData/M9_r1_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M9_r1_spin0998.npz 
	python $< ProductionData/M9_r1_spin0998 $@ 
ProductionData/M9_r1_spin0998grid.dat: ../src/rungyoto ProductionData/M9_r1_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M9_r1_spin0998Gyoto.dat 

ProductionData/M9_r2_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 0.900000 500 0 50 $@ 0 0 0 0.0006671 2 
ProductionData/M9_r2_spin0000.npz: processData.py ProductionData/M9_r2_spin0000.dat 
	python $< ProductionData/M9_r2_spin0000.dat ProductionData/M9_r2_spin0000 
ProductionData/M9_r2_spin0000.radial: ProductionDataVariables.py ProductionData/M9_r2_spin0000.npz 
	python $< ProductionData/M9_r2_spin0000 ProductionData/M9_r2_spin0000.radial 
# Intermediate results 
ProductionData/M9_r2_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M9_r2_spin0000.npz 
	python $< ProductionData/M9_r2_spin0000 $@ 
Plots/M9_r2_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M9_r2_spin0000rho.dat 
	python $< ProductionData/M9_r2_spin0000rho.dat $@ normal 
Plots/M9_r2_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M9_r2_spin0000rho.dat 
	python $< ProductionData/M9_r2_spin0000rho.dat $@ log 
ProductionData/M9_r2_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M9_r2_spin0000.npz 
	python $< ProductionData/M9_r2_spin0000 $@ 
ProductionData/M9_r2_spin0000grid.dat: ../src/rungyoto ProductionData/M9_r2_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M9_r2_spin0000Gyoto.dat 

ProductionData/M9_r2_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 0.900000 500 0.998 50 $@ 0 0 0 0.0006671 2 
ProductionData/M9_r2_spin0998.npz: processData.py ProductionData/M9_r2_spin0998.dat 
	python $< ProductionData/M9_r2_spin0998.dat ProductionData/M9_r2_spin0998 
ProductionData/M9_r2_spin0998.radial: ProductionDataVariables.py ProductionData/M9_r2_spin0998.npz 
	python $< ProductionData/M9_r2_spin0998 ProductionData/M9_r2_spin0998.radial 
# Intermediate results 
ProductionData/M9_r2_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M9_r2_spin0998.npz 
	python $< ProductionData/M9_r2_spin0998 $@ 
Plots/M9_r2_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M9_r2_spin0998rho.dat 
	python $< ProductionData/M9_r2_spin0998rho.dat $@ normal 
Plots/M9_r2_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M9_r2_spin0998rho.dat 
	python $< ProductionData/M9_r2_spin0998rho.dat $@ log 
ProductionData/M9_r2_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M9_r2_spin0998.npz 
	python $< ProductionData/M9_r2_spin0998 $@ 
ProductionData/M9_r2_spin0998grid.dat: ../src/rungyoto ProductionData/M9_r2_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M9_r2_spin0998Gyoto.dat 

ProductionData/M10_r1_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.000000 500 0 50 $@ 0 0 0 0.0006671 1 
ProductionData/M10_r1_spin0000.npz: processData.py ProductionData/M10_r1_spin0000.dat 
	python $< ProductionData/M10_r1_spin0000.dat ProductionData/M10_r1_spin0000 
ProductionData/M10_r1_spin0000.radial: ProductionDataVariables.py ProductionData/M10_r1_spin0000.npz 
	python $< ProductionData/M10_r1_spin0000 ProductionData/M10_r1_spin0000.radial 
# Intermediate results 
ProductionData/M10_r1_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M10_r1_spin0000.npz 
	python $< ProductionData/M10_r1_spin0000 $@ 
Plots/M10_r1_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M10_r1_spin0000rho.dat 
	python $< ProductionData/M10_r1_spin0000rho.dat $@ normal 
Plots/M10_r1_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M10_r1_spin0000rho.dat 
	python $< ProductionData/M10_r1_spin0000rho.dat $@ log 
ProductionData/M10_r1_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M10_r1_spin0000.npz 
	python $< ProductionData/M10_r1_spin0000 $@ 
ProductionData/M10_r1_spin0000grid.dat: ../src/rungyoto ProductionData/M10_r1_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M10_r1_spin0000Gyoto.dat 

ProductionData/M10_r1_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.000000 500 0.998 50 $@ 0 0 0 0.0006671 1 
ProductionData/M10_r1_spin0998.npz: processData.py ProductionData/M10_r1_spin0998.dat 
	python $< ProductionData/M10_r1_spin0998.dat ProductionData/M10_r1_spin0998 
ProductionData/M10_r1_spin0998.radial: ProductionDataVariables.py ProductionData/M10_r1_spin0998.npz 
	python $< ProductionData/M10_r1_spin0998 ProductionData/M10_r1_spin0998.radial 
# Intermediate results 
ProductionData/M10_r1_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M10_r1_spin0998.npz 
	python $< ProductionData/M10_r1_spin0998 $@ 
Plots/M10_r1_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M10_r1_spin0998rho.dat 
	python $< ProductionData/M10_r1_spin0998rho.dat $@ normal 
Plots/M10_r1_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M10_r1_spin0998rho.dat 
	python $< ProductionData/M10_r1_spin0998rho.dat $@ log 
ProductionData/M10_r1_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M10_r1_spin0998.npz 
	python $< ProductionData/M10_r1_spin0998 $@ 
ProductionData/M10_r1_spin0998grid.dat: ../src/rungyoto ProductionData/M10_r1_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M10_r1_spin0998Gyoto.dat 

ProductionData/M10_r2_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.000000 500 0 50 $@ 0 0 0 0.0006671 2 
ProductionData/M10_r2_spin0000.npz: processData.py ProductionData/M10_r2_spin0000.dat 
	python $< ProductionData/M10_r2_spin0000.dat ProductionData/M10_r2_spin0000 
ProductionData/M10_r2_spin0000.radial: ProductionDataVariables.py ProductionData/M10_r2_spin0000.npz 
	python $< ProductionData/M10_r2_spin0000 ProductionData/M10_r2_spin0000.radial 
# Intermediate results 
ProductionData/M10_r2_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M10_r2_spin0000.npz 
	python $< ProductionData/M10_r2_spin0000 $@ 
Plots/M10_r2_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M10_r2_spin0000rho.dat 
	python $< ProductionData/M10_r2_spin0000rho.dat $@ normal 
Plots/M10_r2_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M10_r2_spin0000rho.dat 
	python $< ProductionData/M10_r2_spin0000rho.dat $@ log 
ProductionData/M10_r2_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M10_r2_spin0000.npz 
	python $< ProductionData/M10_r2_spin0000 $@ 
ProductionData/M10_r2_spin0000grid.dat: ../src/rungyoto ProductionData/M10_r2_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M10_r2_spin0000Gyoto.dat 

ProductionData/M10_r2_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.000000 500 0.998 50 $@ 0 0 0 0.0006671 2 
ProductionData/M10_r2_spin0998.npz: processData.py ProductionData/M10_r2_spin0998.dat 
	python $< ProductionData/M10_r2_spin0998.dat ProductionData/M10_r2_spin0998 
ProductionData/M10_r2_spin0998.radial: ProductionDataVariables.py ProductionData/M10_r2_spin0998.npz 
	python $< ProductionData/M10_r2_spin0998 ProductionData/M10_r2_spin0998.radial 
# Intermediate results 
ProductionData/M10_r2_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M10_r2_spin0998.npz 
	python $< ProductionData/M10_r2_spin0998 $@ 
Plots/M10_r2_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M10_r2_spin0998rho.dat 
	python $< ProductionData/M10_r2_spin0998rho.dat $@ normal 
Plots/M10_r2_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M10_r2_spin0998rho.dat 
	python $< ProductionData/M10_r2_spin0998rho.dat $@ log 
ProductionData/M10_r2_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M10_r2_spin0998.npz 
	python $< ProductionData/M10_r2_spin0998 $@ 
ProductionData/M10_r2_spin0998grid.dat: ../src/rungyoto ProductionData/M10_r2_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M10_r2_spin0998Gyoto.dat 

ProductionData/M11_r1_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.100000 500 0 50 $@ 0 0 0 0.0006671 1 
ProductionData/M11_r1_spin0000.npz: processData.py ProductionData/M11_r1_spin0000.dat 
	python $< ProductionData/M11_r1_spin0000.dat ProductionData/M11_r1_spin0000 
ProductionData/M11_r1_spin0000.radial: ProductionDataVariables.py ProductionData/M11_r1_spin0000.npz 
	python $< ProductionData/M11_r1_spin0000 ProductionData/M11_r1_spin0000.radial 
# Intermediate results 
ProductionData/M11_r1_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M11_r1_spin0000.npz 
	python $< ProductionData/M11_r1_spin0000 $@ 
Plots/M11_r1_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M11_r1_spin0000rho.dat 
	python $< ProductionData/M11_r1_spin0000rho.dat $@ normal 
Plots/M11_r1_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M11_r1_spin0000rho.dat 
	python $< ProductionData/M11_r1_spin0000rho.dat $@ log 
ProductionData/M11_r1_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M11_r1_spin0000.npz 
	python $< ProductionData/M11_r1_spin0000 $@ 
ProductionData/M11_r1_spin0000grid.dat: ../src/rungyoto ProductionData/M11_r1_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M11_r1_spin0000Gyoto.dat 

ProductionData/M11_r1_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.100000 500 0.998 50 $@ 0 0 0 0.0006671 1 
ProductionData/M11_r1_spin0998.npz: processData.py ProductionData/M11_r1_spin0998.dat 
	python $< ProductionData/M11_r1_spin0998.dat ProductionData/M11_r1_spin0998 
ProductionData/M11_r1_spin0998.radial: ProductionDataVariables.py ProductionData/M11_r1_spin0998.npz 
	python $< ProductionData/M11_r1_spin0998 ProductionData/M11_r1_spin0998.radial 
# Intermediate results 
ProductionData/M11_r1_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M11_r1_spin0998.npz 
	python $< ProductionData/M11_r1_spin0998 $@ 
Plots/M11_r1_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M11_r1_spin0998rho.dat 
	python $< ProductionData/M11_r1_spin0998rho.dat $@ normal 
Plots/M11_r1_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M11_r1_spin0998rho.dat 
	python $< ProductionData/M11_r1_spin0998rho.dat $@ log 
ProductionData/M11_r1_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M11_r1_spin0998.npz 
	python $< ProductionData/M11_r1_spin0998 $@ 
ProductionData/M11_r1_spin0998grid.dat: ../src/rungyoto ProductionData/M11_r1_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M11_r1_spin0998Gyoto.dat 

ProductionData/M11_r2_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.100000 500 0 50 $@ 0 0 0 0.0006671 2 
ProductionData/M11_r2_spin0000.npz: processData.py ProductionData/M11_r2_spin0000.dat 
	python $< ProductionData/M11_r2_spin0000.dat ProductionData/M11_r2_spin0000 
ProductionData/M11_r2_spin0000.radial: ProductionDataVariables.py ProductionData/M11_r2_spin0000.npz 
	python $< ProductionData/M11_r2_spin0000 ProductionData/M11_r2_spin0000.radial 
# Intermediate results 
ProductionData/M11_r2_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M11_r2_spin0000.npz 
	python $< ProductionData/M11_r2_spin0000 $@ 
Plots/M11_r2_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M11_r2_spin0000rho.dat 
	python $< ProductionData/M11_r2_spin0000rho.dat $@ normal 
Plots/M11_r2_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M11_r2_spin0000rho.dat 
	python $< ProductionData/M11_r2_spin0000rho.dat $@ log 
ProductionData/M11_r2_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M11_r2_spin0000.npz 
	python $< ProductionData/M11_r2_spin0000 $@ 
ProductionData/M11_r2_spin0000grid.dat: ../src/rungyoto ProductionData/M11_r2_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M11_r2_spin0000Gyoto.dat 

ProductionData/M11_r2_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.100000 500 0.998 50 $@ 0 0 0 0.0006671 2 
ProductionData/M11_r2_spin0998.npz: processData.py ProductionData/M11_r2_spin0998.dat 
	python $< ProductionData/M11_r2_spin0998.dat ProductionData/M11_r2_spin0998 
ProductionData/M11_r2_spin0998.radial: ProductionDataVariables.py ProductionData/M11_r2_spin0998.npz 
	python $< ProductionData/M11_r2_spin0998 ProductionData/M11_r2_spin0998.radial 
# Intermediate results 
ProductionData/M11_r2_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M11_r2_spin0998.npz 
	python $< ProductionData/M11_r2_spin0998 $@ 
Plots/M11_r2_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M11_r2_spin0998rho.dat 
	python $< ProductionData/M11_r2_spin0998rho.dat $@ normal 
Plots/M11_r2_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M11_r2_spin0998rho.dat 
	python $< ProductionData/M11_r2_spin0998rho.dat $@ log 
ProductionData/M11_r2_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M11_r2_spin0998.npz 
	python $< ProductionData/M11_r2_spin0998 $@ 
ProductionData/M11_r2_spin0998grid.dat: ../src/rungyoto ProductionData/M11_r2_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M11_r2_spin0998Gyoto.dat 

ProductionData/M12_r1_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.200000 500 0 50 $@ 0 0 0 0.0006671 1 
ProductionData/M12_r1_spin0000.npz: processData.py ProductionData/M12_r1_spin0000.dat 
	python $< ProductionData/M12_r1_spin0000.dat ProductionData/M12_r1_spin0000 
ProductionData/M12_r1_spin0000.radial: ProductionDataVariables.py ProductionData/M12_r1_spin0000.npz 
	python $< ProductionData/M12_r1_spin0000 ProductionData/M12_r1_spin0000.radial 
# Intermediate results 
ProductionData/M12_r1_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M12_r1_spin0000.npz 
	python $< ProductionData/M12_r1_spin0000 $@ 
Plots/M12_r1_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M12_r1_spin0000rho.dat 
	python $< ProductionData/M12_r1_spin0000rho.dat $@ normal 
Plots/M12_r1_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M12_r1_spin0000rho.dat 
	python $< ProductionData/M12_r1_spin0000rho.dat $@ log 
ProductionData/M12_r1_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M12_r1_spin0000.npz 
	python $< ProductionData/M12_r1_spin0000 $@ 
ProductionData/M12_r1_spin0000grid.dat: ../src/rungyoto ProductionData/M12_r1_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M12_r1_spin0000Gyoto.dat 

ProductionData/M12_r1_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.200000 500 0.998 50 $@ 0 0 0 0.0006671 1 
ProductionData/M12_r1_spin0998.npz: processData.py ProductionData/M12_r1_spin0998.dat 
	python $< ProductionData/M12_r1_spin0998.dat ProductionData/M12_r1_spin0998 
ProductionData/M12_r1_spin0998.radial: ProductionDataVariables.py ProductionData/M12_r1_spin0998.npz 
	python $< ProductionData/M12_r1_spin0998 ProductionData/M12_r1_spin0998.radial 
# Intermediate results 
ProductionData/M12_r1_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M12_r1_spin0998.npz 
	python $< ProductionData/M12_r1_spin0998 $@ 
Plots/M12_r1_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M12_r1_spin0998rho.dat 
	python $< ProductionData/M12_r1_spin0998rho.dat $@ normal 
Plots/M12_r1_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M12_r1_spin0998rho.dat 
	python $< ProductionData/M12_r1_spin0998rho.dat $@ log 
ProductionData/M12_r1_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M12_r1_spin0998.npz 
	python $< ProductionData/M12_r1_spin0998 $@ 
ProductionData/M12_r1_spin0998grid.dat: ../src/rungyoto ProductionData/M12_r1_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M12_r1_spin0998Gyoto.dat 

ProductionData/M12_r2_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.200000 500 0 50 $@ 0 0 0 0.0006671 2 
ProductionData/M12_r2_spin0000.npz: processData.py ProductionData/M12_r2_spin0000.dat 
	python $< ProductionData/M12_r2_spin0000.dat ProductionData/M12_r2_spin0000 
ProductionData/M12_r2_spin0000.radial: ProductionDataVariables.py ProductionData/M12_r2_spin0000.npz 
	python $< ProductionData/M12_r2_spin0000 ProductionData/M12_r2_spin0000.radial 
# Intermediate results 
ProductionData/M12_r2_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M12_r2_spin0000.npz 
	python $< ProductionData/M12_r2_spin0000 $@ 
Plots/M12_r2_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M12_r2_spin0000rho.dat 
	python $< ProductionData/M12_r2_spin0000rho.dat $@ normal 
Plots/M12_r2_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M12_r2_spin0000rho.dat 
	python $< ProductionData/M12_r2_spin0000rho.dat $@ log 
ProductionData/M12_r2_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M12_r2_spin0000.npz 
	python $< ProductionData/M12_r2_spin0000 $@ 
ProductionData/M12_r2_spin0000grid.dat: ../src/rungyoto ProductionData/M12_r2_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M12_r2_spin0000Gyoto.dat 

ProductionData/M12_r2_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.200000 500 0.998 50 $@ 0 0 0 0.0006671 2 
ProductionData/M12_r2_spin0998.npz: processData.py ProductionData/M12_r2_spin0998.dat 
	python $< ProductionData/M12_r2_spin0998.dat ProductionData/M12_r2_spin0998 
ProductionData/M12_r2_spin0998.radial: ProductionDataVariables.py ProductionData/M12_r2_spin0998.npz 
	python $< ProductionData/M12_r2_spin0998 ProductionData/M12_r2_spin0998.radial 
# Intermediate results 
ProductionData/M12_r2_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M12_r2_spin0998.npz 
	python $< ProductionData/M12_r2_spin0998 $@ 
Plots/M12_r2_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M12_r2_spin0998rho.dat 
	python $< ProductionData/M12_r2_spin0998rho.dat $@ normal 
Plots/M12_r2_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M12_r2_spin0998rho.dat 
	python $< ProductionData/M12_r2_spin0998rho.dat $@ log 
ProductionData/M12_r2_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M12_r2_spin0998.npz 
	python $< ProductionData/M12_r2_spin0998 $@ 
ProductionData/M12_r2_spin0998grid.dat: ../src/rungyoto ProductionData/M12_r2_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M12_r2_spin0998Gyoto.dat 

ProductionData/M13_r1_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.300000 500 0 50 $@ 0 0 0 0.0006671 1 
ProductionData/M13_r1_spin0000.npz: processData.py ProductionData/M13_r1_spin0000.dat 
	python $< ProductionData/M13_r1_spin0000.dat ProductionData/M13_r1_spin0000 
ProductionData/M13_r1_spin0000.radial: ProductionDataVariables.py ProductionData/M13_r1_spin0000.npz 
	python $< ProductionData/M13_r1_spin0000 ProductionData/M13_r1_spin0000.radial 
# Intermediate results 
ProductionData/M13_r1_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M13_r1_spin0000.npz 
	python $< ProductionData/M13_r1_spin0000 $@ 
Plots/M13_r1_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M13_r1_spin0000rho.dat 
	python $< ProductionData/M13_r1_spin0000rho.dat $@ normal 
Plots/M13_r1_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M13_r1_spin0000rho.dat 
	python $< ProductionData/M13_r1_spin0000rho.dat $@ log 
ProductionData/M13_r1_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M13_r1_spin0000.npz 
	python $< ProductionData/M13_r1_spin0000 $@ 
ProductionData/M13_r1_spin0000grid.dat: ../src/rungyoto ProductionData/M13_r1_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M13_r1_spin0000Gyoto.dat 

ProductionData/M13_r1_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.300000 500 0.998 50 $@ 0 0 0 0.0006671 1 
ProductionData/M13_r1_spin0998.npz: processData.py ProductionData/M13_r1_spin0998.dat 
	python $< ProductionData/M13_r1_spin0998.dat ProductionData/M13_r1_spin0998 
ProductionData/M13_r1_spin0998.radial: ProductionDataVariables.py ProductionData/M13_r1_spin0998.npz 
	python $< ProductionData/M13_r1_spin0998 ProductionData/M13_r1_spin0998.radial 
# Intermediate results 
ProductionData/M13_r1_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M13_r1_spin0998.npz 
	python $< ProductionData/M13_r1_spin0998 $@ 
Plots/M13_r1_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M13_r1_spin0998rho.dat 
	python $< ProductionData/M13_r1_spin0998rho.dat $@ normal 
Plots/M13_r1_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M13_r1_spin0998rho.dat 
	python $< ProductionData/M13_r1_spin0998rho.dat $@ log 
ProductionData/M13_r1_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M13_r1_spin0998.npz 
	python $< ProductionData/M13_r1_spin0998 $@ 
ProductionData/M13_r1_spin0998grid.dat: ../src/rungyoto ProductionData/M13_r1_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M13_r1_spin0998Gyoto.dat 

ProductionData/M13_r2_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.300000 500 0 50 $@ 0 0 0 0.0006671 2 
ProductionData/M13_r2_spin0000.npz: processData.py ProductionData/M13_r2_spin0000.dat 
	python $< ProductionData/M13_r2_spin0000.dat ProductionData/M13_r2_spin0000 
ProductionData/M13_r2_spin0000.radial: ProductionDataVariables.py ProductionData/M13_r2_spin0000.npz 
	python $< ProductionData/M13_r2_spin0000 ProductionData/M13_r2_spin0000.radial 
# Intermediate results 
ProductionData/M13_r2_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M13_r2_spin0000.npz 
	python $< ProductionData/M13_r2_spin0000 $@ 
Plots/M13_r2_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M13_r2_spin0000rho.dat 
	python $< ProductionData/M13_r2_spin0000rho.dat $@ normal 
Plots/M13_r2_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M13_r2_spin0000rho.dat 
	python $< ProductionData/M13_r2_spin0000rho.dat $@ log 
ProductionData/M13_r2_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M13_r2_spin0000.npz 
	python $< ProductionData/M13_r2_spin0000 $@ 
ProductionData/M13_r2_spin0000grid.dat: ../src/rungyoto ProductionData/M13_r2_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M13_r2_spin0000Gyoto.dat 

ProductionData/M13_r2_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.300000 500 0.998 50 $@ 0 0 0 0.0006671 2 
ProductionData/M13_r2_spin0998.npz: processData.py ProductionData/M13_r2_spin0998.dat 
	python $< ProductionData/M13_r2_spin0998.dat ProductionData/M13_r2_spin0998 
ProductionData/M13_r2_spin0998.radial: ProductionDataVariables.py ProductionData/M13_r2_spin0998.npz 
	python $< ProductionData/M13_r2_spin0998 ProductionData/M13_r2_spin0998.radial 
# Intermediate results 
ProductionData/M13_r2_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M13_r2_spin0998.npz 
	python $< ProductionData/M13_r2_spin0998 $@ 
Plots/M13_r2_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M13_r2_spin0998rho.dat 
	python $< ProductionData/M13_r2_spin0998rho.dat $@ normal 
Plots/M13_r2_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M13_r2_spin0998rho.dat 
	python $< ProductionData/M13_r2_spin0998rho.dat $@ log 
ProductionData/M13_r2_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M13_r2_spin0998.npz 
	python $< ProductionData/M13_r2_spin0998 $@ 
ProductionData/M13_r2_spin0998grid.dat: ../src/rungyoto ProductionData/M13_r2_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M13_r2_spin0998Gyoto.dat 

ProductionData/M14_r1_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.400000 500 0 50 $@ 0 0 0 0.0006671 1 
ProductionData/M14_r1_spin0000.npz: processData.py ProductionData/M14_r1_spin0000.dat 
	python $< ProductionData/M14_r1_spin0000.dat ProductionData/M14_r1_spin0000 
ProductionData/M14_r1_spin0000.radial: ProductionDataVariables.py ProductionData/M14_r1_spin0000.npz 
	python $< ProductionData/M14_r1_spin0000 ProductionData/M14_r1_spin0000.radial 
# Intermediate results 
ProductionData/M14_r1_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M14_r1_spin0000.npz 
	python $< ProductionData/M14_r1_spin0000 $@ 
Plots/M14_r1_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M14_r1_spin0000rho.dat 
	python $< ProductionData/M14_r1_spin0000rho.dat $@ normal 
Plots/M14_r1_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M14_r1_spin0000rho.dat 
	python $< ProductionData/M14_r1_spin0000rho.dat $@ log 
ProductionData/M14_r1_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M14_r1_spin0000.npz 
	python $< ProductionData/M14_r1_spin0000 $@ 
ProductionData/M14_r1_spin0000grid.dat: ../src/rungyoto ProductionData/M14_r1_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M14_r1_spin0000Gyoto.dat 

ProductionData/M14_r1_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.400000 500 0.998 50 $@ 0 0 0 0.0006671 1 
ProductionData/M14_r1_spin0998.npz: processData.py ProductionData/M14_r1_spin0998.dat 
	python $< ProductionData/M14_r1_spin0998.dat ProductionData/M14_r1_spin0998 
ProductionData/M14_r1_spin0998.radial: ProductionDataVariables.py ProductionData/M14_r1_spin0998.npz 
	python $< ProductionData/M14_r1_spin0998 ProductionData/M14_r1_spin0998.radial 
# Intermediate results 
ProductionData/M14_r1_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M14_r1_spin0998.npz 
	python $< ProductionData/M14_r1_spin0998 $@ 
Plots/M14_r1_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M14_r1_spin0998rho.dat 
	python $< ProductionData/M14_r1_spin0998rho.dat $@ normal 
Plots/M14_r1_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M14_r1_spin0998rho.dat 
	python $< ProductionData/M14_r1_spin0998rho.dat $@ log 
ProductionData/M14_r1_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M14_r1_spin0998.npz 
	python $< ProductionData/M14_r1_spin0998 $@ 
ProductionData/M14_r1_spin0998grid.dat: ../src/rungyoto ProductionData/M14_r1_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M14_r1_spin0998Gyoto.dat 

ProductionData/M14_r2_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.400000 500 0 50 $@ 0 0 0 0.0006671 2 
ProductionData/M14_r2_spin0000.npz: processData.py ProductionData/M14_r2_spin0000.dat 
	python $< ProductionData/M14_r2_spin0000.dat ProductionData/M14_r2_spin0000 
ProductionData/M14_r2_spin0000.radial: ProductionDataVariables.py ProductionData/M14_r2_spin0000.npz 
	python $< ProductionData/M14_r2_spin0000 ProductionData/M14_r2_spin0000.radial 
# Intermediate results 
ProductionData/M14_r2_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M14_r2_spin0000.npz 
	python $< ProductionData/M14_r2_spin0000 $@ 
Plots/M14_r2_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M14_r2_spin0000rho.dat 
	python $< ProductionData/M14_r2_spin0000rho.dat $@ normal 
Plots/M14_r2_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M14_r2_spin0000rho.dat 
	python $< ProductionData/M14_r2_spin0000rho.dat $@ log 
ProductionData/M14_r2_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M14_r2_spin0000.npz 
	python $< ProductionData/M14_r2_spin0000 $@ 
ProductionData/M14_r2_spin0000grid.dat: ../src/rungyoto ProductionData/M14_r2_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M14_r2_spin0000Gyoto.dat 

ProductionData/M14_r2_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.400000 500 0.998 50 $@ 0 0 0 0.0006671 2 
ProductionData/M14_r2_spin0998.npz: processData.py ProductionData/M14_r2_spin0998.dat 
	python $< ProductionData/M14_r2_spin0998.dat ProductionData/M14_r2_spin0998 
ProductionData/M14_r2_spin0998.radial: ProductionDataVariables.py ProductionData/M14_r2_spin0998.npz 
	python $< ProductionData/M14_r2_spin0998 ProductionData/M14_r2_spin0998.radial 
# Intermediate results 
ProductionData/M14_r2_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M14_r2_spin0998.npz 
	python $< ProductionData/M14_r2_spin0998 $@ 
Plots/M14_r2_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M14_r2_spin0998rho.dat 
	python $< ProductionData/M14_r2_spin0998rho.dat $@ normal 
Plots/M14_r2_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M14_r2_spin0998rho.dat 
	python $< ProductionData/M14_r2_spin0998rho.dat $@ log 
ProductionData/M14_r2_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M14_r2_spin0998.npz 
	python $< ProductionData/M14_r2_spin0998 $@ 
ProductionData/M14_r2_spin0998grid.dat: ../src/rungyoto ProductionData/M14_r2_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M14_r2_spin0998Gyoto.dat 

ProductionData/M15_r1_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.500000 500 0 50 $@ 0 0 0 0.0006671 1 
ProductionData/M15_r1_spin0000.npz: processData.py ProductionData/M15_r1_spin0000.dat 
	python $< ProductionData/M15_r1_spin0000.dat ProductionData/M15_r1_spin0000 
ProductionData/M15_r1_spin0000.radial: ProductionDataVariables.py ProductionData/M15_r1_spin0000.npz 
	python $< ProductionData/M15_r1_spin0000 ProductionData/M15_r1_spin0000.radial 
# Intermediate results 
ProductionData/M15_r1_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M15_r1_spin0000.npz 
	python $< ProductionData/M15_r1_spin0000 $@ 
Plots/M15_r1_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M15_r1_spin0000rho.dat 
	python $< ProductionData/M15_r1_spin0000rho.dat $@ normal 
Plots/M15_r1_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M15_r1_spin0000rho.dat 
	python $< ProductionData/M15_r1_spin0000rho.dat $@ log 
ProductionData/M15_r1_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M15_r1_spin0000.npz 
	python $< ProductionData/M15_r1_spin0000 $@ 
ProductionData/M15_r1_spin0000grid.dat: ../src/rungyoto ProductionData/M15_r1_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M15_r1_spin0000Gyoto.dat 

ProductionData/M15_r1_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.500000 500 0.998 50 $@ 0 0 0 0.0006671 1 
ProductionData/M15_r1_spin0998.npz: processData.py ProductionData/M15_r1_spin0998.dat 
	python $< ProductionData/M15_r1_spin0998.dat ProductionData/M15_r1_spin0998 
ProductionData/M15_r1_spin0998.radial: ProductionDataVariables.py ProductionData/M15_r1_spin0998.npz 
	python $< ProductionData/M15_r1_spin0998 ProductionData/M15_r1_spin0998.radial 
# Intermediate results 
ProductionData/M15_r1_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M15_r1_spin0998.npz 
	python $< ProductionData/M15_r1_spin0998 $@ 
Plots/M15_r1_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M15_r1_spin0998rho.dat 
	python $< ProductionData/M15_r1_spin0998rho.dat $@ normal 
Plots/M15_r1_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M15_r1_spin0998rho.dat 
	python $< ProductionData/M15_r1_spin0998rho.dat $@ log 
ProductionData/M15_r1_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M15_r1_spin0998.npz 
	python $< ProductionData/M15_r1_spin0998 $@ 
ProductionData/M15_r1_spin0998grid.dat: ../src/rungyoto ProductionData/M15_r1_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M15_r1_spin0998Gyoto.dat 

ProductionData/M15_r2_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.500000 500 0 50 $@ 0 0 0 0.0006671 2 
ProductionData/M15_r2_spin0000.npz: processData.py ProductionData/M15_r2_spin0000.dat 
	python $< ProductionData/M15_r2_spin0000.dat ProductionData/M15_r2_spin0000 
ProductionData/M15_r2_spin0000.radial: ProductionDataVariables.py ProductionData/M15_r2_spin0000.npz 
	python $< ProductionData/M15_r2_spin0000 ProductionData/M15_r2_spin0000.radial 
# Intermediate results 
ProductionData/M15_r2_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M15_r2_spin0000.npz 
	python $< ProductionData/M15_r2_spin0000 $@ 
Plots/M15_r2_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M15_r2_spin0000rho.dat 
	python $< ProductionData/M15_r2_spin0000rho.dat $@ normal 
Plots/M15_r2_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M15_r2_spin0000rho.dat 
	python $< ProductionData/M15_r2_spin0000rho.dat $@ log 
ProductionData/M15_r2_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M15_r2_spin0000.npz 
	python $< ProductionData/M15_r2_spin0000 $@ 
ProductionData/M15_r2_spin0000grid.dat: ../src/rungyoto ProductionData/M15_r2_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M15_r2_spin0000Gyoto.dat 

ProductionData/M15_r2_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.500000 500 0.998 50 $@ 0 0 0 0.0006671 2 
ProductionData/M15_r2_spin0998.npz: processData.py ProductionData/M15_r2_spin0998.dat 
	python $< ProductionData/M15_r2_spin0998.dat ProductionData/M15_r2_spin0998 
ProductionData/M15_r2_spin0998.radial: ProductionDataVariables.py ProductionData/M15_r2_spin0998.npz 
	python $< ProductionData/M15_r2_spin0998 ProductionData/M15_r2_spin0998.radial 
# Intermediate results 
ProductionData/M15_r2_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M15_r2_spin0998.npz 
	python $< ProductionData/M15_r2_spin0998 $@ 
Plots/M15_r2_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M15_r2_spin0998rho.dat 
	python $< ProductionData/M15_r2_spin0998rho.dat $@ normal 
Plots/M15_r2_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M15_r2_spin0998rho.dat 
	python $< ProductionData/M15_r2_spin0998rho.dat $@ log 
ProductionData/M15_r2_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M15_r2_spin0998.npz 
	python $< ProductionData/M15_r2_spin0998 $@ 
ProductionData/M15_r2_spin0998grid.dat: ../src/rungyoto ProductionData/M15_r2_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M15_r2_spin0998Gyoto.dat 

ProductionData/M16_r1_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.600000 500 0 50 $@ 0 0 0 0.0006671 1 
ProductionData/M16_r1_spin0000.npz: processData.py ProductionData/M16_r1_spin0000.dat 
	python $< ProductionData/M16_r1_spin0000.dat ProductionData/M16_r1_spin0000 
ProductionData/M16_r1_spin0000.radial: ProductionDataVariables.py ProductionData/M16_r1_spin0000.npz 
	python $< ProductionData/M16_r1_spin0000 ProductionData/M16_r1_spin0000.radial 
# Intermediate results 
ProductionData/M16_r1_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M16_r1_spin0000.npz 
	python $< ProductionData/M16_r1_spin0000 $@ 
Plots/M16_r1_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M16_r1_spin0000rho.dat 
	python $< ProductionData/M16_r1_spin0000rho.dat $@ normal 
Plots/M16_r1_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M16_r1_spin0000rho.dat 
	python $< ProductionData/M16_r1_spin0000rho.dat $@ log 
ProductionData/M16_r1_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M16_r1_spin0000.npz 
	python $< ProductionData/M16_r1_spin0000 $@ 
ProductionData/M16_r1_spin0000grid.dat: ../src/rungyoto ProductionData/M16_r1_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M16_r1_spin0000Gyoto.dat 

ProductionData/M16_r1_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.600000 500 0.998 50 $@ 0 0 0 0.0006671 1 
ProductionData/M16_r1_spin0998.npz: processData.py ProductionData/M16_r1_spin0998.dat 
	python $< ProductionData/M16_r1_spin0998.dat ProductionData/M16_r1_spin0998 
ProductionData/M16_r1_spin0998.radial: ProductionDataVariables.py ProductionData/M16_r1_spin0998.npz 
	python $< ProductionData/M16_r1_spin0998 ProductionData/M16_r1_spin0998.radial 
# Intermediate results 
ProductionData/M16_r1_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M16_r1_spin0998.npz 
	python $< ProductionData/M16_r1_spin0998 $@ 
Plots/M16_r1_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M16_r1_spin0998rho.dat 
	python $< ProductionData/M16_r1_spin0998rho.dat $@ normal 
Plots/M16_r1_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M16_r1_spin0998rho.dat 
	python $< ProductionData/M16_r1_spin0998rho.dat $@ log 
ProductionData/M16_r1_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M16_r1_spin0998.npz 
	python $< ProductionData/M16_r1_spin0998 $@ 
ProductionData/M16_r1_spin0998grid.dat: ../src/rungyoto ProductionData/M16_r1_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M16_r1_spin0998Gyoto.dat 

ProductionData/M16_r2_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.600000 500 0 50 $@ 0 0 0 0.0006671 2 
ProductionData/M16_r2_spin0000.npz: processData.py ProductionData/M16_r2_spin0000.dat 
	python $< ProductionData/M16_r2_spin0000.dat ProductionData/M16_r2_spin0000 
ProductionData/M16_r2_spin0000.radial: ProductionDataVariables.py ProductionData/M16_r2_spin0000.npz 
	python $< ProductionData/M16_r2_spin0000 ProductionData/M16_r2_spin0000.radial 
# Intermediate results 
ProductionData/M16_r2_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M16_r2_spin0000.npz 
	python $< ProductionData/M16_r2_spin0000 $@ 
Plots/M16_r2_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M16_r2_spin0000rho.dat 
	python $< ProductionData/M16_r2_spin0000rho.dat $@ normal 
Plots/M16_r2_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M16_r2_spin0000rho.dat 
	python $< ProductionData/M16_r2_spin0000rho.dat $@ log 
ProductionData/M16_r2_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M16_r2_spin0000.npz 
	python $< ProductionData/M16_r2_spin0000 $@ 
ProductionData/M16_r2_spin0000grid.dat: ../src/rungyoto ProductionData/M16_r2_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M16_r2_spin0000Gyoto.dat 

ProductionData/M16_r2_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.600000 500 0.998 50 $@ 0 0 0 0.0006671 2 
ProductionData/M16_r2_spin0998.npz: processData.py ProductionData/M16_r2_spin0998.dat 
	python $< ProductionData/M16_r2_spin0998.dat ProductionData/M16_r2_spin0998 
ProductionData/M16_r2_spin0998.radial: ProductionDataVariables.py ProductionData/M16_r2_spin0998.npz 
	python $< ProductionData/M16_r2_spin0998 ProductionData/M16_r2_spin0998.radial 
# Intermediate results 
ProductionData/M16_r2_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M16_r2_spin0998.npz 
	python $< ProductionData/M16_r2_spin0998 $@ 
Plots/M16_r2_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M16_r2_spin0998rho.dat 
	python $< ProductionData/M16_r2_spin0998rho.dat $@ normal 
Plots/M16_r2_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M16_r2_spin0998rho.dat 
	python $< ProductionData/M16_r2_spin0998rho.dat $@ log 
ProductionData/M16_r2_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M16_r2_spin0998.npz 
	python $< ProductionData/M16_r2_spin0998 $@ 
ProductionData/M16_r2_spin0998grid.dat: ../src/rungyoto ProductionData/M16_r2_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M16_r2_spin0998Gyoto.dat 

ProductionData/M17_r1_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.700000 500 0 50 $@ 0 0 0 0.0006671 1 
ProductionData/M17_r1_spin0000.npz: processData.py ProductionData/M17_r1_spin0000.dat 
	python $< ProductionData/M17_r1_spin0000.dat ProductionData/M17_r1_spin0000 
ProductionData/M17_r1_spin0000.radial: ProductionDataVariables.py ProductionData/M17_r1_spin0000.npz 
	python $< ProductionData/M17_r1_spin0000 ProductionData/M17_r1_spin0000.radial 
# Intermediate results 
ProductionData/M17_r1_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M17_r1_spin0000.npz 
	python $< ProductionData/M17_r1_spin0000 $@ 
Plots/M17_r1_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M17_r1_spin0000rho.dat 
	python $< ProductionData/M17_r1_spin0000rho.dat $@ normal 
Plots/M17_r1_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M17_r1_spin0000rho.dat 
	python $< ProductionData/M17_r1_spin0000rho.dat $@ log 
ProductionData/M17_r1_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M17_r1_spin0000.npz 
	python $< ProductionData/M17_r1_spin0000 $@ 
ProductionData/M17_r1_spin0000grid.dat: ../src/rungyoto ProductionData/M17_r1_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M17_r1_spin0000Gyoto.dat 

ProductionData/M17_r1_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.700000 500 0.998 50 $@ 0 0 0 0.0006671 1 
ProductionData/M17_r1_spin0998.npz: processData.py ProductionData/M17_r1_spin0998.dat 
	python $< ProductionData/M17_r1_spin0998.dat ProductionData/M17_r1_spin0998 
ProductionData/M17_r1_spin0998.radial: ProductionDataVariables.py ProductionData/M17_r1_spin0998.npz 
	python $< ProductionData/M17_r1_spin0998 ProductionData/M17_r1_spin0998.radial 
# Intermediate results 
ProductionData/M17_r1_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M17_r1_spin0998.npz 
	python $< ProductionData/M17_r1_spin0998 $@ 
Plots/M17_r1_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M17_r1_spin0998rho.dat 
	python $< ProductionData/M17_r1_spin0998rho.dat $@ normal 
Plots/M17_r1_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M17_r1_spin0998rho.dat 
	python $< ProductionData/M17_r1_spin0998rho.dat $@ log 
ProductionData/M17_r1_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M17_r1_spin0998.npz 
	python $< ProductionData/M17_r1_spin0998 $@ 
ProductionData/M17_r1_spin0998grid.dat: ../src/rungyoto ProductionData/M17_r1_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M17_r1_spin0998Gyoto.dat 

ProductionData/M17_r2_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.700000 500 0 50 $@ 0 0 0 0.0006671 2 
ProductionData/M17_r2_spin0000.npz: processData.py ProductionData/M17_r2_spin0000.dat 
	python $< ProductionData/M17_r2_spin0000.dat ProductionData/M17_r2_spin0000 
ProductionData/M17_r2_spin0000.radial: ProductionDataVariables.py ProductionData/M17_r2_spin0000.npz 
	python $< ProductionData/M17_r2_spin0000 ProductionData/M17_r2_spin0000.radial 
# Intermediate results 
ProductionData/M17_r2_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M17_r2_spin0000.npz 
	python $< ProductionData/M17_r2_spin0000 $@ 
Plots/M17_r2_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M17_r2_spin0000rho.dat 
	python $< ProductionData/M17_r2_spin0000rho.dat $@ normal 
Plots/M17_r2_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M17_r2_spin0000rho.dat 
	python $< ProductionData/M17_r2_spin0000rho.dat $@ log 
ProductionData/M17_r2_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M17_r2_spin0000.npz 
	python $< ProductionData/M17_r2_spin0000 $@ 
ProductionData/M17_r2_spin0000grid.dat: ../src/rungyoto ProductionData/M17_r2_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M17_r2_spin0000Gyoto.dat 

ProductionData/M17_r2_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.700000 500 0.998 50 $@ 0 0 0 0.0006671 2 
ProductionData/M17_r2_spin0998.npz: processData.py ProductionData/M17_r2_spin0998.dat 
	python $< ProductionData/M17_r2_spin0998.dat ProductionData/M17_r2_spin0998 
ProductionData/M17_r2_spin0998.radial: ProductionDataVariables.py ProductionData/M17_r2_spin0998.npz 
	python $< ProductionData/M17_r2_spin0998 ProductionData/M17_r2_spin0998.radial 
# Intermediate results 
ProductionData/M17_r2_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M17_r2_spin0998.npz 
	python $< ProductionData/M17_r2_spin0998 $@ 
Plots/M17_r2_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M17_r2_spin0998rho.dat 
	python $< ProductionData/M17_r2_spin0998rho.dat $@ normal 
Plots/M17_r2_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M17_r2_spin0998rho.dat 
	python $< ProductionData/M17_r2_spin0998rho.dat $@ log 
ProductionData/M17_r2_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M17_r2_spin0998.npz 
	python $< ProductionData/M17_r2_spin0998 $@ 
ProductionData/M17_r2_spin0998grid.dat: ../src/rungyoto ProductionData/M17_r2_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M17_r2_spin0998Gyoto.dat 

ProductionData/M18_r1_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.800000 500 0 50 $@ 0 0 0 0.0006671 1 
ProductionData/M18_r1_spin0000.npz: processData.py ProductionData/M18_r1_spin0000.dat 
	python $< ProductionData/M18_r1_spin0000.dat ProductionData/M18_r1_spin0000 
ProductionData/M18_r1_spin0000.radial: ProductionDataVariables.py ProductionData/M18_r1_spin0000.npz 
	python $< ProductionData/M18_r1_spin0000 ProductionData/M18_r1_spin0000.radial 
# Intermediate results 
ProductionData/M18_r1_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M18_r1_spin0000.npz 
	python $< ProductionData/M18_r1_spin0000 $@ 
Plots/M18_r1_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M18_r1_spin0000rho.dat 
	python $< ProductionData/M18_r1_spin0000rho.dat $@ normal 
Plots/M18_r1_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M18_r1_spin0000rho.dat 
	python $< ProductionData/M18_r1_spin0000rho.dat $@ log 
ProductionData/M18_r1_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M18_r1_spin0000.npz 
	python $< ProductionData/M18_r1_spin0000 $@ 
ProductionData/M18_r1_spin0000grid.dat: ../src/rungyoto ProductionData/M18_r1_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M18_r1_spin0000Gyoto.dat 

ProductionData/M18_r1_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.800000 500 0.998 50 $@ 0 0 0 0.0006671 1 
ProductionData/M18_r1_spin0998.npz: processData.py ProductionData/M18_r1_spin0998.dat 
	python $< ProductionData/M18_r1_spin0998.dat ProductionData/M18_r1_spin0998 
ProductionData/M18_r1_spin0998.radial: ProductionDataVariables.py ProductionData/M18_r1_spin0998.npz 
	python $< ProductionData/M18_r1_spin0998 ProductionData/M18_r1_spin0998.radial 
# Intermediate results 
ProductionData/M18_r1_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M18_r1_spin0998.npz 
	python $< ProductionData/M18_r1_spin0998 $@ 
Plots/M18_r1_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M18_r1_spin0998rho.dat 
	python $< ProductionData/M18_r1_spin0998rho.dat $@ normal 
Plots/M18_r1_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M18_r1_spin0998rho.dat 
	python $< ProductionData/M18_r1_spin0998rho.dat $@ log 
ProductionData/M18_r1_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M18_r1_spin0998.npz 
	python $< ProductionData/M18_r1_spin0998 $@ 
ProductionData/M18_r1_spin0998grid.dat: ../src/rungyoto ProductionData/M18_r1_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M18_r1_spin0998Gyoto.dat 

ProductionData/M18_r2_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.800000 500 0 50 $@ 0 0 0 0.0006671 2 
ProductionData/M18_r2_spin0000.npz: processData.py ProductionData/M18_r2_spin0000.dat 
	python $< ProductionData/M18_r2_spin0000.dat ProductionData/M18_r2_spin0000 
ProductionData/M18_r2_spin0000.radial: ProductionDataVariables.py ProductionData/M18_r2_spin0000.npz 
	python $< ProductionData/M18_r2_spin0000 ProductionData/M18_r2_spin0000.radial 
# Intermediate results 
ProductionData/M18_r2_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M18_r2_spin0000.npz 
	python $< ProductionData/M18_r2_spin0000 $@ 
Plots/M18_r2_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M18_r2_spin0000rho.dat 
	python $< ProductionData/M18_r2_spin0000rho.dat $@ normal 
Plots/M18_r2_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M18_r2_spin0000rho.dat 
	python $< ProductionData/M18_r2_spin0000rho.dat $@ log 
ProductionData/M18_r2_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M18_r2_spin0000.npz 
	python $< ProductionData/M18_r2_spin0000 $@ 
ProductionData/M18_r2_spin0000grid.dat: ../src/rungyoto ProductionData/M18_r2_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M18_r2_spin0000Gyoto.dat 

ProductionData/M18_r2_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.800000 500 0.998 50 $@ 0 0 0 0.0006671 2 
ProductionData/M18_r2_spin0998.npz: processData.py ProductionData/M18_r2_spin0998.dat 
	python $< ProductionData/M18_r2_spin0998.dat ProductionData/M18_r2_spin0998 
ProductionData/M18_r2_spin0998.radial: ProductionDataVariables.py ProductionData/M18_r2_spin0998.npz 
	python $< ProductionData/M18_r2_spin0998 ProductionData/M18_r2_spin0998.radial 
# Intermediate results 
ProductionData/M18_r2_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M18_r2_spin0998.npz 
	python $< ProductionData/M18_r2_spin0998 $@ 
Plots/M18_r2_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M18_r2_spin0998rho.dat 
	python $< ProductionData/M18_r2_spin0998rho.dat $@ normal 
Plots/M18_r2_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M18_r2_spin0998rho.dat 
	python $< ProductionData/M18_r2_spin0998rho.dat $@ log 
ProductionData/M18_r2_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M18_r2_spin0998.npz 
	python $< ProductionData/M18_r2_spin0998 $@ 
ProductionData/M18_r2_spin0998grid.dat: ../src/rungyoto ProductionData/M18_r2_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M18_r2_spin0998Gyoto.dat 

ProductionData/M19_r1_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.900000 500 0 50 $@ 0 0 0 0.0006671 1 
ProductionData/M19_r1_spin0000.npz: processData.py ProductionData/M19_r1_spin0000.dat 
	python $< ProductionData/M19_r1_spin0000.dat ProductionData/M19_r1_spin0000 
ProductionData/M19_r1_spin0000.radial: ProductionDataVariables.py ProductionData/M19_r1_spin0000.npz 
	python $< ProductionData/M19_r1_spin0000 ProductionData/M19_r1_spin0000.radial 
# Intermediate results 
ProductionData/M19_r1_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M19_r1_spin0000.npz 
	python $< ProductionData/M19_r1_spin0000 $@ 
Plots/M19_r1_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M19_r1_spin0000rho.dat 
	python $< ProductionData/M19_r1_spin0000rho.dat $@ normal 
Plots/M19_r1_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M19_r1_spin0000rho.dat 
	python $< ProductionData/M19_r1_spin0000rho.dat $@ log 
ProductionData/M19_r1_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M19_r1_spin0000.npz 
	python $< ProductionData/M19_r1_spin0000 $@ 
ProductionData/M19_r1_spin0000grid.dat: ../src/rungyoto ProductionData/M19_r1_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M19_r1_spin0000Gyoto.dat 

ProductionData/M19_r1_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.900000 500 0.998 50 $@ 0 0 0 0.0006671 1 
ProductionData/M19_r1_spin0998.npz: processData.py ProductionData/M19_r1_spin0998.dat 
	python $< ProductionData/M19_r1_spin0998.dat ProductionData/M19_r1_spin0998 
ProductionData/M19_r1_spin0998.radial: ProductionDataVariables.py ProductionData/M19_r1_spin0998.npz 
	python $< ProductionData/M19_r1_spin0998 ProductionData/M19_r1_spin0998.radial 
# Intermediate results 
ProductionData/M19_r1_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M19_r1_spin0998.npz 
	python $< ProductionData/M19_r1_spin0998 $@ 
Plots/M19_r1_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M19_r1_spin0998rho.dat 
	python $< ProductionData/M19_r1_spin0998rho.dat $@ normal 
Plots/M19_r1_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M19_r1_spin0998rho.dat 
	python $< ProductionData/M19_r1_spin0998rho.dat $@ log 
ProductionData/M19_r1_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M19_r1_spin0998.npz 
	python $< ProductionData/M19_r1_spin0998 $@ 
ProductionData/M19_r1_spin0998grid.dat: ../src/rungyoto ProductionData/M19_r1_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M19_r1_spin0998Gyoto.dat 

ProductionData/M19_r2_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.900000 500 0 50 $@ 0 0 0 0.0006671 2 
ProductionData/M19_r2_spin0000.npz: processData.py ProductionData/M19_r2_spin0000.dat 
	python $< ProductionData/M19_r2_spin0000.dat ProductionData/M19_r2_spin0000 
ProductionData/M19_r2_spin0000.radial: ProductionDataVariables.py ProductionData/M19_r2_spin0000.npz 
	python $< ProductionData/M19_r2_spin0000 ProductionData/M19_r2_spin0000.radial 
# Intermediate results 
ProductionData/M19_r2_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M19_r2_spin0000.npz 
	python $< ProductionData/M19_r2_spin0000 $@ 
Plots/M19_r2_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M19_r2_spin0000rho.dat 
	python $< ProductionData/M19_r2_spin0000rho.dat $@ normal 
Plots/M19_r2_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M19_r2_spin0000rho.dat 
	python $< ProductionData/M19_r2_spin0000rho.dat $@ log 
ProductionData/M19_r2_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M19_r2_spin0000.npz 
	python $< ProductionData/M19_r2_spin0000 $@ 
ProductionData/M19_r2_spin0000grid.dat: ../src/rungyoto ProductionData/M19_r2_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M19_r2_spin0000Gyoto.dat 

ProductionData/M19_r2_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.900000 500 0.998 50 $@ 0 0 0 0.0006671 2 
ProductionData/M19_r2_spin0998.npz: processData.py ProductionData/M19_r2_spin0998.dat 
	python $< ProductionData/M19_r2_spin0998.dat ProductionData/M19_r2_spin0998 
ProductionData/M19_r2_spin0998.radial: ProductionDataVariables.py ProductionData/M19_r2_spin0998.npz 
	python $< ProductionData/M19_r2_spin0998 ProductionData/M19_r2_spin0998.radial 
# Intermediate results 
ProductionData/M19_r2_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M19_r2_spin0998.npz 
	python $< ProductionData/M19_r2_spin0998 $@ 
Plots/M19_r2_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M19_r2_spin0998rho.dat 
	python $< ProductionData/M19_r2_spin0998rho.dat $@ normal 
Plots/M19_r2_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M19_r2_spin0998rho.dat 
	python $< ProductionData/M19_r2_spin0998rho.dat $@ log 
ProductionData/M19_r2_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M19_r2_spin0998.npz 
	python $< ProductionData/M19_r2_spin0998 $@ 
ProductionData/M19_r2_spin0998grid.dat: ../src/rungyoto ProductionData/M19_r2_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M19_r2_spin0998Gyoto.dat 

ProductionData/M20_r1_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 2.000000 500 0 50 $@ 0 0 0 0.0006671 1 
ProductionData/M20_r1_spin0000.npz: processData.py ProductionData/M20_r1_spin0000.dat 
	python $< ProductionData/M20_r1_spin0000.dat ProductionData/M20_r1_spin0000 
ProductionData/M20_r1_spin0000.radial: ProductionDataVariables.py ProductionData/M20_r1_spin0000.npz 
	python $< ProductionData/M20_r1_spin0000 ProductionData/M20_r1_spin0000.radial 
# Intermediate results 
ProductionData/M20_r1_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M20_r1_spin0000.npz 
	python $< ProductionData/M20_r1_spin0000 $@ 
Plots/M20_r1_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M20_r1_spin0000rho.dat 
	python $< ProductionData/M20_r1_spin0000rho.dat $@ normal 
Plots/M20_r1_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M20_r1_spin0000rho.dat 
	python $< ProductionData/M20_r1_spin0000rho.dat $@ log 
ProductionData/M20_r1_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M20_r1_spin0000.npz 
	python $< ProductionData/M20_r1_spin0000 $@ 
ProductionData/M20_r1_spin0000grid.dat: ../src/rungyoto ProductionData/M20_r1_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M20_r1_spin0000Gyoto.dat 

ProductionData/M20_r1_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 2.000000 500 0.998 50 $@ 0 0 0 0.0006671 1 
ProductionData/M20_r1_spin0998.npz: processData.py ProductionData/M20_r1_spin0998.dat 
	python $< ProductionData/M20_r1_spin0998.dat ProductionData/M20_r1_spin0998 
ProductionData/M20_r1_spin0998.radial: ProductionDataVariables.py ProductionData/M20_r1_spin0998.npz 
	python $< ProductionData/M20_r1_spin0998 ProductionData/M20_r1_spin0998.radial 
# Intermediate results 
ProductionData/M20_r1_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M20_r1_spin0998.npz 
	python $< ProductionData/M20_r1_spin0998 $@ 
Plots/M20_r1_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M20_r1_spin0998rho.dat 
	python $< ProductionData/M20_r1_spin0998rho.dat $@ normal 
Plots/M20_r1_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M20_r1_spin0998rho.dat 
	python $< ProductionData/M20_r1_spin0998rho.dat $@ log 
ProductionData/M20_r1_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M20_r1_spin0998.npz 
	python $< ProductionData/M20_r1_spin0998 $@ 
ProductionData/M20_r1_spin0998grid.dat: ../src/rungyoto ProductionData/M20_r1_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M20_r1_spin0998Gyoto.dat 

ProductionData/M20_r2_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 2.000000 500 0 50 $@ 0 0 0 0.0006671 2 
ProductionData/M20_r2_spin0000.npz: processData.py ProductionData/M20_r2_spin0000.dat 
	python $< ProductionData/M20_r2_spin0000.dat ProductionData/M20_r2_spin0000 
ProductionData/M20_r2_spin0000.radial: ProductionDataVariables.py ProductionData/M20_r2_spin0000.npz 
	python $< ProductionData/M20_r2_spin0000 ProductionData/M20_r2_spin0000.radial 
# Intermediate results 
ProductionData/M20_r2_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M20_r2_spin0000.npz 
	python $< ProductionData/M20_r2_spin0000 $@ 
Plots/M20_r2_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M20_r2_spin0000rho.dat 
	python $< ProductionData/M20_r2_spin0000rho.dat $@ normal 
Plots/M20_r2_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M20_r2_spin0000rho.dat 
	python $< ProductionData/M20_r2_spin0000rho.dat $@ log 
ProductionData/M20_r2_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M20_r2_spin0000.npz 
	python $< ProductionData/M20_r2_spin0000 $@ 
ProductionData/M20_r2_spin0000grid.dat: ../src/rungyoto ProductionData/M20_r2_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M20_r2_spin0000Gyoto.dat 

ProductionData/M20_r2_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 2.000000 500 0.998 50 $@ 0 0 0 0.0006671 2 
ProductionData/M20_r2_spin0998.npz: processData.py ProductionData/M20_r2_spin0998.dat 
	python $< ProductionData/M20_r2_spin0998.dat ProductionData/M20_r2_spin0998 
ProductionData/M20_r2_spin0998.radial: ProductionDataVariables.py ProductionData/M20_r2_spin0998.npz 
	python $< ProductionData/M20_r2_spin0998 ProductionData/M20_r2_spin0998.radial 
# Intermediate results 
ProductionData/M20_r2_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M20_r2_spin0998.npz 
	python $< ProductionData/M20_r2_spin0998 $@ 
Plots/M20_r2_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M20_r2_spin0998rho.dat 
	python $< ProductionData/M20_r2_spin0998rho.dat $@ normal 
Plots/M20_r2_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M20_r2_spin0998rho.dat 
	python $< ProductionData/M20_r2_spin0998rho.dat $@ log 
ProductionData/M20_r2_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M20_r2_spin0998.npz 
	python $< ProductionData/M20_r2_spin0998 $@ 
ProductionData/M20_r2_spin0998grid.dat: ../src/rungyoto ProductionData/M20_r2_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M20_r2_spin0998Gyoto.dat 

productionr1data: ProductionData/M9_r1_spin0000.dat ProductionData/M9_r1_spin0000.npz ProductionData/M9_r1_spin0000rho.dat ProductionData/M9_r1_spin0000.radial ProductionData/M9_r1_spin0000Gyoto.dat ProductionData/M9_r1_spin0000grid.dat ProductionData/M9_r1_spin0998.dat ProductionData/M9_r1_spin0998.npz ProductionData/M9_r1_spin0998rho.dat ProductionData/M9_r1_spin0998.radial ProductionData/M9_r1_spin0998Gyoto.dat ProductionData/M9_r1_spin0998grid.dat ProductionData/M10_r1_spin0000.dat ProductionData/M10_r1_spin0000.npz ProductionData/M10_r1_spin0000rho.dat ProductionData/M10_r1_spin0000.radial ProductionData/M10_r1_spin0000Gyoto.dat ProductionData/M10_r1_spin0000grid.dat ProductionData/M10_r1_spin0998.dat ProductionData/M10_r1_spin0998.npz ProductionData/M10_r1_spin0998rho.dat ProductionData/M10_r1_spin0998.radial ProductionData/M10_r1_spin0998Gyoto.dat ProductionData/M10_r1_spin0998grid.dat ProductionData/M11_r1_spin0000.dat ProductionData/M11_r1_spin0000.npz ProductionData/M11_r1_spin0000rho.dat ProductionData/M11_r1_spin0000.radial ProductionData/M11_r1_spin0000Gyoto.dat ProductionData/M11_r1_spin0000grid.dat ProductionData/M11_r1_spin0998.dat ProductionData/M11_r1_spin0998.npz ProductionData/M11_r1_spin0998rho.dat ProductionData/M11_r1_spin0998.radial ProductionData/M11_r1_spin0998Gyoto.dat ProductionData/M11_r1_spin0998grid.dat ProductionData/M12_r1_spin0000.dat ProductionData/M12_r1_spin0000.npz ProductionData/M12_r1_spin0000rho.dat ProductionData/M12_r1_spin0000.radial ProductionData/M12_r1_spin0000Gyoto.dat ProductionData/M12_r1_spin0000grid.dat ProductionData/M12_r1_spin0998.dat ProductionData/M12_r1_spin0998.npz ProductionData/M12_r1_spin0998rho.dat ProductionData/M12_r1_spin0998.radial ProductionData/M12_r1_spin0998Gyoto.dat ProductionData/M12_r1_spin0998grid.dat ProductionData/M13_r1_spin0000.dat ProductionData/M13_r1_spin0000.npz ProductionData/M13_r1_spin0000rho.dat ProductionData/M13_r1_spin0000.radial ProductionData/M13_r1_spin0000Gyoto.dat ProductionData/M13_r1_spin0000grid.dat ProductionData/M13_r1_spin0998.dat ProductionData/M13_r1_spin0998.npz ProductionData/M13_r1_spin0998rho.dat ProductionData/M13_r1_spin0998.radial ProductionData/M13_r1_spin0998Gyoto.dat ProductionData/M13_r1_spin0998grid.dat ProductionData/M14_r1_spin0000.dat ProductionData/M14_r1_spin0000.npz ProductionData/M14_r1_spin0000rho.dat ProductionData/M14_r1_spin0000.radial ProductionData/M14_r1_spin0000Gyoto.dat ProductionData/M14_r1_spin0000grid.dat ProductionData/M14_r1_spin0998.dat ProductionData/M14_r1_spin0998.npz ProductionData/M14_r1_spin0998rho.dat ProductionData/M14_r1_spin0998.radial ProductionData/M14_r1_spin0998Gyoto.dat ProductionData/M14_r1_spin0998grid.dat ProductionData/M15_r1_spin0000.dat ProductionData/M15_r1_spin0000.npz ProductionData/M15_r1_spin0000rho.dat ProductionData/M15_r1_spin0000.radial ProductionData/M15_r1_spin0000Gyoto.dat ProductionData/M15_r1_spin0000grid.dat ProductionData/M15_r1_spin0998.dat ProductionData/M15_r1_spin0998.npz ProductionData/M15_r1_spin0998rho.dat ProductionData/M15_r1_spin0998.radial ProductionData/M15_r1_spin0998Gyoto.dat ProductionData/M15_r1_spin0998grid.dat ProductionData/M16_r1_spin0000.dat ProductionData/M16_r1_spin0000.npz ProductionData/M16_r1_spin0000rho.dat ProductionData/M16_r1_spin0000.radial ProductionData/M16_r1_spin0000Gyoto.dat ProductionData/M16_r1_spin0000grid.dat ProductionData/M16_r1_spin0998.dat ProductionData/M16_r1_spin0998.npz ProductionData/M16_r1_spin0998rho.dat ProductionData/M16_r1_spin0998.radial ProductionData/M16_r1_spin0998Gyoto.dat ProductionData/M16_r1_spin0998grid.dat ProductionData/M17_r1_spin0000.dat ProductionData/M17_r1_spin0000.npz ProductionData/M17_r1_spin0000rho.dat ProductionData/M17_r1_spin0000.radial ProductionData/M17_r1_spin0000Gyoto.dat ProductionData/M17_r1_spin0000grid.dat ProductionData/M17_r1_spin0998.dat ProductionData/M17_r1_spin0998.npz ProductionData/M17_r1_spin0998rho.dat ProductionData/M17_r1_spin0998.radial ProductionData/M17_r1_spin0998Gyoto.dat ProductionData/M17_r1_spin0998grid.dat ProductionData/M18_r1_spin0000.dat ProductionData/M18_r1_spin0000.npz ProductionData/M18_r1_spin0000rho.dat ProductionData/M18_r1_spin0000.radial ProductionData/M18_r1_spin0000Gyoto.dat ProductionData/M18_r1_spin0000grid.dat ProductionData/M18_r1_spin0998.dat ProductionData/M18_r1_spin0998.npz ProductionData/M18_r1_spin0998rho.dat ProductionData/M18_r1_spin0998.radial ProductionData/M18_r1_spin0998Gyoto.dat ProductionData/M18_r1_spin0998grid.dat ProductionData/M19_r1_spin0000.dat ProductionData/M19_r1_spin0000.npz ProductionData/M19_r1_spin0000rho.dat ProductionData/M19_r1_spin0000.radial ProductionData/M19_r1_spin0000Gyoto.dat ProductionData/M19_r1_spin0000grid.dat ProductionData/M19_r1_spin0998.dat ProductionData/M19_r1_spin0998.npz ProductionData/M19_r1_spin0998rho.dat ProductionData/M19_r1_spin0998.radial ProductionData/M19_r1_spin0998Gyoto.dat ProductionData/M19_r1_spin0998grid.dat ProductionData/M20_r1_spin0000.dat ProductionData/M20_r1_spin0000.npz ProductionData/M20_r1_spin0000rho.dat ProductionData/M20_r1_spin0000.radial ProductionData/M20_r1_spin0000Gyoto.dat ProductionData/M20_r1_spin0000grid.dat ProductionData/M20_r1_spin0998.dat ProductionData/M20_r1_spin0998.npz ProductionData/M20_r1_spin0998rho.dat ProductionData/M20_r1_spin0998.radial ProductionData/M20_r1_spin0998Gyoto.dat ProductionData/M20_r1_spin0998grid.dat 

productionr1plots: Plots/M9_r1_spin0000rho.pdf Plots/M9_r1_spin0000logrho.pdf Plots/M9_r1_spin0998rho.pdf Plots/M9_r1_spin0998logrho.pdf Plots/M10_r1_spin0000rho.pdf Plots/M10_r1_spin0000logrho.pdf Plots/M10_r1_spin0998rho.pdf Plots/M10_r1_spin0998logrho.pdf Plots/M11_r1_spin0000rho.pdf Plots/M11_r1_spin0000logrho.pdf Plots/M11_r1_spin0998rho.pdf Plots/M11_r1_spin0998logrho.pdf Plots/M12_r1_spin0000rho.pdf Plots/M12_r1_spin0000logrho.pdf Plots/M12_r1_spin0998rho.pdf Plots/M12_r1_spin0998logrho.pdf Plots/M13_r1_spin0000rho.pdf Plots/M13_r1_spin0000logrho.pdf Plots/M13_r1_spin0998rho.pdf Plots/M13_r1_spin0998logrho.pdf Plots/M14_r1_spin0000rho.pdf Plots/M14_r1_spin0000logrho.pdf Plots/M14_r1_spin0998rho.pdf Plots/M14_r1_spin0998logrho.pdf Plots/M15_r1_spin0000rho.pdf Plots/M15_r1_spin0000logrho.pdf Plots/M15_r1_spin0998rho.pdf Plots/M15_r1_spin0998logrho.pdf Plots/M16_r1_spin0000rho.pdf Plots/M16_r1_spin0000logrho.pdf Plots/M16_r1_spin0998rho.pdf Plots/M16_r1_spin0998logrho.pdf Plots/M17_r1_spin0000rho.pdf Plots/M17_r1_spin0000logrho.pdf Plots/M17_r1_spin0998rho.pdf Plots/M17_r1_spin0998logrho.pdf Plots/M18_r1_spin0000rho.pdf Plots/M18_r1_spin0000logrho.pdf Plots/M18_r1_spin0998rho.pdf Plots/M18_r1_spin0998logrho.pdf Plots/M19_r1_spin0000rho.pdf Plots/M19_r1_spin0000logrho.pdf Plots/M19_r1_spin0998rho.pdf Plots/M19_r1_spin0998logrho.pdf Plots/M20_r1_spin0000rho.pdf Plots/M20_r1_spin0000logrho.pdf Plots/M20_r1_spin0998rho.pdf Plots/M20_r1_spin0998logrho.pdf 

productionr2data: ProductionData/M9_r2_spin0000.dat ProductionData/M9_r2_spin0000.npz ProductionData/M9_r2_spin0000rho.dat ProductionData/M9_r2_spin0000.radial ProductionData/M9_r2_spin0000Gyoto.dat ProductionData/M9_r2_spin0000grid.dat ProductionData/M9_r2_spin0998.dat ProductionData/M9_r2_spin0998.npz ProductionData/M9_r2_spin0998rho.dat ProductionData/M9_r2_spin0998.radial ProductionData/M9_r2_spin0998Gyoto.dat ProductionData/M9_r2_spin0998grid.dat ProductionData/M10_r2_spin0000.dat ProductionData/M10_r2_spin0000.npz ProductionData/M10_r2_spin0000rho.dat ProductionData/M10_r2_spin0000.radial ProductionData/M10_r2_spin0000Gyoto.dat ProductionData/M10_r2_spin0000grid.dat ProductionData/M10_r2_spin0998.dat ProductionData/M10_r2_spin0998.npz ProductionData/M10_r2_spin0998rho.dat ProductionData/M10_r2_spin0998.radial ProductionData/M10_r2_spin0998Gyoto.dat ProductionData/M10_r2_spin0998grid.dat ProductionData/M11_r2_spin0000.dat ProductionData/M11_r2_spin0000.npz ProductionData/M11_r2_spin0000rho.dat ProductionData/M11_r2_spin0000.radial ProductionData/M11_r2_spin0000Gyoto.dat ProductionData/M11_r2_spin0000grid.dat ProductionData/M11_r2_spin0998.dat ProductionData/M11_r2_spin0998.npz ProductionData/M11_r2_spin0998rho.dat ProductionData/M11_r2_spin0998.radial ProductionData/M11_r2_spin0998Gyoto.dat ProductionData/M11_r2_spin0998grid.dat ProductionData/M12_r2_spin0000.dat ProductionData/M12_r2_spin0000.npz ProductionData/M12_r2_spin0000rho.dat ProductionData/M12_r2_spin0000.radial ProductionData/M12_r2_spin0000Gyoto.dat ProductionData/M12_r2_spin0000grid.dat ProductionData/M12_r2_spin0998.dat ProductionData/M12_r2_spin0998.npz ProductionData/M12_r2_spin0998rho.dat ProductionData/M12_r2_spin0998.radial ProductionData/M12_r2_spin0998Gyoto.dat ProductionData/M12_r2_spin0998grid.dat ProductionData/M13_r2_spin0000.dat ProductionData/M13_r2_spin0000.npz ProductionData/M13_r2_spin0000rho.dat ProductionData/M13_r2_spin0000.radial ProductionData/M13_r2_spin0000Gyoto.dat ProductionData/M13_r2_spin0000grid.dat ProductionData/M13_r2_spin0998.dat ProductionData/M13_r2_spin0998.npz ProductionData/M13_r2_spin0998rho.dat ProductionData/M13_r2_spin0998.radial ProductionData/M13_r2_spin0998Gyoto.dat ProductionData/M13_r2_spin0998grid.dat ProductionData/M14_r2_spin0000.dat ProductionData/M14_r2_spin0000.npz ProductionData/M14_r2_spin0000rho.dat ProductionData/M14_r2_spin0000.radial ProductionData/M14_r2_spin0000Gyoto.dat ProductionData/M14_r2_spin0000grid.dat ProductionData/M14_r2_spin0998.dat ProductionData/M14_r2_spin0998.npz ProductionData/M14_r2_spin0998rho.dat ProductionData/M14_r2_spin0998.radial ProductionData/M14_r2_spin0998Gyoto.dat ProductionData/M14_r2_spin0998grid.dat ProductionData/M15_r2_spin0000.dat ProductionData/M15_r2_spin0000.npz ProductionData/M15_r2_spin0000rho.dat ProductionData/M15_r2_spin0000.radial ProductionData/M15_r2_spin0000Gyoto.dat ProductionData/M15_r2_spin0000grid.dat ProductionData/M15_r2_spin0998.dat ProductionData/M15_r2_spin0998.npz ProductionData/M15_r2_spin0998rho.dat ProductionData/M15_r2_spin0998.radial ProductionData/M15_r2_spin0998Gyoto.dat ProductionData/M15_r2_spin0998grid.dat ProductionData/M16_r2_spin0000.dat ProductionData/M16_r2_spin0000.npz ProductionData/M16_r2_spin0000rho.dat ProductionData/M16_r2_spin0000.radial ProductionData/M16_r2_spin0000Gyoto.dat ProductionData/M16_r2_spin0000grid.dat ProductionData/M16_r2_spin0998.dat ProductionData/M16_r2_spin0998.npz ProductionData/M16_r2_spin0998rho.dat ProductionData/M16_r2_spin0998.radial ProductionData/M16_r2_spin0998Gyoto.dat ProductionData/M16_r2_spin0998grid.dat ProductionData/M17_r2_spin0000.dat ProductionData/M17_r2_spin0000.npz ProductionData/M17_r2_spin0000rho.dat ProductionData/M17_r2_spin0000.radial ProductionData/M17_r2_spin0000Gyoto.dat ProductionData/M17_r2_spin0000grid.dat ProductionData/M17_r2_spin0998.dat ProductionData/M17_r2_spin0998.npz ProductionData/M17_r2_spin0998rho.dat ProductionData/M17_r2_spin0998.radial ProductionData/M17_r2_spin0998Gyoto.dat ProductionData/M17_r2_spin0998grid.dat ProductionData/M18_r2_spin0000.dat ProductionData/M18_r2_spin0000.npz ProductionData/M18_r2_spin0000rho.dat ProductionData/M18_r2_spin0000.radial ProductionData/M18_r2_spin0000Gyoto.dat ProductionData/M18_r2_spin0000grid.dat ProductionData/M18_r2_spin0998.dat ProductionData/M18_r2_spin0998.npz ProductionData/M18_r2_spin0998rho.dat ProductionData/M18_r2_spin0998.radial ProductionData/M18_r2_spin0998Gyoto.dat ProductionData/M18_r2_spin0998grid.dat ProductionData/M19_r2_spin0000.dat ProductionData/M19_r2_spin0000.npz ProductionData/M19_r2_spin0000rho.dat ProductionData/M19_r2_spin0000.radial ProductionData/M19_r2_spin0000Gyoto.dat ProductionData/M19_r2_spin0000grid.dat ProductionData/M19_r2_spin0998.dat ProductionData/M19_r2_spin0998.npz ProductionData/M19_r2_spin0998rho.dat ProductionData/M19_r2_spin0998.radial ProductionData/M19_r2_spin0998Gyoto.dat ProductionData/M19_r2_spin0998grid.dat ProductionData/M20_r2_spin0000.dat ProductionData/M20_r2_spin0000.npz ProductionData/M20_r2_spin0000rho.dat ProductionData/M20_r2_spin0000.radial ProductionData/M20_r2_spin0000Gyoto.dat ProductionData/M20_r2_spin0000grid.dat ProductionData/M20_r2_spin0998.dat ProductionData/M20_r2_spin0998.npz ProductionData/M20_r2_spin0998rho.dat ProductionData/M20_r2_spin0998.radial ProductionData/M20_r2_spin0998Gyoto.dat ProductionData/M20_r2_spin0998grid.dat 

productionr2plots: Plots/M9_r2_spin0000rho.pdf Plots/M9_r2_spin0000logrho.pdf Plots/M9_r2_spin0998rho.pdf Plots/M9_r2_spin0998logrho.pdf Plots/M10_r2_spin0000rho.pdf Plots/M10_r2_spin0000logrho.pdf Plots/M10_r2_spin0998rho.pdf Plots/M10_r2_spin0998logrho.pdf Plots/M11_r2_spin0000rho.pdf Plots/M11_r2_spin0000logrho.pdf Plots/M11_r2_spin0998rho.pdf Plots/M11_r2_spin0998logrho.pdf Plots/M12_r2_spin0000rho.pdf Plots/M12_r2_spin0000logrho.pdf Plots/M12_r2_spin0998rho.pdf Plots/M12_r2_spin0998logrho.pdf Plots/M13_r2_spin0000rho.pdf Plots/M13_r2_spin0000logrho.pdf Plots/M13_r2_spin0998rho.pdf Plots/M13_r2_spin0998logrho.pdf Plots/M14_r2_spin0000rho.pdf Plots/M14_r2_spin0000logrho.pdf Plots/M14_r2_spin0998rho.pdf Plots/M14_r2_spin0998logrho.pdf Plots/M15_r2_spin0000rho.pdf Plots/M15_r2_spin0000logrho.pdf Plots/M15_r2_spin0998rho.pdf Plots/M15_r2_spin0998logrho.pdf Plots/M16_r2_spin0000rho.pdf Plots/M16_r2_spin0000logrho.pdf Plots/M16_r2_spin0998rho.pdf Plots/M16_r2_spin0998logrho.pdf Plots/M17_r2_spin0000rho.pdf Plots/M17_r2_spin0000logrho.pdf Plots/M17_r2_spin0998rho.pdf Plots/M17_r2_spin0998logrho.pdf Plots/M18_r2_spin0000rho.pdf Plots/M18_r2_spin0000logrho.pdf Plots/M18_r2_spin0998rho.pdf Plots/M18_r2_spin0998logrho.pdf Plots/M19_r2_spin0000rho.pdf Plots/M19_r2_spin0000logrho.pdf Plots/M19_r2_spin0998rho.pdf Plots/M19_r2_spin0998logrho.pdf Plots/M20_r2_spin0000rho.pdf Plots/M20_r2_spin0000logrho.pdf Plots/M20_r2_spin0998rho.pdf Plots/M20_r2_spin0998logrho.pdf 

ProductionData/M10_r1_spin0000.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.000000 500 0 50 $@ 0 0 0 0.0006671 1 
ProductionData/M10_r1_spin0000.npz: processData.py ProductionData/M10_r1_spin0000.dat 
	python $< ProductionData/M10_r1_spin0000.dat ProductionData/M10_r1_spin0000 
ProductionData/M10_r1_spin0000.radial: ProductionDataVariables.py ProductionData/M10_r1_spin0000.npz 
	python $< ProductionData/M10_r1_spin0000 ProductionData/M10_r1_spin0000.radial 
# Intermediate results 
ProductionData/M10_r1_spin0000rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M10_r1_spin0000.npz 
	python $< ProductionData/M10_r1_spin0000 $@ 
Plots/M10_r1_spin0000rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M10_r1_spin0000rho.dat 
	python $< ProductionData/M10_r1_spin0000rho.dat $@ normal 
Plots/M10_r1_spin0000logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M10_r1_spin0000rho.dat 
	python $< ProductionData/M10_r1_spin0000rho.dat $@ log 
ProductionData/M10_r1_spin0000Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M10_r1_spin0000.npz 
	python $< ProductionData/M10_r1_spin0000 $@ 
ProductionData/M10_r1_spin0000grid.dat: ../src/rungyoto ProductionData/M10_r1_spin0000Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M10_r1_spin0000Gyoto.dat 

ProductionData/M10_r1_spin0300.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.000000 500 0.3 50 $@ 0 0 0 0.0006671 1 
ProductionData/M10_r1_spin0300.npz: processData.py ProductionData/M10_r1_spin0300.dat 
	python $< ProductionData/M10_r1_spin0300.dat ProductionData/M10_r1_spin0300 
ProductionData/M10_r1_spin0300.radial: ProductionDataVariables.py ProductionData/M10_r1_spin0300.npz 
	python $< ProductionData/M10_r1_spin0300 ProductionData/M10_r1_spin0300.radial 
# Intermediate results 
ProductionData/M10_r1_spin0300rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M10_r1_spin0300.npz 
	python $< ProductionData/M10_r1_spin0300 $@ 
Plots/M10_r1_spin0300rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M10_r1_spin0300rho.dat 
	python $< ProductionData/M10_r1_spin0300rho.dat $@ normal 
Plots/M10_r1_spin0300logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M10_r1_spin0300rho.dat 
	python $< ProductionData/M10_r1_spin0300rho.dat $@ log 
ProductionData/M10_r1_spin0300Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M10_r1_spin0300.npz 
	python $< ProductionData/M10_r1_spin0300 $@ 
ProductionData/M10_r1_spin0300grid.dat: ../src/rungyoto ProductionData/M10_r1_spin0300Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M10_r1_spin0300Gyoto.dat 

ProductionData/M10_r1_spin0600.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.000000 500 0.6 50 $@ 0 0 0 0.0006671 1 
ProductionData/M10_r1_spin0600.npz: processData.py ProductionData/M10_r1_spin0600.dat 
	python $< ProductionData/M10_r1_spin0600.dat ProductionData/M10_r1_spin0600 
ProductionData/M10_r1_spin0600.radial: ProductionDataVariables.py ProductionData/M10_r1_spin0600.npz 
	python $< ProductionData/M10_r1_spin0600 ProductionData/M10_r1_spin0600.radial 
# Intermediate results 
ProductionData/M10_r1_spin0600rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M10_r1_spin0600.npz 
	python $< ProductionData/M10_r1_spin0600 $@ 
Plots/M10_r1_spin0600rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M10_r1_spin0600rho.dat 
	python $< ProductionData/M10_r1_spin0600rho.dat $@ normal 
Plots/M10_r1_spin0600logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M10_r1_spin0600rho.dat 
	python $< ProductionData/M10_r1_spin0600rho.dat $@ log 
ProductionData/M10_r1_spin0600Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M10_r1_spin0600.npz 
	python $< ProductionData/M10_r1_spin0600 $@ 
ProductionData/M10_r1_spin0600grid.dat: ../src/rungyoto ProductionData/M10_r1_spin0600Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M10_r1_spin0600Gyoto.dat 

ProductionData/M10_r1_spin0998.dat: ../src/growthsampling 
	OMP_NUM_THREADS=8 $< 1200000 1.000000 500 0.998 50 $@ 0 0 0 0.0006671 1 
ProductionData/M10_r1_spin0998.npz: processData.py ProductionData/M10_r1_spin0998.dat 
	python $< ProductionData/M10_r1_spin0998.dat ProductionData/M10_r1_spin0998 
ProductionData/M10_r1_spin0998.radial: ProductionDataVariables.py ProductionData/M10_r1_spin0998.npz 
	python $< ProductionData/M10_r1_spin0998 ProductionData/M10_r1_spin0998.radial 
# Intermediate results 
ProductionData/M10_r1_spin0998rho.dat: ../PostProcessing/calculateDensity.py ProductionData/M10_r1_spin0998.npz 
	python $< ProductionData/M10_r1_spin0998 $@ 
Plots/M10_r1_spin0998rho.pdf: ../PostProcessing/plotDensity.py ProductionData/M10_r1_spin0998rho.dat 
	python $< ProductionData/M10_r1_spin0998rho.dat $@ normal 
Plots/M10_r1_spin0998logrho.pdf: ../PostProcessing/plotDensity.py ProductionData/M10_r1_spin0998rho.dat 
	python $< ProductionData/M10_r1_spin0998rho.dat $@ log 
ProductionData/M10_r1_spin0998Gyoto.dat: ../Tests/convertDataToGyoto.py ProductionData/M10_r1_spin0998.npz 
	python $< ProductionData/M10_r1_spin0998 $@ 
ProductionData/M10_r1_spin0998grid.dat: ../src/rungyoto ProductionData/M10_r1_spin0998Gyoto.dat 
	OMP_NUM_THREADS=8 $< runge_kutta_fehlberg78 200 30 1 200 $@ 1 ProductionData/M10_r1_spin0998Gyoto.dat 

production3ddata: ProductionData/M10_r1_spin0000.dat ProductionData/M10_r1_spin0000.npz ProductionData/M10_r1_spin0000rho.dat ProductionData/M10_r1_spin0000.radial ProductionData/M10_r1_spin0000grid.dat ProductionData/M10_r1_spin0300.dat ProductionData/M10_r1_spin0300.npz ProductionData/M10_r1_spin0300rho.dat ProductionData/M10_r1_spin0300.radial ProductionData/M10_r1_spin0300grid.dat ProductionData/M10_r1_spin0600.dat ProductionData/M10_r1_spin0600.npz ProductionData/M10_r1_spin0600rho.dat ProductionData/M10_r1_spin0600.radial ProductionData/M10_r1_spin0600grid.dat ProductionData/M10_r1_spin0998.dat ProductionData/M10_r1_spin0998.npz ProductionData/M10_r1_spin0998rho.dat ProductionData/M10_r1_spin0998.radial ProductionData/M10_r1_spin0998grid.dat 

production3dplots: Plots/M10_r1_spin0000rho.pdf Plots/M10_r1_spin0000logrho.pdf Plots/M10_r1_spin0300rho.pdf Plots/M10_r1_spin0300logrho.pdf Plots/M10_r1_spin0600rho.pdf Plots/M10_r1_spin0600logrho.pdf Plots/M10_r1_spin0998rho.pdf Plots/M10_r1_spin0998logrho.pdf 

###################################################################################################################################
###################################################################################################################################
###################################################################################################################################

development: 

all: data datap plotdata plots

testplot: ../Tests/createComparisonPlot.py ProductionData/M10_r1_spin0998.npz ProductionData/M10_r1_spin0998grid.dat
	python ../Tests/createComparisonPlot.py ProductionData/M10_r1_spin0998 ProductionData/M10_r1_spin0998grid.dat /tmp/test.pdf

productionfigures: ProductionDataVariables.py  ProductionDensityProfileSpin.py  ProductionFigureFit.py  ProductionFigureGrowthOfAnnihilation.py  ProductionFigureRho.py
	python ProductionDataVariables.py
	python ProductionDensityProfileSpin.py
	python ProductionFigureFit.py
	python ProductionFigureGrowthOfAnnihilation.py false normal
	python ProductionFigureRho.py





