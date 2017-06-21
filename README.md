# MARBLE
Project to investigate dark matter spikes near rotating black holes in full general relativity

# Dependencies:
- Gyoto
- boost
- openmp


# INSTALLATION (Paris)
cd Paris
less INSTALL # Read installation instructions

# INSTALLATION:
./bootstrap.sh
# Configure with any flags
CXXFLAGS="-O3 -march=native" ./configure --prefix=${HOME}/MARBLE
cd src && make -j install && cd -

# Run tests:
cd Tests && make && cd -

# Run production run
cd PostProcessing && make productionr2data && cd -


