#!/bin/bash
########################################################################
####  script for paralization of quantile estimator simulation of Cope
####  processes
########################################################################
# go to folder where the simulation script is stored
cd ..

# Loop over parallelisations
for parallel_count in {1..16}
do
    matlab -nodesktop -nosplash -r "Sim_EstimCopeQuantilesSNR('$parallel_count'); exit" &
done
