#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=12:mem=190gb

module load anaconda3/personal

source activate py27                                                                                            #Load python 2.7 environment

################# 1. Run simulations ###################

cp -r $WORK/msms $TMPDIR/.                                                                                      #Copy msms
cp -r $WORK/Work/Scripts/HPC_test_preprocessing/Simulations/HPC_create_index.py $TMPDIR/.                       #Copy script


for epoch in 1 2 3
    do
        cp -r $WORK/Work/Scripts/HPC_test_preprocessing/Simulations/epoch"$epoch".sh $TMPDIR/.                  #Copy sh epoch-specific scripts.
        mkdir -p NREPL="$NREPL"/Simulations_"$epoch"epoch                                                       #Create main output folder (this is what will be copied back in the las line) and three subdirectories for the simulations
        sh epoch"$epoch".sh $NREPL

        ################# 2. Create Indexes ###################

        python HPC_create_index.py NREPL="$NREPL"/Simulations_"$epoch"epoch                                             #Create index in each simulations folder (1 per epoch)
    done

################# 3. Copy results back ###################

cp -r $TMPDIR/NREPL="$NREPL" $WORK/Work/Data/test_preprocessing/.            #Copy back results
