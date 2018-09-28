#!/bin/bash


# RUNNING:  (sh MAIN_CNN.sh NREPL batch_size epochs filters kernel_size pooling_size) --> sh MAIN_CNN.sh 10000 32 10 32 3 2


for epoch in 1 2 3
do
    for colour_coding_method in $( echo M_m )
    do
        for ordering_option in $( echo no_ordering; echo order_rows; echo order_cols; echo order_rows_and_cols; )
        do
            for nb_classes in 41 #5 10 20
            do
                dfn=06        # Variance of the Gaussian used to define true labels (only this option so far. The 06CNN41.py script can be modifies (also the name e.g. 03CNN41) if this option is to be changed)
                qsub -v NREPL=$1,epoch=$epoch,colour_coding_method=$colour_coding_method,ordering_option=$ordering_option,batch_size=$2,epochs=$3,filters=$4,kernel_size=$5,pooling_size=$6,nb_classes=$nb_classes,dfn=$dfn HPC_CNN_arg-epoch.sh
            done
        done
    done
done
