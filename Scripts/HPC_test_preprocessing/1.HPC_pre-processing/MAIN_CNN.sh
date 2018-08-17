#!/bin/bash

# Requirements: previous running of [qsub -v NREPL=1000 HPC_main_simulations.sh] and [sh MAIN_IMAGES.sh]


# sh MAIN_CNN.sh NREPL batch_size epochs filters kernel_size pooling_size --> sh MAIN_CNN.sh 1000 3 20 32 3 2

NREPL=$1
#2 3
for epoch in 1
    do
        #; echo A_D )
        for colour_coding_method in $( echo M_m )
            do
                for ordering_option in $( echo order_rows; )  #$( echo no_ordering; echo order_rows; echo order_cols; echo order_rows_and_cols; )
                    do
                        qsub -v NREPL=$NREPL,epoch=$epoch,colour_coding_method=$colour_coding_method,ordering_option=$ordering_option,batch_size=$2,epochs=$3,filters=$4,kernel_size=$5,pooling_size=$6 HPC_CNN_arg-epoch.sh
                    done
            done
    done
