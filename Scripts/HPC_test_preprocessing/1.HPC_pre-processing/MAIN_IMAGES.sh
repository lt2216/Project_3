#!/bin/bash

# This fragments the *task* into 8 HPC large-memory jobs. *task*: create indexes and images from simulations in a combinatorial way for all the pre-processing options under assessment.

# Requirements: previous running of qsub -v NREPL=1000 HPC_main_simulations.sh


# sh MAIN_IMAGES.sh 10000

for epoch in 1
    do
        colour_coding_method_number=1
        #; echo A_D )
        for colour_coding_method in $( echo M_m )
            do
                ordering_option_number=1
                for ordering_option in $( echo order_rows_and_cols; ) #no_ordering; echo order_rows; echo order_cols; echo order_rows_and_cols; )
                    do
                        qsub -v NREPL=$1,epoch=$epoch,colour_coding_method=$colour_coding_method,colour_coding_method_number=$colour_coding_method_number,ordering_option=$ordering_option,ordering_option_number=$ordering_option_number HPC_images_arg-epoch.sh
                        ordering_option_number=$(($ordering_option_number + 1))
                    done
                colour_coding_method_number=$(($colour_coding_method_number + 1))
            done
    done
