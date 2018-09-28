#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=12:mem=190gb

module load anaconda3/personal

source activate py27                                                                                            #Load python 2.7 environment


######### Working locally in the HPC server: copy information from home directory to the HPC allocated computers #########


mkdir -p NREPL="$NREPL"/Simulations_"$epoch"epoch
mkdir -p NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method

cp -r $WORK/Work/Data/test_preprocessing/NREPL="$NREPL"/Simulations_"$epoch"epoch $TMPDIR/NREPL="$NREPL"/.
cp -r $WORK/Work/Data/test_preprocessing/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/$ordering_option $TMPDIR/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/.
cp -r $WORK/Work/Scripts/HPC_test_preprocessing/"$dfn"CNN"$nb_classes".py $TMPDIR/.


######### Set


channel_options=1
#raw_resize_option1_number=1        # 1: raw , 5: 128x128                                                                                      #Fixed image_creation options
raw_resize_option2_number=1         # 1: mean

IN_PATH=NREPL="$NREPL"/Simulations_"$epoch"epoch
OUT_PATH=NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/$ordering_option
mkdir -p OUT_PATH


################# 4. CNN ###################

#argv[12]:
#training_options=1 --> Categorical
#training_options=3 --> PMF

raw_resize_option1_number=5
for raw_resize_option in $( echo 128x128; ) #; echo raw_mean )   # raw_resize_option1_number=1: raw (raw_resize_options_2=1 only applies if raw is chosen, and 1 means 'mean' size), raw_resize_option_number=5: 128x128
do
        python "$dfn"CNN"$nb_classes".py $IN_PATH $OUT_PATH $channel_options $raw_resize_option1_number $raw_resize_option2_number $batch_size $epochs $filters $kernel_size $pooling_size $raw_resize_option 3


        ###### Copy back the results ######
        PWD=$TMPDIR/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/$ordering_option/CNN_"$raw_resize_option"_PMF$nb_classes
        OUT=$WORK/Work/Data/test_preprocessing/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/$ordering_option/"$dfn"CNN_SORTING_"$nb_classes"_$ordering_option

        cp -r $PWD $OUT


        raw_resize_option1_number=$(($raw_resize_option1_number - 4)) # Jump to the next option.
done
