#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=12:mem=190gb

module load anaconda3/personal

source activate py27                                                                                            #Load python 2.7 environment


mkdir -p NREPL="$NREPL"/Simulations_"$epoch"epoch
mkdir -p NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method

echo "Start copying simulations: " `date`
cp -r $WORK/Work/Data/test_preprocessing/NREPL="$NREPL"/Simulations_"$epoch"epoch $TMPDIR/NREPL="$NREPL"/.
echo "End copying simulations: " `date`
echo "Start copying images: " `date`
cp -r $WORK/Work/Data/test_preprocessing/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/$ordering_option $TMPDIR/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/.
echo "End copying images: " `date`

                                    #Copy python scripts
cp -r $WORK/Work/Scripts/HPC_test_preprocessing/HPC_run_cnn_CrossEntropy_copy_5.py $TMPDIR/.
#cp -r $WORK/Work/Scripts/HPC_test_preprocessing/HPC_run_cnn_CrossEntropy_copy_5_Leaky.py $TMPDIR/.
#cp -r $WORK/Work/Scripts/HPC_test_preprocessing/HPC_run_cnn_CrossEntropy_copy_41.py $TMPDIR/.
#cp -r $WORK/Work/Scripts/methods_comparison/ANN/HPC_run_ANN_KL.py $TMPDIR/.

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

        IMAGES_FOLDER_IN_WORK_DIR=$WORK/Work/Data/test_preprocessing/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/$ordering_option/.

        echo "Start cnn_cat: "`date`
        python HPC_run_cnn_CrossEntropy_copy_5.py $IN_PATH $OUT_PATH $channel_options $raw_resize_option1_number $raw_resize_option2_number $batch_size $epochs $filters $kernel_size $pooling_size $raw_resize_option 1                #Run CNN (24 cores)
        echo "End cnn_cat: "`date`

        PWD=$TMPDIR/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/$ordering_option/CNN_"$raw_resize_option"_Categorical
        cp -r $PWD $WORK/Work/Data/test_preprocessing/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/$ordering_option/CNN_"$raw_resize_option"_Categorical_$batch_size



        #echo "Start cnn_cat: "`date`
        #python HPC_run_cnn_CrossEntropy_copy_5_Leaky.py $IN_PATH $OUT_PATH $channel_options $raw_resize_option1_number $raw_resize_option2_number $batch_size $epochs $filters $kernel_size $pooling_size $raw_resize_option 1                #Run CNN (24 cores)
        #echo "End cnn_cat: "`date`

        #PWD=$TMPDIR/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/$ordering_option/CNN_"$raw_resize_option"_Categorical
        #cp -r $PWD $WORK/Work/Data/test_preprocessing/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/$ordering_option/CNN_"$raw_resize_option"_Categorical_Leaky



        #echo "Start cnn_cat: "`date`
        #python HPC_run_cnn_CrossEntropy_copy_41.py $IN_PATH $OUT_PATH $channel_options $raw_resize_option1_number $raw_resize_option2_number $batch_size $epochs $filters $kernel_size $pooling_size $raw_resize_option 1                #Run CNN (24 cores)
        #echo "End cnn_cat: "`date`

        #PWD=$TMPDIR/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/$ordering_option/CNN_"$raw_resize_option"_Categorical
        #cp -r $PWD $WORK/Work/Data/test_preprocessing/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/$ordering_option/CNN_"$raw_resize_option"_Categorical_41



        #echo "Start cnn_KL: "`date`
        #python HPC_run_cnn_CrossEntropy_copy_5.py $IN_PATH $OUT_PATH $channel_options $raw_resize_option1_number $raw_resize_option2_number $batch_size $epochs $filters $kernel_size $pooling_size $raw_resize_option 3            #Run CNN (24 cores)
        #echo "End cnn_KL: "`date`

        #PWD=$TMPDIR/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/$ordering_option/CNN_"$raw_resize_option"_PMF
        #cp -r $PWD $IMAGES_FOLDER_IN_WORK_DIR


        #echo "Start ann_cat: "`date`
        #python HPC_run_ANN_CrossEntropy.py $IN_PATH $OUT_PATH $channel_options $raw_resize_option1_number $raw_resize_option2_number $batch_size $epochs $filters $kernel_size $pooling_size $raw_resize_option 1            #Run CNN (24 cores)
        #echo "End ann_cat: "`date`

        #PWD=$TMPDIR/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/$ordering_option/ANN_"$raw_resize_option"_KLD
        #cp -r $PWD $IMAGES_FOLDER_IN_WORK_DIR


        #echo "Start ann_cat: "`date`
        #python HPC_run_ANN_CrossEntropy.py $IN_PATH $OUT_PATH $channel_options $raw_resize_option1_number $raw_resize_option2_number $batch_size $epochs $filters $kernel_size $pooling_size $raw_resize_option 3            #Run CNN (24 cores)
        #echo "End ann_cat: "`date`

        #PWD=$TMPDIR/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/$ordering_option/ANN_"$raw_resize_option"_KLD
        #cp -r $PWD $IMAGES_FOLDER_IN_WORK_DIR



        raw_resize_option1_number=$(($raw_resize_option1_number - 4))

    done
