#!/bin/bash
#PBS -l walltime=15:00:00
#PBS -l select=1:ncpus=12:mem=190gb

module load anaconda3/personal

source activate py27                                                                                            #Load python 2.7 environment


################# 1. Copy simulations ###################

mkdir -p $TMPDIR/NREPL="$NREPL"

echo "Start copying simulations: "`date`
cp -r $WORK/Work/Data/test_preprocessing/NREPL="$NREPL"/Simulations_"$epoch"epoch $TMPDIR/NREPL="$NREPL"/.
echo "End copying simulations: "`date`


################# 2. Create Images ###################

cp -r $WORK/Work/Scripts/HPC_test_preprocessing/1.HPC_pre-processing/HPC_create_imgs_copy.py $TMPDIR/.                             #Copy python scripts                                                                                           #Copy msms


IN_PATH=NREPL="$NREPL"/Simulations_"$epoch"epoch
OUT_PATH=NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/$ordering_option

mkdir -p $OUT_PATH                                                                                                      #Create an output folder

data_type=1
channel_options=1
#raw_resize_option1_number=1        # 1: raw , 5: 128x128                                                                                      #Fixed image_creation options
raw_resize_option2_number=1        # 1: mean
app_thes=2
feature=2

echo "Start creating images: "`date`

python HPC_create_imgs_copy.py $data_type $IN_PATH $OUT_PATH $ordering_option_number $app_thes $colour_coding_method_number          #Create images

echo "End creating images: "`date`


################# 3. Copy results back ###################

echo "Start copying images back: "`date`

WORK_DIR=$WORK/Work/Data/test_preprocessing/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method

mkdir -p $WORK_DIR
cp -r $TMPDIR/NREPL="$NREPL"/Images_"$epoch"epoch/$colour_coding_method/$ordering_option $WORK_DIR/.           #Copy back results

echo "End copying images back: "`date`
