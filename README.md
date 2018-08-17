Setting up the environment:

1. Clone repository.
2. Set working directory to 'Scripts/HPC_test_preprocessing/1.HPC_pre-processing'

\\
Pipeline:

1. Simulations for the 1,2 and 3 epochs models

run: 
```
HPC_main_simulations.sh NREPL
```
(NREPL is an integer equal to the number of replicates. One simulations folder will be placed in the Data folder for each epoch).


2. Generate training, validation and testing datasets (proportions: 0.5, 0.25, 0.25), with 4 different pre-processing options per epoch model. You can choose just one of them by easily editing the for loop inside the MAIN_IMAGES.sh script (in 'Scripts/HPC_test_preprocessing/1.HPC_pre-processing'). (One image folder will be placed in the Data folder for each epoch)

run: 
```
sh MAIN_IMAGES.sh NREPL
```


3. Train a CNN and test it. Edit the for loop for selecting the training dataset from the 4 sorting options previously generated (inside the same folder than the dataset produced in the previous step).

run: 
```
sh MAIN_CNN.sh NREPL batch_size epochs filters kernel_size pooling_size
```
