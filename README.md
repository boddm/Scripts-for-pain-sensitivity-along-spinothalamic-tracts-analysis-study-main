# Spinothalamic Tract Microstructure Analysis

MATLAB scripts used in the study **“Spinothalamic Tract Microstructure Is Associated with Pain Sensitivity Across Modalities: A Combined Brain-Spinal Cord Diffusion Imaging Study”**

## Requirements

- MATLAB
- Statistics and Machine Learning Toolbox
- 8 GB RAM or more

The scripts have been tested with MATLAB R2025a-2026a on Windows, macOS, and Linux.

## Setup

Download this repository, open MATLAB in the project folder, and run:

```matlab
addpath(genpath('1_PLSC_Script'))
addpath(genpath('2_Prediction_Script'))
```

Sample datasets are provided in `0_Data/`.

## Run the analyses

### PLSC analysis

Edit the data settings in `1_PLSC_Script/myPLS_inputs.m`, then run:

```matlab
cd 1_PLSC_Script
myPLS_main
```

### PLSC + PLSR prediction

Edit the data and output paths in the main script, then run:

```matlab
cd '2_Prediction_Script/1_PLSC+PLSR/CPM_PLSC_PLSR'
CPM_PLSC_PLSR_real
```

Additional prediction and final-model scripts are available in:

```text
2_Prediction_Script/2_PLSC+Mean+ElasticNet+PLSR/
```

## Data format

Input data should be stored in a `.mat` file. The imaging and behavioral matrices must contain the same number of subjects.

```matlab
save('data.mat', 'brain_data', 'beh_data', 'diagnosis')
```

- `brain_data`: subjects × imaging features
- `beh_data`: subjects × behavioral variables
- `diagnosis`: optional group labels

Results are saved as `.mat` files in the output folder. Some scripts also display figures automatically.

## Citation

Please cite the related paper when using these scripts:

```text
[Citation will be added after publication]
```

## License and contact

MIT License\
Questions: <boddm123@gmail.com>
