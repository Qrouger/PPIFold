# PPIFold
Automatised pipeline for massive PPI prediction and figure creation.

This python package is use to generate PPI prediction on lot of system, without lost time in generate initial file and sort result. 
### Requirements

- Python >3.9
- AlphaPulldown 1.0.4 https://github.com/KosinskiLab/AlphaPulldown with singularity image (v0.4)
- SignalP5 https://services.healthtech.dtu.dk/services/SignalP-5.0/9-Downloads.php

### Installation

*Installation steps*

### Pipeline
 
![Pipeline](Pipeline.PNG)

### Arguments

Obligatory
 --txt_name Name of the file who contains all Uniprot ID
 --data_dir Path to the directory with all alphafold database
 --dir_alpha_analysis Path to singularity image

Optionnal
--make_multimers This argument is True by default, if you just want to make feature you have to set it on False
--env_feature The name of the conda environment need to make feature with AlphaPulldown, set by default on None
--env_multimer The name of the conda environment need to make multimer with AlphaPulldown, set by default on None
--max_aa The maximum lenght of a model generable by your GPU (in function of Vram), set by default on 2500 (24 Go)
--use_signalP Use SignalP if your protéins can be periplasmic, set by default on True
