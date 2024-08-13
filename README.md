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

### Initial File

The initial file need to be Uniprot ID, separate by a comma and in a ".txt" file.

UniprotID1,UniprotID2,UniprotID3...

### Arguments

Mandatory

 --txt_name Name of the file who contains all Uniprot ID <br>
 --data_dir Path to the directory with all alphafold database <br>
 --dir_alpha_analysis Path to singularity image <br>

Optional

--make_multimers This argument is True by default, if you just want make feature you have to set it on False <br>
--env_feature The name of the conda environment need to make feature with AlphaPulldown, set by default on None <br>
--env_multimer The name of the conda environment need to make multimer with AlphaPulldown, set by default on None <br>
--max_aa The maximum lenght of a model generable by your GPU (in function of Vram), set by default on 2400 (24 Go) <br>
--use_signalP Use SignalP if your protéins can be periplasmic, set by default on True <br>

### Result

This pipeline have a cutoff on PAE, iQ-score and hiQ-score. An alert is set for proteins with MSA depth lower than 100 sequences.<br>
Result for proteins with less than 100 sequences in MSA is not accurate for validate and invalidate predict PPI.
