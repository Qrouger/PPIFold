# PPIFold
Automatised pipeline for massive PPI prediction and figure creation.

This python package is use to generate PPI prediction on lot of system, without lost time in generate initial file and sort result.
It predict the better homo-oligomer for protein and better interface to interact with specific proteins.<br>
This allow to predict massive multimer complexe with lot of PPI.

## Requirements

- Python >=3.10
- AlphaPulldown 1.0.4 https://github.com/KosinskiLab/AlphaPulldown with singularity image (v0.4)
- SignalP5 https://services.healthtech.dtu.dk/services/SignalP-5.0/9-Downloads.php

## Installation

You need to install AlphaPulldown 1.0.4 with his database, SignalpP5 and python 3.10.

git clone https://github.com/Qrouger/PPIFold.git

## Pipeline

![pipeline](https://github.com/user-attachments/assets/21dc8eab-5322-4f00-942f-bdac4d723b72)


## Initial File

The initial file need to be Uniprot ID, separate by a comma and in a ".txt" file.

Ex : UniprotID1,UniprotID2,UniprotID3...

The file conf.txt need to be filed with all Paths.

Path_Uniprot_ID : Path and name of the UniprotID file.<br>
Path_AlphaFold_Data : Path of the AlphaFold data base (default on ./alphadata).<br>
Path_Singularity_Image : Path of the singularity image.<br>
Path_Pickle_Feature : Path of your feature folder (default on ./faeture).<br>

## Arguments<br>
python main.py --make_multimers Boolean --env_multimer String --env_feature String --max_aa Integer --use_signalP Boolean --org String

Optional arguments

--make_multimers This argument is True by default, if you just want make feature and MSA generation you have to set it on False <br>
--env_feature The name of the conda environment need to make feature with AlphaPulldown, set by default on None <br>
--env_multimer The name of the conda environment need to make multimer with AlphaPulldown, set by default on None <br>
--max_aa The maximum lenght of a model generable by your GPU (in function of Vram), set by default on 2400 (24 Go) <br>
--use_signalP Use SignalP if your proteins can be periplasmic, set by default on True <br>
--org If you use SignalP, you can select the organism (gram-, gram+, arch or euch), set by default on gram- <br>

Good use :<br>
Save all your pickles files on the same directory. 

## Result

This pipeline have a cutoff on PAE (10), iQ-score (35) and hiQ-score (50).


### Figures
**MSA depth<br>**
All aligned homologous sequences.

**Residue interaction table<br>**
Table of distance into two atoms of two chains.

**Distogram<br>**
Distance map between each atoms of each chains.

**Interaction network<br>**
**iQ-Score heatmap<br>**

### Other Files
**OOM_int.txt<br>**
A text file who contains too large interaction in function of --max_aa.<br>

**Bad_MSA.txt<br>**
A text file who contains proteins with MSA depth lower than 100 sequences.<br>
/!\Result for proteins with less than 100 sequences in MSA is not accurate for validate or invalidate predict PPI.<br>

**table.cyt<br>**
A file to make network manualy on Cytoscape.

Packages versions :<br>
numpy: v2.0.1<br>
matplotlib: v3.9.1<br>
Bio: v1.7.1<br>
networkx: v3.3<br>

