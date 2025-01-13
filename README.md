# PPIFold
Automatised pipeline for massive PPI prediction and figure creation.

PPIFold it's a tool for analysis Protein-Protein Interaction from AlphaPullDown, with an automatised pre and post-processing.
It's uses to generate PPI prediction on lot of systems, without losing time in generating initial file and sort results.
It predicts the better homo-oligomer for protein and better interface to interact with specific proteins.
This allows to predict massive multimer complex with lot of PPI.

## Requirements
- AlphaFold data base
- Anaconda
- SignalP5
- Singularity and Singularity Image

## Installation
Installation AlphaFold database :<br>
```bash
sudo apt install aria2
git clone https://github.com/deepmind/alphafold.git
cd alphafold
scripts/download_all_data.sh /<Directory></Directory> > download.log 2> download_all.log
```

Install SignalP5 (optional) :<br>

https://services.healthtech.dtu.dk/services/SignalP-5.0/9-Downloads.php<br>
```bash
tar -xvzf signalp-5.0b.Linux.tar.gz
cd signalp-5.0b/
cp bin/signalp /usr/local/bin
sudo cp -r lib/* /usr/local/lib
```

> [!NOTE]
> If you don't want to use SignalP use --use_signalP False

Installation Singularity :<br>

https://docs.sylabs.io/guides/3.0/user-guide/installation.html#install-on-linux

Download Singularity image (score generation) :<br>

https://github.com/KosinskiLab/AlphaPulldown?tab=readme-ov-file#03-installation-for-the-downstream-analysis-tools

Installation PPIFold :<br>
```bash
conda create -n PPIFold -c omnia -c bioconda -c conda-forge python==3.11 openmm==8.0 pdbfixer==1.9 kalign2 hhsuite hmmer modelcif networkx
conda activate PPIFold
pip install PPIFold
pip install -U "jax[cuda12]"
```
## Pipeline

![Rouger_2024_Figure-1](https://github.com/user-attachments/assets/09bacfc4-103c-4910-aeea-2753e0ccca33)


## Initial File

Inital file can be set with Uniprot ID, sequence fasta or both.<br>
UniprotID need to be in same line separate by comma.<br>

Ex : <br>
UniprotID1,UniprotID2,UniprotID3...<br>

Fasta sequence need to start with ">" follow by the protein name.<br>

Ex : <br>
\>Name<br>
MFKRSGSLSLALMSSFCSSSLATPLSSAEFDHVARKCAPSVATSTLAAIAK<br>
VESRFDPLAIHDNTTGETLHWQDHTQATQVVRHRLDARHSLDVGLMQINSR<br>
NFSMLGLTPDGALKACPSLSAAANMLKSRYAGGETIDEKQIALRRAISAYN<br>
TGNFIRGFANGYVRKVETAAQSLVPALIEPPQDDHKALKSEDTWDVWGSYQ<br>
RRSQEDGVGGSIAPQPPDQDNGKSADDNQVLFDLY<br>

This file need to be a ".txt" file.<br>

The file conf.txt need to contains all Paths.

Path_Uniprot_ID : Path and name of the UniprotID file.<br>
Path_AlphaFold_Data : Path of the AlphaFold data base (default on ./alphadata).<br>
Path_Singularity_Image : Path of the singularity image.<br>
Path_Pickle_Feature : Path of your feature folder (default on ./feature).<br>

## Arguments<br>
```bash
PPIFold --make_multimers Boolean --max_aa Integer --use_signalP Boolean --org String
```
Optional arguments

--make_multimers This argument is True by default, if you just want make feature and MSA generation you have to set it on False <br>
--max_aa The maximum lenght of a model generable by your GPU (in function of Vram), set by default on 2000 (24 Go) <br>
--use_signalP Use SignalP if your proteins can be periplasmic, set by default on True <br>
--org If you use SignalP, you can select the organism (gram-, gram+, arch or euk), set by default on gram- <br>

>[!TIP]
> Save all your pickles files on the same directory. 

## Result

This pipeline have a cutoff on PAE (10), iQ-score (35) and hiQ-score (50).


### Figures
**MSA depth<br>**
All aligned homologous sequences for O50333.<br>
<img alt="MSA Depth" src="https://github.com/user-attachments/assets/d31c276a-ac54-4b22-a305-531d30e8c270" width="400"/><br>

Axe y is the number of homologous sequence, axe x is the posistions on the sequence. Colour represent the sequence identity.

**Residue interaction table<br>**
Table of distance into two atoms of O50331 and O5333.<br>
<img width="400" alt="residue_interaction_table" src="https://github.com/user-attachments/assets/ffdddb90-db4b-42ab-bdfe-41f09eae98f4" /><br>

Chains are differents proteins, two residue in contact are specified as well as their distances. Distances are calculated from center of mass the residues. Distance has a threshold of 10 and PAE of 5.

**Distogram<br>**
Distance map between each atoms of O50331 and O5333.<br>
<img width="400" alt="Distogram" src="https://github.com/user-attachments/assets/42476dbd-7c90-4a3c-b95f-84ee6c495b34" /><br>

Axes x and y are proteins in interaction, pixels into black square represent intra residue distances and pixels outside inter residue distances.
Colour represent distance in ångström, blue color represent a short distance between two residue and yellow a big distance.

**Interaction network<br>**
PPI network with iQ-score and homo-oligomers (hiQ-score) predictions.<br>

<img src="https://github.com/user-attachments/assets/8af56e6d-3548-452a-b720-4c4d6c2dac68" alt="interaction_network" width="400"/><br>
This network represent interaction into R388 proteins. Each interaction is represented by a line connecting two proteins, colored with the corresponding iQ-score. A loop on a protein indicates the better homo-oligomers with the higher hiQ-score.

**iQ-Score heatmap<br>**
Heatmap of iQ-score between each PPI.<br>
<img src="https://github.com/user-attachments/assets/cf1b8d62-45e8-41a7-a149-e4e466ba251c" alt="iQ_score_heatmap" width="400"/><br>
Colour represent iQ-score, better iQ-score is represent by a lighter color. 
The black boxes represent either bad PAE, homo-oligomer or too large total proteins lenght.

### Generated Files
**OOM_int.txt<br>**
A text file who contains too large interaction in function of --max_aa.<br>

**Bad_MSA.txt<br>**
A text file who contains proteins with MSA depth lower than 100 sequences.<br>
> [!WARNING]
> Result for proteins with less than 100 sequences in MSA is not accurate for validate or invalidate predict PPI.<br>

**table.cyt<br>**
A file to make network manualy on Cytoscape.

**.pdb file<br>**
Model structure, color residue in interaction with the B-factor.

> [!WARNING]
> Two pipelines cannot be launched simultaneously on the same PC.
