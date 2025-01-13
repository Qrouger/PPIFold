# PPIFold
Automated pipeline for massive PPI prediction and figure creation.

PPIFold is a tool for analyzing Protein-Protein Interactions from AlphaPullDown, with automated pre- and post-processing. It is used to generate PPI predictions for multiple systems without wasting time on generating initial files and sorting results. It predicts the best homo-oligomer for a protein and the best interface for interacting with specific proteins. This allows the prediction of massive multimeric complexes with numerous PPIs.

## Requirements
- AlphaFold data base
- Anaconda
- SignalP5
- Singularity and Singularity Image

## Installation
Installation of AlphaFold database :<br>
```bash
sudo apt install aria2
git clone https://github.com/deepmind/alphafold.git
cd alphafold
scripts/download_all_data.sh /<Directory></Directory> > download.log 2> download_all.log
```

Installation of SignalP5 (optional) :<br>

https://services.healthtech.dtu.dk/services/SignalP-5.0/9-Downloads.php<br>
```bash
tar -xvzf signalp-5.0b.Linux.tar.gz
cd signalp-5.0b/
cp bin/signalp /usr/local/bin
sudo cp -r lib/* /usr/local/lib
```

> [!NOTE]
> If you do not want to use SignalP, set --use_signalP to False.

Installation of Singularity :<br>

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

The initial file can be set using a UniProt ID, a FASTA sequence, or both.<br>
UniProt IDs need to be on the same line, separated by commas.<br>

Ex : <br>
UniprotID1,UniprotID2,UniprotID3...<br>

The FASTA sequence needs to start with ">", followed by the protein name. <br>

Ex : <br>
\>Name<br>
MFKRSGSLSLALMSSFCSSSLATPLSSAEFDHVARKCAPSVATSTLAAIAK<br>
VESRFDPLAIHDNTTGETLHWQDHTQATQVVRHRLDARHSLDVGLMQINSR<br>
NFSMLGLTPDGALKACPSLSAAANMLKSRYAGGETIDEKQIALRRAISAYN<br>
TGNFIRGFANGYVRKVETAAQSLVPALIEPPQDDHKALKSEDTWDVWGSYQ<br>
RRSQEDGVGGSIAPQPPDQDNGKSADDNQVLFDLY<br>

This file needs to be a ".txt" file.<br>

The conf.txt file needs to contains all path.

Path_Uniprot_ID : Path and name of the UniprotID file.<br>
Path_AlphaFold_Data : Path to the AlphaFold database (default on ./alphadata).<br>
Path_Singularity_Image : Path to the singularity image.<br>
Path_Pickle_Feature : Path to the feature folder (default on ./feature).<br>

## Arguments<br>
```bash
PPIFold --make_multimers Boolean --max_aa Integer --use_signalP Boolean --org String
```
Optional arguments

--make_multimers This argument is set on True by default, if you just want to make feature and MSA generation you have to set it on False <br>
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

**Shallow_MSA.txt<br>**
A text file who contains proteins with MSA depth lower than 100 sequences.<br>
> [!WARNING]
> Result for proteins with less than 100 sequences in MSA is not accurate for validate or invalidate predict PPI.<br>

**table.cyt<br>**
A file to make network manualy on Cytoscape.

**.pdb file<br>**
Model structure, color residue in interaction with the B-factor.

> [!WARNING]
> Two pipelines cannot be launched simultaneously on the same PC.
