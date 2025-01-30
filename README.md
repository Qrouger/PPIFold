# PPIFold
Automated pipeline for massive PPI prediction and figure creation.

PPIFold is a tool for analyzing Protein-Protein Interactions from [AlphaPulldown](https://github.com/KosinskiLab/AlphaPulldown#alphapulldown-version-200), with automated pre- and post-processing. It is used to generate PPI predictions for multiple systems without wasting time on generating initial files and sorting results. It predicts the best homo-oligomer for a protein and the best interface for interacting with specific proteins. This allows for the prediction of massive multimeric complexes with numerous PPIs.

## Requirements
- AlphaFold data base
- Anaconda
- SignalP5 (optional)
- Singularity and Singularity Image

## Installations
Installation of AlphaFold database :<br>
```bash
sudo apt install aria2
git clone https://github.com/deepmind/alphafold.git
cd alphafold
scripts/download_all_data.sh /<Directory></Directory> > download.log 2> download_all.log
```

SignalP5 installation (optional) :<br>

https://services.healthtech.dtu.dk/services/SignalP-5.0/9-Downloads.php<br>
```bash
tar -xvzf signalp-5.0b.Linux.tar.gz
cd signalp-5.0b/
cp bin/signalp /usr/local/bin
sudo cp -r lib/* /usr/local/lib
```

> [!NOTE]
> If you do not want to use SignalP, set --use_signalP to False.

Singularity installation :<br>

https://docs.sylabs.io/guides/3.0/user-guide/installation.html#install-on-linux

Download Singularity image (score generation) :<br>

https://github.com/KosinskiLab/AlphaPulldown?tab=readme-ov-file#03-installation-for-the-downstream-analysis-tools

PPIFold installation :<br>
```bash
conda create -n PPIFold -c omnia -c bioconda -c conda-forge python==3.11 openmm==8.0 pdbfixer==1.9 kalign2 hhsuite hmmer modelcif networkx
conda activate PPIFold
pip install PPIFold
pip install -U "jax[cuda12]"
```
## Pipeline

![Rouger_2024_Figure-1](https://github.com/user-attachments/assets/09bacfc4-103c-4910-aeea-2753e0ccca33)


## Initial Files

You need two intial file :

**test.txt**<br>
This file needs to be a ".txt" file.<br>
The initial file can be set up using UniProt IDs, FASTA sequences, or both.<br>
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

**conf.txt**<br>
The conf.txt file needs to contains all path.<br>

Path_Uniprot_ID : Path and name of the initial file.<br>
Path_AlphaFold_Data : Path to the AlphaFold database (default on ./alphadata).<br>
Path_Singularity_Image : Path and name of the singularity image.<br>
Path_Pickle_Feature : Path to your feature folder (default on ./feature).<br>

## Arguments<br>
To use PPIFold, simply run the PPIFold command in the folder containing conf.txt and test.txt.<br>
```bash
PPIFold --make_multimers Boolean --max_aa Integer --use_signalP Boolean --org String
```
Optional arguments

--make_multimers This argument is set to True by default. If you only want to generate features and MSA, you need to set it to False <br>
--max_aa The maximum length of a model that can be generated by your GPU (depending on VRAM), set to 2000 by default (24 GB) <br>
--use_signalP Use SignalP if your proteins can be periplasmic, set to True by default <br>
--org If you use SignalP, you can select the organism (gram-, gram+, arch or euk), set to Gram- by default <br>

>[!TIP]
> Save all your pickle files in the same directory. 

## Results

This pipeline have a cutoff on PAE (10), iQ-score (35) and hiQ-score (50).

### Figures
**MSA depth<br>**
All aligned homologous sequences for O50333.<br>
<img alt="MSA Depth" src="https://github.com/user-attachments/assets/d31c276a-ac54-4b22-a305-531d30e8c270" width="400"/><br>

The y-axis represents the number of homologous sequences, the x-axis represents the positions in the sequence. The color represents the sequence identity.

**Residue interaction table<br>**
Table of distance between two atoms of O50331 and O5333.<br>
<img width="400" alt="residue_interaction_table" src="https://github.com/user-attachments/assets/ffdddb90-db4b-42ab-bdfe-41f09eae98f4" /><br>

Chains represent different proteins. Two residues in contact are specified, along with their distances. Distances are calculated from the center of mass of the residues. The distance threshold is 10 angtroms, and the PAE is 5.

**Distogram<br>**
Distance map between each atom of O50331 and O5333.<br>
<img width="400" alt="Distogram" src="https://github.com/user-attachments/assets/42476dbd-7c90-4a3c-b95f-84ee6c495b34" /><br>

The x and y axes represent interacting proteins. Pixels inside the black squares represent intra-protein residue distances, while pixels outside represent inter-protein residue distances. The color represents the distance in angstroms: blue indicates a short distance between two residues, and yellow indicates a large distance.

**Interaction network<br>**
Protein-protein interaction network with iQ-score and homo-oligomers (hiQ-score) predictions.<br>

<img src="https://github.com/user-attachments/assets/8af56e6d-3548-452a-b720-4c4d6c2dac68" alt="interaction_network" width="400"/><br>
This network represents interactions between R388 proteins. Each interaction is represented by a line connecting two proteins, colored according to the corresponding iQ-score. A loop on a protein indicates the best homo-oligomers with the highest hiQ-score.

**iQ-Score heatmap<br>**
Heatmap of iQ-score between each PPI.<br>
<img src="https://github.com/user-attachments/assets/cf1b8d62-45e8-41a7-a149-e4e466ba251c" alt="iQ_score_heatmap" width="400"/><br>
Color represents the iQ-score, with a better iQ-score indicated by a lighter color. The black boxes represent either poor PAE, homo-oligomers, or overly large total protein length.

### Generated Files
**OOM_int.txt<br>**
A text file containing interactions that are too large, based on --max_aa.<br>

**Shallow_MSA.txt<br>**
A text file containing proteins with an MSA depth lower than 100 sequences.<br>
> [!WARNING]
> Results for proteins with fewer than 100 sequences in the MSA are not accurate for validating or invalidating predicted PPIs.

**table.cyt<br>**
A file for manually generating a network in Cytoscape.<br>

**_summary.signalp5<br>**
A file who resume signal peptides for all proteins.<br>


**.pdb file<br>**
Model structure, with residues colored according to their interaction. <br>
