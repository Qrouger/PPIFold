# ColabFold - v1.5.5

For details of what was changed in v1.5, see [change log](https://github.com/sokrypton/ColabFold/wiki/v1.5.0)!

<p align="center"><img src="https://github.com/sokrypton/ColabFold/raw/main/.github/ColabFold_Marv_Logo.png" height="250"/></p>

### Making Protein folding accessible to all via Google Colab!

| Notebooks                                                                                                                                        | monomers | complexes | mmseqs2 | jackhmmer | templates |
| :----------------------------------------------------------------------------------------------------------------------------------------------- | -------- | --------- | ------- | --------- | --------- |
| [AlphaFold2_mmseqs2](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)                                    | Yes      | Yes       | Yes     | No        | Yes       |
| [AlphaFold2_batch](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/batch/AlphaFold2_batch.ipynb)                          | Yes      | Yes       | Yes     | No        | Yes       |
| [AlphaFold2](https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb) (from Deepmind)                    | Yes      | Yes       | No      | Yes       | No        |
| [relax_amber](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/beta/relax_amber.ipynb) (relax input structure)             |          |           |         |           |           |
| [ESMFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/ESMFold.ipynb)                                                  | Yes      | Maybe     | No      | No        | No        |
|                                                                                                                                                  |
| **BETA (in development) notebooks**                                                                                                              |          |           |         |           |           |
| [RoseTTAFold2](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/RoseTTAFold2.ipynb)                                        | Yes      | Yes       | Yes     | No        | WIP       |
| [OmegaFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/beta/omegafold.ipynb)                                         | Yes      | Maybe     | No      | No        | No        |
| [AlphaFold2_advanced_v2](https://colab.research.google.com/github/sokrypton/ColabDesign/blob/gamma/af/examples/predict.ipynb) (new experimental notebook)                  | Yes      | Yes       | Yes     | No        | Yes       |

Check the wiki page [old retired notebooks](https://github.com/sokrypton/ColabFold/wiki/Old-retired-notebooks) for unsupported notebooks.

### FAQ
- Where can I chat with other ColabFold users?
  - See our [Discord](https://discord.gg/gna8maru7d) channel!
- Can I use the models for **Molecular Replacement**?
  - Yes, but be **CAREFUL**, the bfactor column is populated with pLDDT confidence values (higher = better). Phenix.phaser expects a "real" bfactor, where (lower = better). See [post](https://twitter.com/cheshireminima/status/1423929241675120643) from Claudia Millán.
- What is the maximum length?
  - Limits depends on free GPU provided by Google-Colab `fingers-crossed`
  - For GPU: `Tesla T4` or `Tesla P100` with ~16G the max length is ~2000
  - For GPU: `Tesla K80` with ~12G the max length is ~1000
  - To check what GPU you got, open a new code cell and type `!nvidia-smi`
- Is it okay to use the MMseqs2 MSA server (`cf.run_mmseqs2`) on a local computer?
  - You can access the server from a local computer if you queries are serial from a single IP. Please do not use multiple computers to query the server.
- Where can I download the databases used by ColabFold?
  - The databases are available at [colabfold.mmseqs.com](https://colabfold.mmseqs.com)
- I want to render my own images of the predicted structures, how do I color by pLDDT?
  - In pymol for AlphaFold structures: `spectrum b, red_yellow_green_cyan_blue, minimum=50, maximum=90`
  - If you want to use AlphaFold Colours (credit: Konstantin Korotkov)
    ```python
    set_color n0, [0.051, 0.341, 0.827]
    set_color n1, [0.416, 0.796, 0.945]
    set_color n2, [0.996, 0.851, 0.212]
    set_color n3, [0.992, 0.490, 0.302]
    color n0, b < 100; color n1, b < 90
    color n2, b < 70;  color n3, b < 50
    ```
  - In pymol for RoseTTAFold structures: `spectrum b, red_yellow_green_cyan_blue, minimum=0.5, maximum=0.9`
- What is the difference between the AlphaFold2_advanced and AlphaFold2_mmseqs2 (_batch) notebook for complex prediction?
  - We currently have two different ways to predict protein complexes: (1) using the AlphaFold2 model with residue index jump and (2) using the AlphaFold2-multimer model. AlphaFold2_advanced supports (1) and AlphaFold2_mmseqs2 (_batch) (2).
- What is the difference between localcolabfold and the pip installable colabfold_batch?
  -  [LocalColabFold](https://github.com/YoshitakaMo/localcolabfold) is an installer script designed to make ColabFold functionality available on local users' machines. It supports wide range of operating systems, such as Windows 10 or later (using Windows Subsystem for Linux 2), macOS, and Linux.
- Is there a way to amber-relax structures without having to rerun alphafold/colabfold from scratch?
  - Yes, see this [notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/beta/relax_amber.ipynb).
- Where can I find the old notebooks that were previously developed and are now retired?
  - You can find the list of retired notebooks in the [old retired notebooks](https://github.com/sokrypton/ColabFold/wiki/Old-retired-notebooks) wiki page.
- Where can I find the history of MSA Server Databases used in ColabFold?
  - You can view the database version history on the [MSA Server Database History](https://github.com/sokrypton/ColabFold/wiki/MSA-Server-Database-History) wiki page.

### Running locally
For instructions on how to install ColabFold locally refer to [localcolabfold](https://github.com/YoshitakaMo/localcolabfold) or see our [wiki](https://github.com/sokrypton/ColabFold/wiki/Running-ColabFold-in-Docker) on how to run ColabFold within Docker.

### Generating MSAs for small scale local structure/complex predictions using the MSA server

When you pass a FASTA or CSV file containing your sequences to `colabfold_batch` it will automatically query the public MSA server to generate MSAs. You might want to split this into two steps for better GPU resource utilization:

```
# Query the MSA server and predict the structure on local GPU in one go:
colabfold_batch input_sequences.fasta out_dir

# Split querying MSA server and GPU predictions into two steps
colabfold_batch input_sequences.fasta out_dir --msa-only
colabfold_batch input_sequences.fasta out_dir
```

### Generating MSAs for large scale structure/complex predictions

First create a directory for the databases on a disk with sufficient storage (940GB (!)). Depending on where you are, this will take a couple of hours:

Note: [MMseqs2 `71dd32ec43e3ac4dabf111bbc4b124f1c66a85f1` (May 28, 2023)](https://github.com/soedinglab/MMseqs2/archive/71dd32ec43e3ac4dabf111bbc4b124f1c66a85f1.zip) is used to create the databases and perform sequece search in the ColabFold MSA server. Please use this version if you want to obtain the same MSAs as the server.

```shell
MMSEQS_NO_INDEX=1 ./setup_databases.sh /path/to/db_folder
```

If MMseqs2 is not installed in your `PATH`, add `--mmseqs <path to mmseqs>` to your `mmseqs` in `colabfold_search`:

```shell
# This needs a lot of CPU
colabfold_search --mmseqs /path/to/bin/mmseqs input_sequences.fasta /path/to/db_folder msas
# This needs a GPU
colabfold_batch msas predictions
```

This will create intermediate folder `msas` that contains all input multiple sequence alignments formated as a3m files and a `predictions` folder with all predicted pdb,json and png files.

The procedure above disables MMseqs2 preindexing of the various ColabFold databases by setting the `MMSEQS_NO_INDEX=1` environment variable before calling the database setup script. For most use-cases of `colabfold_search` precomputing the index is not required and might hurt search speed. The precomputed index is necessary for fast response times of the ColabFold server, where the whole database is permamently kept in memory. In any case the batch searches will require a machine with about 128GB RAM or, if the databases are to be kept permamently in RAM, with over 1TB RAM.

In some cases using precomputed database can still be useful. For the following cases, call the `setup_databases.sh` script without the `MMSEQS_NO_INDEX` environment variable:

(0) As mentioned above, if you want to set-up a server.

(1) If the precomputed index is stored on a very fast storage system (e.g., NVMe-SSDs) it might be faster to read the index from disk than computing in on the fly.  In this case, the search should be performed on the same machine that called `setup_databases.sh` since the precomputed index is created to fit within the given main memory size. Additionaly, pass the `--db-load-mode 0` option to make sure the database is read once from the storage system before use.

(2) Fast single query searches require the full index (the `.idx` files) to be kept in memory. This can be done with e.g. by using [vmtouch](https://github.com/hoytech/vmtouch). Thus, this type of search requires a machine with at least 768GB to 1TB RAM for the ColabfoldDB. If the index is present in memory, use the `--db-load-mode 2` parameter in `colabfold_search` to avoid index loading overhead.

If no index was created (`MMSEQS_NO_INDEX=1` was set), then `--db-load-mode` does not do anything and can be ignored.

### Tutorials & Presentations
- ColabFold Tutorial presented at the Boston Protein Design and Modeling Club. [[video]](https://www.youtube.com/watch?v=Rfw7thgGTwI) [[slides]](https://docs.google.com/presentation/d/1mnffk23ev2QMDzGZ5w1skXEadTe54l8-Uei6ACce8eI).

### Projects based on ColabFold or helpers

- [Run ColabFold on your local computer](https://github.com/YoshitakaMo/localcolabfold) by Yoshitaka Moriwaki
- [ColabFold/AlphaFold2 for protein structure predictions for Discoba species](https://github.com/zephyris/discoba_alphafold) by Richard John Wheeler
- [Cloud-based molecular simulations for everyone](https://github.com/pablo-arantes/Making-it-rain) by Pablo R. Arantes, Marcelo D. Polêto, Conrado Pedebos and Rodrigo Ligabue-Braun
- [getmoonbear is a webserver to predict protein structures](https://www.getmoonbear.com/AlphaFold2) by Stephanie Zhang and Neil Deshmukh
- [ColabFold/AlphaFold2 IDR complex prediction](https://github.com/normandavey/AlphaFold2-IDR-complex-prediction) by Balint Meszaros
- [ColabFold/AlphaFold2 (Phenix version) for macromolecular structure determination](https://colab.research.google.com/github/phenix-project/Colabs/blob/main/alphafold2/AlphaFold2.ipynb) by Tom Terwilliger
- [AlphaPickle: making AlphaFold2/ColabFold outputs interpretable](https://colab.research.google.com/github/mattarnoldbio/alphapickle/blob/main/AlphaPickle.ipynb) by Matt Arnold

### Acknowledgments
- We would like to thank the [RoseTTAFold](https://github.com/RosettaCommons/RoseTTAFold) and [AlphaFold](https://github.com/deepmind/alphafold) team for doing an excellent job open sourcing the software.
- Also credit to [David Koes](https://github.com/dkoes) for his awesome [py3Dmol](https://3dmol.csb.pitt.edu/) plugin, without whom these notebooks would be quite boring!
- A colab by Sergey Ovchinnikov (@sokrypton), Milot Mirdita (@milot_mirdita) and Martin Steinegger (@thesteinegger).

### How do I reference this work?

- Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S and Steinegger M. ColabFold: Making protein folding accessible to all. <br />
  Nature Methods (2022) doi: [10.1038/s41592-022-01488-1](https://www.nature.com/articles/s41592-022-01488-1)
- If you’re using **AlphaFold**, please also cite: <br />
  Jumper et al. "Highly accurate protein structure prediction with AlphaFold." <br />
  Nature (2021) doi: [10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2)
- If you’re using **AlphaFold-multimer**, please also cite: <br />
  Evans et al. "Protein complex prediction with AlphaFold-Multimer." <br />
  biorxiv (2021) doi: [10.1101/2021.10.04.463034v1](https://www.biorxiv.org/content/10.1101/2021.10.04.463034v1)
- If you are using **RoseTTAFold**, please also cite: <br />
  Minkyung et al. "Accurate prediction of protein structures and interactions using a three-track neural network." <br />
  Science (2021) doi: [10.1126/science.abj8754](https://doi.org/10.1126/science.abj8754)

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.5123296.svg)](https://doi.org/10.5281/zenodo.5123296)

-----------------
**OLD Updates**
```diff
  31Jul2023: 2023/07/31: The ColabFold MSA server is back to normal
             It was using older DB (UniRef30 2202/PDB70 220313) from 27th ~8:30 AM CEST to 31st ~11:10 AM CEST.
  27Jul2023: ColabFold MSA server issue:
             We are using the backup server with old databases
             (UniRef30 2202/PDB70 220313) starting from ~8:30 AM CEST until we resolve the issue.
             Resolved on 31Jul2023 ~11:10 CEST.
  12Jun2023: New databases! UniRef30 updated to 2302 and PDB to 230517.
             We now use PDB100 instead of PDB70 (see notes in the [main](https://colabfold.com) notebook).
  12Jun2023: We introduced a new default pairing strategy:
             Previously, for multimer predictions with more than 2 chains,
             we only pair if all sequences taxonomically match ("complete" pairing).
             The new default "greedy" strategy pairs any taxonomically matching subsets.
  30Apr2023: Amber is working again in our ColabFold Notebook
  29Apr2023: Amber is not working in our Notebook due to Colab update
  18Feb2023: v1.5.2 - fixing: fixing memory leak for large proteins
                    - fixing: --use_dropout (random seed was not changing between recycles)
  06Feb2023: v1.5.1 - fixing: --save-all/--save-recycles
  04Feb2023: v1.5.0 - ColabFold updated to use AlphaFold v2.3.1!
  03Jan2023: The MSA server's faulty hardware from 12/26 was replaced.
             There were intermittent failures on 12/26 and 1/3. Currently,
             there are no known issues. Let us know if you experience any.
  10Oct2022: Bugfix: random_seed was not being used for alphafold-multimer.
             Same structure was returned regardless of defined seed. This
             has been fixed!
  13Jul2022: We have set up a new ColabFold MSA server provided by Korean
             Bioinformation Center. It provides accelerated MSA generation,
             we updated the UniRef30 to 2022_02 and PDB/PDB70 to 220313.
  11Mar2022: We use in default AlphaFold-multimer-v2 weights for complex modeling.
             We also offer the old complex modes "AlphaFold-ptm" or "AlphaFold-multimer-v1"
  04Mar2022: ColabFold now uses a much more powerful server for MSAs and searches through the ColabFoldDB instead of BFD/MGnify.
             Please let us know if you observe any issues.
  26Jan2022: AlphaFold2_mmseqs2, AlphaFold2_batch and colabfold_batch's multimer complexes predictions are
             now in default reranked by iptmscore*0.8+ptmscore*0.2 instead of ptmscore
  16Aug2021: WARNING - MMseqs2 API is undergoing upgrade, you may see error messages.
  17Aug2021: If you see any errors, please report them.
  17Aug2021: We are still debugging the MSA generation procedure...
  20Aug2021: WARNING - MMseqs2 API is undergoing upgrade, you may see error messages.
             To avoid Google Colab from crashing, for large MSA we did -diff 1000 to get
             1K most diverse sequences. This caused some large MSA to degrade in quality,
             as sequences close to query were being merged to single representive.
             We are working on updating the server (today) to fix this, by making sure
             that both diverse and sequences close to query are included in the final MSA.
             We'll post update here when update is complete.
  21Aug2021  The MSA issues should now be resolved! Please report any errors you see.
             In short, to reduce MSA size we filter (qsc > 0.8, id > 0.95) and take 3K
             most diverse sequences at different qid (sequence identity to query) intervals
             and merge them. More specifically 3K sequences at qid at (0→0.2),(0.2→0.4),
             (0.4→0.6),(0.6→0.8) and (0.8→1). If you submitted your sequence between
             16Aug2021 and 20Aug2021, we recommend submitting again for best results!
  21Aug2021  The use_templates option in AlphaFold2_mmseqs2 is not properly working. We are
             working on fixing this. If you are not using templates, this does not affect the
             the results. Other notebooks that do not use_templates are unaffected.
  21Aug2021  The templates issue is resolved!
  11Nov2021  [AlphaFold2_mmseqs2] now uses Alphafold-multimer for complex (homo/hetero-oligomer) modeling.
             Use [AlphaFold2_advanced] notebook for the old complex prediction logic.
  11Nov2021  ColabFold can be installed locally using pip!
  14Nov2021  Template based predictions works again in the Alphafold2_mmseqs2 notebook.
  14Nov2021  WARNING "Single-sequence" mode in AlphaFold2_mmseqs2 and AlphaFold2_batch was broken
             starting 11Nov2021. The MMseqs2 MSA was being used regardless of selection.
  14Nov2021  "Single-sequence" mode is now fixed.
  20Nov2021  WARNING "AMBER" mode in AlphaFold2_mmseqs2 and AlphaFold2_batch was broken
             starting 11Nov2021. Unrelaxed proteins were returned instead.
  20Nov2021  "AMBER" is fixed thanks to Kevin Pan
```
-----------------
