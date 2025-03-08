""" All useful functions for PPIFold

    Author: Quentin Rouger
"""
import numpy as np
import pickle
import csv
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.special import softmax
import os
import networkx as nx
from math import *
import pandas as pd
import json
import gzip
import string
import seaborn
from Bio import PDB

from .File_proteins import *

def define_path() :
    """
    Extract all paths from conf.txt and store them in a dictionary.
    
    Parameters:
    ----------

    Returns:
    ----------
    path_dict : dict
    """
    path_dict = dict()
    with open("conf.txt", "r") as file :
        for lines in file :
            path_name = lines.split(":")[0].strip().strip("\n")
            path_dir = lines.split(":")[1].strip().strip("\n")
            path_dict[path_name] = path_dir
    for path_key in path_dict.keys() :
        if path_key == "Path_Singularity_Image" and ".sif" not in path_dict[path_key] :
            print("You need to specify the name of the Singularity image in Path_Singularity_Image.")
            exit()
        if len(path_dict[path_key]) == 0 :
            print (f'Path to {path_key} file is empty')
            if path_key == "Path_Uniprot_ID" or path_key == "Path_Singularity_Image" :
                exit()
            elif path_key == "Path_AlphaFold_Data" :
                print("set by default on ./alphadata")
                path_dict[path_key] = "./alphadata"
            elif path_key == "Path_Pickle_Feature" :
                print("set by default on ./feature")
                path_dict[path_key] = "./feature"
    return(path_dict)

def remove_SP (file, org) :
    """
    Create a new FASTA file without the signal peptide using SignalP.

    Parameters:
    ----------
    file : object of class File_proteins
    org : string

    Returns:
    ----------
    """
    final_file = str()
    SP_signal = 0
    prot_SP = dict()
    Prot_Signal_string = str()
    fasta_file = file.get_fasta_file()
    cmd = "signalp -fasta " + fasta_file + " -org " + org
    os.system(cmd)
    file_signalp = fasta_file.replace(".fasta","_summary.signalp5")
    with open(file_signalp,"r") as fh :
        for line in fh :
            new_line = line.split("\t")
            if new_line[1] != "OTHER" and new_line[0][0] != "#" :
                prot_SP[new_line[0]] = new_line[len(new_line)-1][11:13]
    new_fasta_dict = dict()
    with open(fasta_file, "r") as fa_file :
        for line2 in fa_file :
            new_line2 = line2
            if SP_signal == 0 and line2[0] != ">" :
                new_line2 = line2
                new_fasta_dict[save_key] = line2
            if int(SP_signal) > 0 :
                new_line2 = line2[int(SP_signal)-1:len(line2)]
                new_fasta_dict[save_key] = line2[int(SP_signal)+1:len(line2)]
                SP_signal = 0
            if line2[0] == ">" :
                save_key = line2[1:len(line2)-1]
                if str(line2[1:len(line2)-1]) in prot_SP.keys() :
                    SP_signal = prot_SP[line2[1:len(line2)-1]]
            final_file = final_file + new_line2
    file.set_proteins_sequence(new_fasta_dict)
    cmd2 = "rm " + fasta_file
    os.system(cmd2)
    with open(fasta_file, "w") as new_file2 :
        new_file2.write(final_file)

def create_feature (file, data_dir, Path_Pickle_Feature, mmseq) :
    """
    Launch command to generate features.

    Parameters:
    ----------
    file : object of class File_proteins
    data_dir : string
    Path_Pickle_Feature : string
    mmseq : boolean
    
    Returns:
    ----------
    """
    fasta_file = file.get_fasta_file()
    cmd = f"create_individual_features.py --fasta_paths=./{fasta_file} \--data_dir={data_dir} \--save_msa_files=True \--output_dir={Path_Pickle_Feature} \--max_template_date=2024-05-02 \--skip_existing=True \--use_mmseqs2={mmseq}"
    os.system(cmd)

def Make_all_MSA_coverage (file, Path_Pickle_Feature) :
    """
    Generating MSA coverage for all proteins.

    Parameters:
    ----------
    file : object of class File_proteins
    Path_Pickle_Feature : string
        
    Returns:
    ----------
    """
    old_proteins = file.get_proteins()
    new_proteins = file.get_new_pickle() 
    shallow_MSA = str()
    for prot in new_proteins :
        pre_feature_dict = pickle.load(open(f'{Path_Pickle_Feature}/{prot}.pkl','rb'))
        feature_dict = pre_feature_dict.feature_dict
        msa = feature_dict['msa']
        seqid = (np.array(msa[0] == msa).mean(-1))
        seqid_sort = seqid.argsort()
        non_gaps = (msa != 21).astype(float)
        non_gaps[non_gaps == 0] = np.nan
        final = non_gaps[seqid_sort] * seqid[seqid_sort, None]
        plt.figure(figsize=(14, 4), dpi=100)
        plt.subplot(1, 2, 1)
        plt.title(f"Sequence coverage ({prot})")
        plt.imshow(final, interpolation='nearest', aspect='auto', cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
        plt.plot((msa != 21).sum(0), color='black')
        plt.xlim(-0.5, msa.shape[1] - 0.5)
        plt.ylim(-0.5, msa.shape[0] - 0.5)
        plt.colorbar(label="Sequence identity to query", )
        plt.xlabel("Positions")
        plt.ylabel("Sequences")
        plt.savefig(f"{Path_Pickle_Feature}/{prot+('_' if prot else '')}coverage.pdf")
        plt.close()
    for prot in old_proteins : #just write shallow_MSA.txt
        pre_feature_dict = pickle.load(open(f'{Path_Pickle_Feature}/{prot}.pkl','rb'))
        feature_dict = pre_feature_dict.feature_dict
        msa = feature_dict['msa']
        if len(msa) <= 100 :
            shallow_MSA += prot + " : " + str(len(msa)) + " sequences\n"
    with open("shallow_MSA.txt", "w") as MSA_file :
        MSA_file.write(shallow_MSA)

def generate_APD_script (file, max_aa) :
    """
    Write two local scripts to use AlphaPullDown. These scripts should be written based on the maximum number of amino acids.

    Parameters:
    ----------
    file : object of class File_proteins
    max_aa : integer

    Returns:
    ----------
    """
    all_vs_all_script = str()
    homo_oligo_script = str()
    OOM_int = str()
    proteins = file.get_proteins()
    lenght_prot = file.get_lenght_prot()
    for index_protein in range(len(proteins)) :
        lenght = lenght_prot[proteins[index_protein]]
        for index2_protein in range(index_protein+1,len(proteins)) :
            int_lenght = lenght + lenght_prot[proteins[index2_protein]]
            if int_lenght >= max_aa :
                OOM_int = OOM_int + proteins[index_protein] + ";" +  proteins[index2_protein]+ "\n"
            elif os.path.exists(f"./result_all_vs_all/{proteins[index_protein]}_and_{proteins[index2_protein]}/ranked_0.pdb") == False and os.path.exists(f"./result_all_vs_all/{proteins[index2_protein]}_and_{proteins[index_protein]}/ranked_0.pdb") == False: #make interaction if doesn't exist and is not too long
                all_vs_all_script = all_vs_all_script + proteins[index_protein] + ";" +  proteins[index2_protein]+ "\n"
            else :
                pass
        lenght_homo = lenght
        homo_dir = proteins[index_protein]
        for nbr_homo in range(2,21) :
           lenght_homo += lenght
           homo_dir += "_and_" + proteins[index_protein]
           if lenght_homo >= max_aa :
               OOM_int = OOM_int + proteins[index_protein] + "," + str(nbr_homo) + "\n"
           elif os.path.exists(f"./result_homo_oligo/{homo_dir}/ranked_0.pdb") == False and len(homo_dir) < 4096 : #homo_oligo is too long to create a directory
               homo_oligo_script = homo_oligo_script + proteins[index_protein] + "," + str(nbr_homo) + "\n"
           else :
               pass
    with open("homo_oligo.txt", "w") as homo_file:
       homo_file.write(homo_oligo_script)
    with open("all_vs_all.txt", "w") as all_file:
       all_file.write(all_vs_all_script)
    with open("OOM_int.txt", "w") as OOM_file :
       OOM_file.write(OOM_int)

### Generating Multimers

def Make_all_vs_all (data_dir, Path_Pickle_Feature) :
    """
    Use Alphapulldown script to generate all versus all interactions.

    Parameters:
    ----------
    data_dir : string
    Path_Pickle_Feature : string

    Returns:
    ----------
    """
    os.environ['XLA_PYTHON_CLIENT_PREALLOCATE'] = 'false'
    os.environ['TF_FORCE_UNIFIED_MEMORY'] = 'true'
    os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '3.2'
    os.environ['XLA_FLAGS'] = '--xla_gpu_enable_triton_gemm=false'
    cmd1 =f"run_multimer_jobs.py --mode=custom \--num_cycle=3 \--num_predictions_per_model=1 \--compress_result_pickles=True \--output_path=./result_all_vs_all \--data_dir={data_dir} \--protein_lists=all_vs_all.txt \--monomer_objects_dir={Path_Pickle_Feature} \--remove_keys_from_pickles=False"
    os.system(cmd1)

def add_iQ_score (dir_alpha) :
    """
    Generate iQ_score for all interactions.

    Parameters:
    ----------
    dir_alpha : string

    Returns:
    ----------
    """    
    if os.path.isdir("./result_all_vs_all") == True :
       cmd4 = f"singularity exec --no-home --bind result_all_vs_all:/mnt {dir_alpha} run_get_good_pae.sh --output_dir=/mnt --cutoff=10"
       os.system(cmd4)
       with open("result_all_vs_all/predictions_with_good_interpae.csv", "r") as file1 :
          reader = csv.DictReader(file1)
          all_lines = "jobs,pi_score,iptm_ptm,pDockQ,iQ_score\n"
          for row in reader :
             job = row['jobs']
             if '_and_' in job and row['pi_score'] != 'No interface detected' :
                iQ_score = ((float(row['pi_score'])+2.63)/5.26)*40+float(row['iptm_ptm'])*30+float(row['mpDockQ/pDockQ'])*30
                line =f'{row["jobs"]},{row["pi_score"]},{row["iptm_ptm"]},{row["mpDockQ/pDockQ"]},{str(iQ_score)}\n'
                all_lines = all_lines + line
       with open("result_all_vs_all/new_predictions_with_good_interpae.csv", "w") as file2 :
          file2.write(all_lines)
    else : #allow to do PPIFold for only one protein
       print("all_vs_all directory is empty")
        
def create_out_fig (file) :
    """
    Generate result figure for validate interaction (iQ_score) and better interaction (hiQ_score).

    Parameters:
    ----------
    file : object of class File_proteins

    Returns:
    ----------
    """
    iQ_score_dict = file.get_iQ_score_dict()
    for interaction in iQ_score_dict.keys() :
        if float(iQ_score_dict[interaction]) >= 50 : #Plot figure of interest just for interesting interactions
            job1 = interaction[0] + "_and_" + interaction[1]
            plot_Distogram("./result_all_vs_all/" + job1) #need distogram key in pickle file
            make_table_res_int(file, "./result_all_vs_all/" + job1)
    hiQ_score_dict = file.get_hiQ_score_dict()
    for homo_oligo in hiQ_score_dict.keys() :
        if float(hiQ_score_dict[homo_oligo][0]) >= 50 :
            job2 = homo_oligo #job2 = homo_oligo + "_homo_" + str(hiQ_score_dict[homo_oligo][1]) + "er" # wait AFPD homo release
            for count in range(1,hiQ_score_dict[homo_oligo][1]) :
                job2 += "_and_" + homo_oligo
            plot_Distogram("./result_homo_oligo/" + job2) #need distogram key in pickle file
            make_table_res_int(file, "./result_homo_oligo/" + job2)

def make_table_res_int (file, path_int) :
    """
    Generate a table of residues in interactions.

    Parameters:
    ----------
    file : object of class File_proteins
    path_int : string

    Returns:
    ----------
    """
    ranking_results = json.load(open(os.path.join(f'{path_int}/ranking_debug.json')))
    best_model = ranking_results["order"][0]
    
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', path_int + "/ranked_0.pdb")
    dict_int = dict()
    int_already_know = dict()
    proteins = path_int.split('/')[2].split('_and_')
    color_res = dict()
    color_res[proteins[0]] = set()
    color_res[proteins[1]] = set()
    atom_possible_contact = ["C","CA","CB"] #["O","OH","NH2","NH1","OG","NE2","ND2","NZ","NE","N","OE1","OE2","OD2","OG1"] #hydrogen bond
    with open(os.path.join(f'{path_int}/result_{best_model}.pkl.gz'), 'rb') as gz_file :
       pickle_dict = pickle.load(gzip.open(gz_file))
       if "distogram" not in pickle_dict.keys() or "predicted_aligned_error" not in pickle_dict.keys() :
          for model in structure:
             list_chain = model.get_list()
             max_chain = len(list_chain)
             if max_chain > 3 : #for homo-oligo bigger than 3 make interface for third first chains with all other
                max_chain = 3
             for index1 in range(0,max_chain) :
                chain1 = list_chain[index1] #B and C
                for residue1 in chain1 : #number of residue (len of the chain)
                   for atom1 in residue1 : #type of the atom
                      if atom1.get_id() in atom_possible_contact :
                          for index2 in range(index1+1,len(list_chain)) :
                              chain2 = list_chain[index2]
                              for residue2 in chain2 :
                                  for atom2 in residue2 :
                                      if atom2.get_id() in atom_possible_contact :
                                          distance = atom1 - atom2
                                          if distance <= 6 : #and atom1.bfactor >=70 and atom2.bfactor >= 70 : #filtered on pLDDT and distance, be stringent to avoid false residue interaction (or maybe use PAE ?)
                                              res_int = chain1.get_id()+":"+residue1.get_resname()+" "+str(residue1.get_id()[1])," "+chain2.get_id()+":"+residue2.get_resname()+" "+str(residue2.get_id()[1])
                                              res_num = str(residue1.get_id()[1]), str(residue2.get_id()[1])
                                              if chain1.get_id()+chain2.get_id() in dict_int.keys() : #to make different table for different interaction
                                                  if res_int in int_already_know.keys() and int_already_know[res_int] > str(distance) :
                                                      dict_int[chain1.get_id()+chain2.get_id()].remove([res_int[0][2:5]+":"+res_num[0]," "+res_int[1][3:6]+":"+res_num[1]," "+str(int_already_know[res_int])])
                                                      dict_int[chain1.get_id()+chain2.get_id()].append([res_int[0][2:5]+":"+res_num[0]," "+res_int[1][3:6]+":"+res_num[1]," "+str(distance)])
                                                      int_already_know[res_int] = str(distance)
                                                      color_res[proteins[0]].add(res_num[0])
                                                      color_res[proteins[1]].add(res_num[1])
                                                  elif res_int in int_already_know.keys() and int_already_know[res_int] < str(distance) : #skip double interaction with differents atoms
                                                      pass
                                                  else :
                                                      dict_int[chain1.get_id()+chain2.get_id()].append([res_int[0][2:5]+":"+res_num[0]," "+res_int[1][3:6]+":"+res_num[1]," "+str(distance)])
                                                      int_already_know[res_int] = str(distance)
                                                      color_res[proteins[0]].add(res_num[0])
                                                      color_res[proteins[1]].add(res_num[1])
                                              else :
                                                  dict_int[chain1.get_id()+chain2.get_id()] = [["Chain "+chain1.get_id()," Chain "+chain2.get_id()," Distance Ä"]]
                                                  dict_int[chain1.get_id()+chain2.get_id()].append([res_int[0][2:5]+":"+res_num[0]," "+res_int[1][3:6]+":"+res_num[1]," "+str(distance)])
                                                  int_already_know[res_int] = str(distance)
                                                  color_res[proteins[0]].add(res_num[0])
                                                  color_res[proteins[1]].add(res_num[1])
                                          else :
                                              pass
       else : #last version of APD
          lenght_prot = file.get_lenght_prot()
          seq_prot = file.get_proteins_sequence()
          names = path_int.split("/")[2].split("_and_")
          chains = path_int.split("/")[2]
          dict_int = dict()
          color_res = dict()
          color_res[names[0]] = set()
          color_res[names[1]] = set()
          with open(os.path.join(f'{path_int}/result_{best_model}.pkl.gz'), 'rb') as gz_file :
             pickle_dict = pickle.load(gzip.open(gz_file))
             pae_mtx = pickle_dict['predicted_aligned_error']#take PAE
             bin_edges = pickle_dict["distogram"]["bin_edges"]#take distogram for distance
             bin_edges = np.insert(bin_edges, 0, 0)
             distogram_softmax = softmax(pickle_dict["distogram"]["logits"], axis=2)
             dist = np.sum(np.multiply(distogram_softmax, bin_edges), axis=2) #center of the residue
             dict_int[chains] = [[names[0]," "+names[1]," Distance Ä"," PAE score"]]
             for line in range(lenght_prot[names[0]],len(dist)) :
                hori_index = 0
                for distance in dist[line] :
                   hori_index += 1
                   if hori_index < lenght_prot[names[0]] :
                      if distance <= 10 :
                         if pae_mtx[line][hori_index] < 5 :
                             residue1 = seq_prot[names[0]][hori_index+1]
                             residue2 = seq_prot[names[1]][line-lenght_prot[names[0]]+1]
                             dict_int[chains].append([residue1+":"+str(hori_index+1)," "+residue2+":"+str(line-lenght_prot[names[0]]+1)," "+str(distance), " "+str(pae_mtx[line][hori_index])]) #+1 to match with pdb model
                             color_res[names[0]].add(str(hori_index+1))
                             color_res[names[1]].add(str(line-lenght_prot[names[0]]+1))
#    file.define_interface(dict_interface[chains],names) #update interaction interface
#    color_int_residues(path_int,color_res,names) #color residue in interaction on the pdb
#    fileout = chains+"_res_int.csv"
#    np_table = np.array(dict_interface[chains])
#    with open(f"{path_int}/"+fileout, "w", newline="") as file :
#        mywriter = csv.writer(file, delimiter=",")
#        mywriter.writerows(np_table) 
    
    index_homo = 1
    int_names = path_int.split("/")[2]
    residues_at_interface = dict()
    residues_at_interface[int_names] = []
    for chains in dict_int.keys() :
        index_homo += 1
        fileout = chains+"_res_int.csv"
        np_table = np.array(dict_int[chains])
        with open(f"{path_int}/"+fileout, "w", newline="") as csv_table :
             mywriter = csv.writer(csv_table, delimiter=",")
             mywriter.writerows(np_table)
        del dict_int[chains][0] #delete title of each col
        for interaction in dict_int[chains] :
            if interaction not in residues_at_interface[int_names] :
                residues_at_interface[int_names].append(interaction)
    print("Write residue table")
    names = path_int.split("/")[2].split("_and_")
    if residues_at_interface[int_names] != [] : #can arrive if it don't fin atom with distance < 10 or PAE < 10
       file.define_interface(residues_at_interface[int_names],names) #update interaction interface
       color_int_residues(path_int,color_res,proteins) #color residue in interaction on the pdb

def plot_Distogram (job) :
    """
    Generate distogram, only for best models.

    Parameters:
    ----------
    job : string
    
    Returns:
    ----------
    """
    ranking_results = json.load(open(os.path.join(f'{job}/ranking_debug.json')))
    best_model = ranking_results["order"][0]
    with open(os.path.join(f'{job}/result_{best_model}.pkl.gz'), 'rb') as gz_file :
        results = pickle.load(gzip.open(gz_file))
        if "distogram" in results.keys() : #avoid error from APD release 
           bin_edges = results["distogram"]["bin_edges"]
           bin_edges = np.insert(bin_edges, 0, 0)
           distogram_softmax = softmax(results["distogram"]["logits"], axis=2)
           dist = np.sum(np.multiply(distogram_softmax, bin_edges), axis=2)
           np.savetxt(f"{job}/result_{best_model}.pkl.dmap", dist)
           lenght_list = []
           for seq in results["seqs"] :
              lenght_list.append(len(seq))
           print("make Distogram")
           initial_lenght = 0
           fig, ax = plt.subplots()
           d = ax.imshow(dist)
           plt.colorbar(d, ax=ax, fraction=0.046, pad=0.04)
           ax.title.set_text("Distance map")
           for index in range(len(lenght_list)-1) :
              initial_lenght += lenght_list[index]
              ax.axhline(initial_lenght, color="black", linewidth=1.5)
              ax.axvline(initial_lenght, color="black", linewidth=1.5)
           plt.savefig(f"{job}/result_{best_model}.dmap.png", dpi=600)
           plt.close()
            
def Make_homo_oligo (data_dir, Path_Pickle_Feature) :
    """
    Use Alphapulldown script to generate all homo-oligomer.

    Parameters:
    ----------
    data_dir : string
    Path_Pickle_Feature : string
        
    Returns:
    ----------
    """
    os.environ['XLA_PYTHON_CLIENT_PREALLOCATE'] = 'false'
    os.environ['TF_FORCE_UNIFIED_MEMORY'] = 'true'
    os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '3.2'
    os.environ['XLA_FLAGS'] = '--xla_gpu_enable_triton_gemm=false'
    cmd1 =f"run_multimer_jobs.py --mode=custom \--num_cycle=3 \--num_predictions_per_model=1 \--compress_result_pickles=True \--output_path=./result_homo_oligo \--data_dir={data_dir} \--protein_lists=homo_oligo.txt \--monomer_objects_dir={Path_Pickle_Feature} \--remove_keys_from_pickles=False"
    os.system(cmd1)

def add_hiQ_score (dir_alpha) :
    """
    Generate a score table for all homo-oligomer.

    Parameters:
    ----------
    dir_alpha : string
        
    Returns:
    ----------
    """
    cmd4 = f"singularity exec --no-home --bind result_homo_oligo:/mnt {dir_alpha} run_get_good_pae.sh --output_dir=/mnt --cutoff=10"
    os.system(cmd4)
    with open("./result_homo_oligo/predictions_with_good_interpae.csv", "r") as file1 :
        reader = csv.DictReader(file1)
        all_lines = "jobs,pi_score,iptm_ptm,hiQ_score\n"
        all_homo = dict()
        save_pi_score = dict()
        for row in reader :
            job = row['jobs']
            #if 'homo' in job and row['pi_score'] != 'No interface detected' : #need AFPD release with homo-oligo ####_homo_2er
            if row['pi_score'] != 'No interface detected' :
                if job not in all_homo.keys() :
                    all_homo[job] = (row['pi_score'],1,row)
                    save_pi_score[job] = [float(row['pi_score'])]
                else :
                    save_pi_score[job].append(float(row['pi_score']))
                    sum_pi_score = float(all_homo[job][0]) + float(row['pi_score'])
                    sum_int = all_homo[job][1] + 1
                    all_homo[job] = (sum_pi_score,sum_int,row)
    for key in all_homo.keys() :
        row = all_homo[key][2]
        #number_oligo = row["jobs"].split("_")[2].replace("er","") #wait AFPD release homo_oligo
        number_oligo = len(row["jobs"].split("_and_"))
        if len(save_pi_score[key]) > int(number_oligo) : #if model have more interface than number of homo-oligomerization
            new_sum_pi_score = 0
            save_pi_score[key].sort(reverse=True)
            for index in range(0,int(number_oligo)) :
                new_sum_pi_score += save_pi_score[key][index]
                hiQ_score = (((float(new_sum_pi_score)/int(number_oligo))+2.63)/5.26)*60+float(row['iptm_ptm'])*40 #cause iptm_ptm are always same for each homo of same protein
            line =f'{key},{str(float(new_sum_pi_score)/int(number_oligo))},{row["iptm_ptm"]},{str(hiQ_score)}\n'
            all_lines += line
        else :
            hiQ_score = (((float(all_homo[key][0])/all_homo[key][1])+2.63)/5.26)*60+float(row['iptm_ptm'])*40 #cause iptm_ptm is always same for each homo of same protein
            line =f'{key},{str(float(all_homo[key][0])/all_homo[key][1])},{row["iptm_ptm"]},{str(hiQ_score)}\n'
            all_lines += line
    with open("./result_homo_oligo/new_predictions_with_good_interpae.csv", "w") as file2 :
        file2.write(all_lines)

def generate_interaction_network (file) :
    """
    Generate interaction network, colored by iQ_score.
   
    Parameters:
    ----------
    file : object of File_proteins class

    Returns:
    ----------
    """
    valid_interactions = list()
    iQ_score_dict = file.get_iQ_score_dict()
    for interactions in iQ_score_dict.keys() :
        names = [interactions[0], interactions[1]]
        if names not in [x[0] for x in valid_interactions] and float(iQ_score_dict[interactions]) >= 50 :
            valid_interactions.append([names, float(iQ_score_dict[interactions])])
    hiQ_score_dict = file.get_hiQ_score_dict()
    for homo_oligomer in hiQ_score_dict.keys() :
        if float(hiQ_score_dict[homo_oligomer][0]) >= 50 :
            valid_interactions.append([[homo_oligomer,homo_oligomer], hiQ_score_dict[homo_oligomer][1]])
    int_graph = nx.Graph()
    list_inter_score = list()
    prots = set()
    dict_name = file.get_names()
    with open('table.cyt', 'w') as f :
        f.write(('source,targer,interaction,score\n'))
        for inter, score in valid_interactions :
            if inter[0] in dict_name.keys() :
               inter0 = inter[0]+f"({dict_name[inter[0]]})" #set uniprotID with the name of protein
            else :
               inter0 = inter[0]
            if inter[1] in dict_name.keys() :
               inter1 = inter[1]+f"({dict_name[inter[1]]})"
            else :
               inter1 = inter[1]
            f.write(f'{inter0},{inter1},pp,{round(score,2)}\n')
            prots.add(inter0)
            prots.add(inter1)
            list_inter_score.append((inter0,inter1,round(score,2)))
    prots = list(prots)
    int_graph.add_nodes_from(prots)
    int_graph.add_weighted_edges_from(list_inter_score)
    fig, ax = plt.subplots(figsize=(10, 6), dpi = 300)
    pos = nx.spring_layout(int_graph, k = len(prots)+50, scale = 3, seed = 8)
    nx.draw_networkx_nodes(int_graph,pos)
    nx.draw_networkx_edges(int_graph,pos, node_size = 100, edgelist = int_graph.edges, width = 1, style = "solid")
    nx.draw_networkx_labels(int_graph, pos, font_size = 10, font_family = "sans-serif", hide_ticks = 'True')
    edge_labels = nx.get_edge_attributes(int_graph, "weight")
    homo_label_dict = dict()
    color_label = dict()
    for prot_int in edge_labels.keys() :
        if len(str(edge_labels[prot_int])) < 3 :
            homo_label_dict[prot_int] = edge_labels[prot_int]
        else :
            color_label[prot_int] = edge_labels[prot_int]
    selected_weights = {edge: edge_labels[edge] for edge in color_label.keys()}
    if len(selected_weights) != 0 :
       norm = mcolors.Normalize(vmin=min(selected_weights.values()), vmax=max(selected_weights.values()))
       cmap = plt.cm.coolwarm
       edge_colors = {}
       for edge in int_graph.edges():
          if edge in color_label.keys() or (edge[1], edge[0]) in color_label.keys():
             edge_colors[edge] = cmap(norm(edge_labels[edge]))
          else:
             edge_colors[edge] = 'gray'
       edge_colors_list = [edge_colors[edge] for edge in int_graph.edges()]
       nx.draw_networkx_edge_labels(int_graph, pos, homo_label_dict, verticalalignment = 'bottom')
       nx.draw(int_graph, pos, edge_color=edge_colors_list, node_color='lightblue', node_size=500, ax=ax)
       sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
       sm.set_array([])
       plt.colorbar(sm, ax=ax, label="iQ_score")
       plt.savefig("network.png")
       plt.close()
    else :
        print("no PPI interactions for network")

def generate_heatmap (file) :
    """
    Generate a heatmap of interaction scores between all proteins.
   
    Parameters:
    ----------
    file : object of File_proteins class

    Returns:
    ----------
    """
    iQ_score_data_matrix = list()
    iptm_ptm_data_matrix = list()
    iQ_score_dict = file.get_iQ_score_dict()
    proteins_list = file.get_proteins()
    proteins_name = file.get_names()
    index_prot = list()
    iptm_ptm_dict = dict()
    with open("result_all_vs_all/new_predictions_with_good_interpae.csv", "r") as file1 :
        reader = csv.DictReader(file1)
        for row in reader :
            iptm_ptm_dict[row['jobs']] = row["iptm_ptm"]
    for protein in proteins_list :
        if protein in proteins_name.keys() :
          index_prot.append(protein+"_"+proteins_name[protein])
        else :
          index_prot.append(protein)
    for protein1 in proteins_list :
        iQ_score_line = list()
        iptm_ptm_line = list()
        for protein2 in proteins_list :
            if protein1 == protein2 :
                iQ_score_line.append(0)
                iptm_ptm_line.append(0)
            else :
                if (protein1,protein2) in iQ_score_dict.keys() :
                    iQ_score_line.append(float(iQ_score_dict[(protein1,protein2)]))
                    iptm_ptm_line.append(float(iptm_ptm_dict[protein1+"_and_"+protein2]))
                elif (protein2,protein1) in iQ_score_dict.keys() :
                    iQ_score_line.append(float(iQ_score_dict[(protein2,protein1)]))
                    iptm_ptm_line.append(float(iptm_ptm_dict[protein2+"_and_"+protein1]))
                else :
                    iQ_score_line.append(0)
                    iptm_ptm_line.append(0)
                    print (protein1 + " and " + protein2 + " are not in score table or have bad inter PAE")
        iQ_score_data_matrix.append(iQ_score_line)
        iptm_ptm_data_matrix.append(iptm_ptm_line)
    iQ_score_complet_matrix = pd.DataFrame(iQ_score_data_matrix,index = index_prot, columns = index_prot)
    iptm_ptm_complet_matrix = pd.DataFrame(iptm_ptm_data_matrix,index = index_prot, columns = index_prot)
    ax2 = seaborn.heatmap(iQ_score_complet_matrix, cbar_kws = {'label' : 'iQ_score'}, xticklabels=True, yticklabels=True)
    ax2.figure.tight_layout()
    plt.savefig("iQ_score_heatmap.png")
    plt.close()
    ax2 = seaborn.heatmap(iptm_ptm_complet_matrix, cbar_kws = {'label' : 'iptm_ptm'}, xticklabels=True, yticklabels=True)
    ax2.figure.tight_layout()
    plt.savefig("iptm_ptm_heatmap.png")
    plt.close()

def cluster_interface (file) :
    """
    Compare interface in function of smaller interface and classify them with a letter representing the interface.
   
    Parameters:
    ----------
    file : object of File_proteins class

    Returns:
    ----------
    interface_dict : dict
    """
    alphabet = string.ascii_lowercase
    interface_dict = file.get_interface_dict()
    for proteins in interface_dict.keys() :
        already_inter = list()
        interface_dict[proteins] = sorted(interface_dict[proteins], key=lambda x : len(x)) #sorted all interface in function of number of residues
        for interface1 in range(len(interface_dict[proteins])) :
            if interface1 == 0 : #if it's the first interface, define a
                interface_dict[proteins][interface1].insert(0,alphabet[0])
                already_inter.append(alphabet[0])
            for interface2 in range(interface1+1,len(interface_dict[proteins])) :
                list_inter = list(set(interface_dict[proteins][interface1]).intersection(set(interface_dict[proteins][interface2])))
                simi_inter = len(list_inter)/(len(set(interface_dict[proteins][interface1]).union(set(interface_dict[proteins][interface2])))-3) #indice jaccard # -3 just to remove interface 'a' and uniprotID from .union()
                if simi_inter < 0.35: #create a new interface
                    if interface_dict[proteins][interface2][0] in already_inter : #Don't create new interface if it already has one
                        pass
                    else :
                        interface_dict[proteins][interface2].insert(0,alphabet[interface2])
                        already_inter.append(alphabet[interface2])
                else : #if interfaces got more than 0.50 of same residues, it's the same interface
                    interface_dict[proteins][interface2].insert(0,interface_dict[proteins][interface1][0])
    return(interface_dict)

def color_int_residues(pdb_path, residues_to_color, names) :
    """
    Color residues in interaction in a PDB file.
   
    Parameters:
    ----------
    pdb_path : string
    residues_to_color : dict
    names : string
    
    Returns:
    ----------
    """
    name_prot = names[0]
    save_line = str()
    chain1 = "B"
    with open(f'{pdb_path}/ranked_0.pdb', 'r') as file :
        for line in file:
            if line.startswith("ATOM") :
                chain2 = line[21]
                if chain1 != chain2 :
                   name_prot = names[1] #use new dict to color atoms
                res_num = line[22:26].strip()
                if res_num in residues_to_color[name_prot] : #change B-factor in color interaction residue
                    line = line[:60] + " 100  " + line[66:]
                else :
                    line = line[:60] + " 0    " + line[66:]
                chain1 = line[21]
            save_line += line
    with open(f'{pdb_path}/ranked_0.pdb', 'w') as writer:
        writer.write(save_line)

def plot_sequence_interface (file, cluster_dict) :
   """
   Generated figures for interface in one sequence.

   Parameters:
   ----------
   file : object of File_proteins class
   cluster_dict : dict

   Returns:
   """
   if not os.path.exists("./interface_fig/") :
      os.makedirs("./interface_fig/")
   sequence_dict = file.get_proteins_sequence()
   dict_inter = file.get_interface_dict()
   all_color= ['red','green', 'blue', 'orange', 'purple', 'cyan', 'magenta', 'yellow', 'pink', 'brown','lime', 'indigo', 'violet', 'turquoise', 'teal', 'crimson', 'gold', 'salmon', 'plum', 'chartreuse']
   for uniprotID_main in dict_inter.keys() :
      sequence = sequence_dict[uniprotID_main]
      indice_color = -1
      interface_done = dict()
      index_to_color = dict()
      uniprot_id_interface = dict()
      for interaction in dict_inter[uniprotID_main] : #list of residue + interface + UniprotID in interaction
         if interaction[0] not in interface_done.keys() : #if it's a new interface
            indice_color += 1
            interface_done[interaction[0]] = all_color[indice_color]
            uniprot_id_interface[interaction[len(interaction)-1]] = all_color[indice_color]
            for aa_to_color in interaction :
               if " " in aa_to_color :
                  if aa_to_color.split(" ")[1] not in index_to_color.keys() :
                     index_to_color[aa_to_color.split(" ")[1]] = [all_color[indice_color]]
                  if aa_to_color.split(" ")[1] in index_to_color.keys() and all_color[indice_color] not in index_to_color[aa_to_color.split(" ")[1]] : #add two colour if it's in two interface
                     index_to_color[aa_to_color.split(" ")[1]].append(all_color[indice_color])
               else : #for seconde residue table
                  if aa_to_color not in index_to_color.keys() :
                     index_to_color[aa_to_color] = [all_color[indice_color]]
                  if aa_to_color in index_to_color.keys() and all_color[indice_color] not in index_to_color[aa_to_color] : #add two colour if it's in two interface
                     index_to_color[aa_to_color].append(all_color[indice_color])
         else :
            uniprot_id_interface[interaction[len(interaction)-1]] = all_color[indice_color]
            for aa_to_color in interaction :
               if " " in aa_to_color :
                  if aa_to_color.split(" ")[1] not in index_to_color.keys() :
                     index_to_color[aa_to_color.split(" ")[1]] = [all_color[indice_color]]
                  if aa_to_color.split(" ")[1] in index_to_color.keys() and all_color[indice_color] not in index_to_color[aa_to_color.split(" ")[1]] : #add two colour if it's in two interface
                     index_to_color[aa_to_color.split(" ")[1]].append(all_color[indice_color])
               else : #for seconde residue table
                  if aa_to_color not in index_to_color.keys() :
                     index_to_color[aa_to_color] = [all_color[indice_color]]
                  if aa_to_color in index_to_color.keys() and all_color[indice_color] not in index_to_color[aa_to_color] : #add two colour if it's in two interface
                     index_to_color[aa_to_color].append(all_color[indice_color])
      line_adjust = 150 #max aa per line
      n_lines = (len(sequence) + line_adjust - 1) // line_adjust
      fig, ax = plt.subplots(figsize=(line_adjust / 4, n_lines*1.5)) #Adjust figsize
      for line_index in range(0, len(sequence), line_adjust) :
         sub_sequence = sequence[line_index:line_index + line_adjust]
         y_pos = -line_index // line_adjust * 1.5
         for i in range(len(sub_sequence)) :
            aa = sub_sequence[i]
            total_index = line_index + i
            if str(total_index + 1) in index_to_color.keys() :
               colors = index_to_color[str(total_index + 1)]
               height = 0.5 / len(colors)
               for color_index, color in enumerate(colors) :
                  ax.add_patch(plt.Rectangle((i, y_pos +color_index * height), 1, height, color=color))
               ax.text(i + 0.5, y_pos + 0.25, aa, ha='center', va='center', color='white')
            else :
               ax.add_patch(plt.Rectangle((i, y_pos), 1, 0.6, color="white"))
               ax.text(i + 0.5, y_pos + 0.25, aa, ha='center', va='center', color='black')
            if (total_index+1) % 10 == 0 or i == 0:
               ax.text(i + 0.5, y_pos + 0.5, str(total_index+1), ha='center', va='center', color='black', fontsize=7)
      for index_neigh, neigh in enumerate(uniprot_id_interface) :
         ax.text(index_neigh * 6, -n_lines * 2, neigh, ha='center', va='center', color=uniprot_id_interface[neigh], fontsize=8)
      ax.text(-2, 0.25, uniprotID_main, ha='right', va='center', color='black', fontsize=10, fontweight='bold')
      ax.set_xlim(0, line_adjust)
      ax.set_ylim(-n_lines*2, 1)  #Adjust high
      ax.axis('off')
      plt.savefig("./interface_fig/"+uniprotID_main+"_interface_fig.png", dpi=300, bbox_inches='tight')

def recover_prot_sequence(file, path_pkl) :
   list_proteins = file.get_proteins()
   new_dict_sequence = dict()
   for protein in list_proteins :
      with open(os.path.join(f'{path_pkl}/{protein}.pkl'), 'rb') as pkl_file :
         pickle_dict = pickle.load(pkl_file)
         new_dict_sequence[protein] = pickle_dict.sequence
   file.set_proteins_sequence(new_dict_sequence)


