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
            if int(SP_signal) > 0 :
                new_line2 = line2[int(SP_signal)+1:len(line2)]
                new_fasta_dict[save_key] = line2[int(SP_signal)+1:len(line2)]
                SP_signal = 0
            if line2[0] == ">" :
                if str(line2[1:len(line2)-1]) in prot_SP.keys() :
                    SP_signal = prot_SP[line2[1:len(line2)-1]]
                    save_key = line2[1:len(line2)-1]
            final_file = final_file + new_line2
    if len(new_fasta_dict) > 0 :
        file.find_prot_lenght(new_fasta_dict)
    cmd2 = "rm " + fasta_file
    os.system(cmd2)
    with open(fasta_file, "w") as new_file2 :
        new_file2.write(final_file)

def create_feature (file, data_dir, Path_Pickle_Feature) :
    """
    Launch command to generate features.

    Parameters:
    ----------
    file : object of class File_proteins
    data_dir : string
    Path_Pickle_Feature : string

    Returns:
    ----------
    """
    fasta_file = file.get_fasta_file()
    cmd = f"create_individual_features.py --fasta_paths=./{fasta_file} \--data_dir={data_dir} \--save_msa_files=True \--output_dir={Path_Pickle_Feature} \--max_template_date=2024-05-02 \--skip_existing=True \--use_mmseqs2=True"
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
    if len(new_proteins) > 0 : #if you have new proteins, generate MSA Depth
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
            plt.imshow(final,
                interpolation='nearest', aspect='auto',
                cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
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
                elif os.path.exists(f"./result_all_vs_all/{proteins[index_protein]}_and_{proteins[index2_protein]}/ranked_0.pdb") == False : #make interaction if doesn't exist and is not too long
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
    cmd1 = "export XLA_PYTHON_CLIENT_PREALLOCATE=false"
    cmd2 = "TF_FORCE_UNIFIED_MEMORY=true"
    cmd3 = "XLA_CLIENT_MEM_FRACTION=3.2"
    cmd4 =f"run_multimer_jobs.py --mode=custom \--num_cycle=3 \--num_predictions_per_model=1 \--compress_result_pickles=True \--output_path=./result_all_vs_all \--data_dir={data_dir} \--protein_lists=all_vs_all.txt \--monomer_objects_dir={Path_Pickle_Feature} \--noremove_keys_from_pickles"
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)

def add_indice_Q (dir_alpha) :
    """
    Generate indice_Q for all interactions.

    Parameters:
    ----------
    dir_alpha : string

    Returns:
    ----------
    """
    cmd4 = f"singularity exec --no-home --bind result_all_vs_all:/mnt {dir_alpha}/fold_analysis_latest.sif run_get_good_pae.sh --output_dir=/mnt --cutoff=10"
    os.system(cmd4)
    with open("result_all_vs_all/predictions_with_good_interpae.csv", "r") as file1 :
        reader = csv.DictReader(file1)
        all_lines = "jobs,pi_score,iptm_ptm,pDockQ,indice_Q\n"
        for row in reader :
            job = row['jobs']
            if '_and_' in job and row['pi_score'] != 'No interface detected' :
                indice_Q = ((float(row['pi_score'])+2.63)/5.26)*40+float(row['iptm_ptm'])*30+float(row['mpDockQ/pDockQ'])*30
                line =f'{row["jobs"]},{row["pi_score"]},{row["iptm_ptm"]},{row["mpDockQ/pDockQ"]},{str(indice_Q)}\n'
                all_lines = all_lines + line
    with open("result_all_vs_all/new_predictions_with_good_interpae.csv", "w") as file2 :
        file2.write(all_lines)

def create_out_fig (file) :
    """
    Generate result figure for validate interaction (indice_Q) and better interaction (indice_hQ).

    Parameters:
    ----------
    file : object of class File_proteins

    Returns:
    ----------
    """
    indice_Q_dict = file.get_indice_Q_dict()
    for interaction in indice_Q_dict.keys() :
        if float(indice_Q_dict[interaction]) >= 50 : #Plot figure of interest just for interesting interactions
            job1 = interaction[0] + "_and_" + interaction[1]
            plot_Distogram("./result_all_vs_all/" + job1) #need distogram key in pickle file
            make_table_res_int("./result_all_vs_all/" + job1)
    indice_hQ_dict = file.get_indice_hQ_dict()
    for homo_oligo in indice_hQ_dict.keys() :
        if float(indice_hQ_dict[homo_oligo][0]) >= 50 :
            job2 = homo_oligo #job2 = homo_oligo + "_homo_" + str(indice_hQ_dict[homo_oligo][1]) + "er" # wait AFPD homo release
            for count in range(1,indice_hQ_dict[homo_oligo][1]) :
                job2 += "_and_" + homo_oligo
            plot_Distogram("./result_homo_oligo/" + job2) #need distogram key in pickle file
            make_table_res_int("./result_homo_oligo/" + job2)

def make_table_res_int (path_int) :
        """
        Generate a table of residues in interactions.

        Parameters:
        ----------
        path_int : string

        Returns:
        ----------

        """
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', path_int + "/ranked_0.pdb")
        dict_int = dict()
        int_already_know = dict()
        proteins = path_int.split('/')[2].split('_and_')
        color_res = dict()
        color_res[names[0]] = set()
        color_res[names[1]] = set()
        atom_possible_contact = ["O","OH","NH2","NH1","OG","NE2","ND2","NZ","NE","N","OE1","OE2","OD2","OG1"] #hydrogen bond
        for model in structure:
            list_chain = model.get_list()
            for index1 in range(0,len(list_chain)) :
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
                                            if distance <= 3.64 and atom1.bfactor >=70 and atom2.bfactor >= 70 : #filtered on pLDDT and distance, be stringent to avoid false residue interaction (or maybe use PAE ?)
                                                res_int = chain1.get_id()+":"+residue1.get_resname()+" "+str(residue1.get_id()[1])," "+chain2.get_id()+":"+residue2.get_resname()+" "+str(residue2.get_id()[1])
                                                if chain1.get_id()+chain2.get_id() in dict_int.keys() : #to make different table for different interaction
                                                    if res_int in int_already_know.keys() and int_already_know[res_int] > str(distance) :
                                                        dict_int[chain1.get_id()+chain2.get_id()].remove([res_int[0]," "+res_int[1]," "+str(int_already_know[res_int])])
                                                        dict_int[chain1.get_id()+chain2.get_id()].append([res_int[0]," "+res_int[1]," "+str(distance)])
                                                        int_already_know[res_int] = str(distance)
                                                        color_res[names[0]].add(res_int[0])
                                                        color_res[names[1]].add(res_int[1])
                                                    elif res_int in int_already_know.keys() and int_already_know[res_int] < str(distance) : #skip double interaction with differents atoms
                                                        pass
                                                    else :
                                                        dict_int[chain1.get_id()+chain2.get_id()].append([res_int[0]," "+res_int[1]," "+str(distance)])
                                                        int_already_know[res_int] = str(distance)
                                                        color_res[names[0]].add(res_int[0])
                                                        color_res[names[1]].add(res_int[1])
                                                else :
                                                    dict_int[chain1.get_id()+chain2.get_id()] = [["Chain "+chain1.get_id()," Chain "+chain2.get_id()," Distance Ä"]]
                                                    dict_int[chain1.get_id()+chain2.get_id()].append([res_int[0]," "+res_int[1]," "+str(distance)])
                                                    int_already_know[res_int] = str(distance)
                                                    color_res[names[0]].add(res_int[0])
                                                    color_res[names[1]].add(res_int[1])
                                            else :
                                                pass
        for chains in dict_int.keys() :
            fileout = chains+"_res_int.csv"
            np_table = np.array(dict_int[chains])
            with open(f"{path_int}/"+fileout, "w", newline="") as file :
                mywriter = csv.writer(file, delimiter=",")
                mywriter.writerows(np_table)
            print("Write table")
        color_int_residues(path_int,color_res,names) #color residue in interaction on the pdb
    
#def make_table_res_int (file, path_int) : #need key distogram in pickle file
#    """
#    Generate a table of residue in interactions.
#
#    Parameters:
#    ----------
#    file : object of class File_proteins
#    path_int : string
#
#    Returns:
#    ----------
#    """
#    ranking_results = json.load(open(os.path.join(f'{path_int}/ranking_debug.json')))
#    best_model = ranking_results["order"][0]
#    lenght_prot = file.get_lenght_prot()
#    seq_prot = file.get_proteins_sequence()
#    names = path_int.split("/")[2].split("_and_")
#    chains = path_int.split("/")[2]
#    dict_interface = dict()
#    color_res = dict()
#    color_res[names[0]] = set()
#    color_res[names[1]] = set()
#    with open(os.path.join(f'{path_int}/result_{best_model}.pkl.gz'), 'rb') as gz_file :
#        pickle_dict = pickle.load(gzip.open(gz_file))
#        pae_mtx = pickle_dict['predicted_aligned_error']#take PAE
#        bin_edges = pickle_dict["distogram"]["bin_edges"]#take distogram for distance
#        bin_edges = np.insert(bin_edges, 0, 0)
#        distogram_softmax = softmax(pickle_dict["distogram"]["logits"], axis=2)
#        dist = np.sum(np.multiply(distogram_softmax, bin_edges), axis=2) #center of the residue
#        dict_interface[chains] = [[names[0]," "+names[1]," Distance Ä"," PAE score"]]
#        for line in range(lenght_prot[names[1]],len(dist)) :
#            hori_index = 0
#            for distance in dist[line] :
#                hori_index += 1
#                if hori_index < lenght_prot[names[1]] :
#                    if distance <= 10 :
#                        if pae_mtx[line][hori_index] < 5 :
#                            residue1 = seq_prot[names[0]][line-lenght_prot[names[1]]]
#                            residue2 = seq_prot[names[1]][hori_index]
#                            dict_interface[chains].append([residue1+" "+str(line-lenght_prot[names[1]]+1)," "+residue2+" "+str(hori_index+1)," "+str(distance), " "+str(pae_mtx[line][hori_index])]) #+1 to match with pdb model
#                            color_res[names[0]].add(line-lenght_prot[names[1]]+1)
#                            color_res[names[1]].add(hori_index+1)
#    file.define_interface(dict_interface[chains],names) #update interaction interface
#    color_int_residues(path_int,color_res,names) #color residue in interaction on the pdb
#    fileout = chains+"_res_int.csv"
#    np_table = np.array(dict_interface[chains])
#    with open(f"{path_int}/"+fileout, "w", newline="") as file :
#        mywriter = csv.writer(file, delimiter=",")
#        mywriter.writerows(np_table)

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
    cmd1 = "export XLA_PYTHON_CLIENT_PREALLOCATE=false"
    cmd2 = "TF_FORCE_UNIFIED_MEMORY=true"
    cmd3 = "XLA_CLIENT_MEM_FRACTION=3.2"
    cmd4 =f"run_multimer_jobs.py --mode=custom \--num_cycle=3 \--num_predictions_per_model=1 \--compress_result_pickles=True \--output_path=./result_homo_oligo \--data_dir={data_dir} \--protein_lists=homo_oligo.txt \--monomer_objects_dir={Path_Pickle_Feature} \--noremove_keys_from_pickles"
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)

def add_indice_hQ (dir_alpha) :
    """
    Generate a score table for all homo-oligomer.

    Parameters:
    ----------
    dir_alpha : string
        
    Returns:
    ----------
    """
    cmd4 = f"singularity exec --no-home --bind result_homo_oligo:/mnt {dir_alpha}/fold_analysis_latest.sif run_get_good_pae.sh --output_dir=/mnt --cutoff=10"
    os.system(cmd4)
    with open("./result_homo_oligo/predictions_with_good_interpae.csv", "r") as file1 :
        reader = csv.DictReader(file1)
        all_lines = "jobs,pi_score,iptm_ptm,indice_hQ\n"
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
        number_oligo = len(row["jobs"].split("_"))
        if len(save_pi_score[key]) > int(number_oligo) : #if model have more interface than number of homo-oligomerization
            new_sum_pi_score = 0
            save_pi_score[key].sort(reverse=True)
            for index in range(0,int(number_oligo)) :
                new_sum_pi_score += save_pi_score[key][index]
                indice_hQ = (((float(new_sum_pi_score)/int(number_oligo))+2.63)/5.26)*60+float(row['iptm_ptm'])*40 #cause iptm_ptm are always same for each homo of same protein
            line =f'{key},{str(all_homo[key][0])},{row["iptm_ptm"]},{str(indice_hQ)}\n'
            all_lines += line
        else :
            indice_hQ = (((float(all_homo[key][0])/all_homo[key][1])+2.63)/5.26)*60+float(row['iptm_ptm'])*40 #cause iptm_ptm is always same for each homo of same protein
            line =f'{key},{str(all_homo[key][0])},{row["iptm_ptm"]},{str(indice_hQ)}\n'
            all_lines += line
    with open("./result_homo_oligo/new_predictions_with_good_interpae.csv", "w") as file2 :
        file2.write(all_lines)

def generate_interaction_network (file) :
    """
    Generate interaction network, colored by indice_Q.
   
    Parameters:
    ----------
    file : object of File_proteins class

    Returns:
    ----------
    """
    valid_interactions = list()
    indice_Q_dict = file.get_indice_Q_dict()
    for interactions in indice_Q_dict.keys() :
        names = [interactions[0], interactions[1]]
        if names not in [x[0] for x in valid_interactions] and float(indice_Q_dict[interactions]) >= 50 :
            valid_interactions.append([names, float(indice_Q_dict[interactions])])
    indice_hQ_dict = file.get_indice_hQ_dict()
    for homo_oligomer in indice_hQ_dict.keys() :
        if float(indice_hQ_dict[homo_oligomer][0]) >= 50 :
            valid_interactions.append([[homo_oligomer,homo_oligomer], indice_hQ_dict[homo_oligomer][1]])
    int_graph = nx.Graph()
    list_inter_score = list()
    prots = set()
    dict_name = file.get_names()
    with open('table.cyt', 'w') as f :
        f.write(('source,targer,interaction,score\n'))
        for inter, score in valid_interactions :
            inter0 = inter[0]+f"({dict_name[inter[0]]})" #set uniprotID with the name of protein
            inter1 = inter[1]+f"({dict_name[inter[1]]})"
            f.write(f'{inter0},{inter1},pp,{round(score,2)}\n')
            prots.add(inter0)
            prots.add(inter1)
            list_inter_score.append((inter0,inter1,round(score,2)))
    prots = list(prots)
    int_graph.add_nodes_from(prots)
    int_graph.add_weighted_edges_from(list_inter_score)
    fig, ax = plt.subplots(figsize=(10, 6))
    pos = nx.spring_layout(int_graph, k = len(prots)+50, scale = 3, seed=8)
    nx.draw_networkx_nodes(int_graph,pos)
    nx.draw_networkx_edges(int_graph,pos, node_size=100, edgelist=int_graph.edges, width=1, style="dashed")
    nx.draw_networkx_labels(int_graph, pos, font_size=10, font_family="sans-serif", hide_ticks = 'True')
    edge_labels = nx.get_edge_attributes(int_graph, "weight")
    homo_label_dict = dict()
    color_label = dict()
    for prot_int in edge_labels.keys() :
        if len(str(edge_labels[prot_int])) < 3 :
            homo_label_dict[prot_int] = edge_labels[prot_int]
        else :
            color_label[prot_int] = edge_labels[prot_int]
    selected_weights = {edge: edge_labels[edge] for edge in color_label.keys()}
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
    plt.colorbar(sm, ax=ax, label="indice_Q")
    plt.savefig("network.png")
    plt.close()

def generate_heatmap (file) :
    """
    Generate a heatmap of interaction scores between all proteins.
   
    Parameters:
    ----------
    file : object of File_proteins class

    Returns:
    ----------
    """
    indice_Q_data_matrix = list()
    iptm_ptm_data_matrix = list()
    indice_Q_dict = file.get_indice_Q_dict()
    proteins_list = file.get_proteins()
    proteins_name = file.get_names()
    index_prot = list()
    iptm_ptm_dict = dict()
    with open("result_all_vs_all/new_predictions_with_good_interpae.csv", "r") as file1 :
        reader = csv.DictReader(file1)
        for row in reader :
            iptm_ptm_dict[row['jobs']] = row["iptm_ptm"]
    for protein in proteins_list :
        index_prot.append(protein+"_"+proteins_name[protein])
    for protein1 in proteins_list :
        indice_Q_line = list()
        iptm_ptm_line = list()
        for protein2 in proteins_list :
            if protein1 == protein2 :
                indice_Q_line.append(0)
                iptm_ptm_line.append(0)
            else :
                if (protein1,protein2) in indice_Q_dict.keys() :
                    indice_Q_line.append(float(indice_Q_dict[(protein1,protein2)]))
                    iptm_ptm_line.append(float(iptm_ptm_dict[protein1+"_and_"+protein2]))
                elif (protein2,protein1) in indice_Q_dict.keys() :
                    indice_Q_line.append(float(indice_Q_dict[(protein2,protein1)]))
                    iptm_ptm_line.append(float(iptm_ptm_dict[protein2+"_and_"+protein1]))
                else :
                    indice_Q_line.append(0)
                    iptm_ptm_line.append(0)
                    print (protein1 + " and " + protein2 + " are not in score table or have bad inter PAE")
        indice_Q_data_matrix.append(indice_Q_line)
        iptm_ptm_data_matrix.append(iptm_ptm_line)
    indice_Q_complet_matrix = pd.DataFrame(indice_Q_data_matrix,index = index_prot, columns = index_prot)
    iptm_ptm_complet_matrix = pd.DataFrame(iptm_ptm_data_matrix,index = index_prot, columns = index_prot)
    ax2 = seaborn.heatmap(indice_Q_complet_matrix, cbar_kws = {'label' : 'indice_Q'}, xticklabels=True, yticklabels=True)
    ax2.figure.tight_layout()
    plt.savefig("indice_Q_heatmap.png")
    plt.close()
    ax2 = seaborn.heatmap(iptm_ptm_complet_matrix, cbar_kws = {'label' : 'iptm_ptm'}, xticklabels=True, yticklabels=True)
    ax2.figure.tight_layout()
    plt.savefig("iptm_ptm_heatmap.png")
    plt.close()

###Probably add to generate_interaction_network ???
def redef_interface (file) :
    """
    Compare interface in function of smaller interface and classify them with a letter representing the interface.
   
    Parameters:
    ----------
    file : object of File_proteins class

    Returns:
    ----------
    """
    alphabet = string.ascii_lowercase
    interface_dict = file.get_interface_dict()
    for proteins in interface_dict.keys() :
        already_inter = list()
        interface_dict[proteins] = sorted(interface_dict[proteins]) #sorted all interface in function of number of resiudes
        for interface1 in range(len(interface_dict[proteins])) :
            if interface1 == 0 : #if it's the first interface, define a
                interface_dict[proteins][interface1].insert(0,alphabet[0])
                already_inter.append(alphabet[0])
            for interface2 in range(interface1+1,len(interface_dict[proteins])) :
                list_inter = list(set(interface_dict[proteins][interface1]).intersection(set(interface_dict[proteins][interface2])))
                simi_inter = len(list_inter)/(len(set(interface_dict[proteins][interface1]).union(set(interface_dict[proteins][interface2])))-3) #indice jaccard # -3 just to remove interface 'a' and uniprotID from .union()
                print(simi_inter, interface1, interface2)
                if simi_inter < 0.50 : #create a new interface
                    if interface_dict[proteins][interface2][0] in already_inter : #Don't create new interface if it already has one
                        pass
                    else :
                        interface_dict[proteins][interface2].insert(0,alphabet[interface2])
                        already_inter.append(alphabet[interface2])
                        print(already_inter)
                else : #if interfaces got more than 0.50 of same residues, it's the same interface
                    interface_dict[proteins][interface2].insert(0,interface_dict[proteins][interface1][0])
    print(interface_dict)
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
    name_prot = names[1]
    save_line = str()
    chain1 = "B"
    residues_to_color
    with open(f'{pdb_path}/ranked_0.pdb', 'r') as file :
        for line in file:
            if line.startswith("ATOM") :
                chain2 = line[21]
                if chain1 == chain2 :
                    res_num = int(line[22:26].strip())
                    if res_num in residues_to_color[name_prot] : #change B-factor in color interaction residue
                        line = line[:60] + " 100  " + line[66:]
                    else :
                        line = line[:60] + " 0    " + line[66:]
                else :
                    name_prot = names[0] #use new dict to color atoms
                    res_num = int(line[22:26].strip())
                    if res_num in residues_to_color[name_prot] :
                        line = line[:60] + " 100  " + line[66:]
                    else :
                        line = line[:60] + " 0    " + line[66:]
                chain1 = line[21]
            save_line += line
    with open(f'{pdb_path}/ranked_0.pdb', 'w') as writer:
        writer.write(save_line)

