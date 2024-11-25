import numpy as np
import pickle
import csv
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from Bio import PDB
import copy
import glob
import logging
from scipy.special import softmax
import os
import networkx as nx
from math import *
import seaborn
import pandas as pd
from adjustText import adjust_text

from File_proteins import *

def define_path() :
    """
    Take all path from conf.txt, and set up on dictionnary.
    
    Parameters:
    ----------

    Returns:
    ----------
    path_dict : dict
    """
    path_dict = dict()
    with open("conf.txt", "r") as file :
        while True:
            lines = file.readline()
            if not lines:
                break
            else :
                if lines[0:17] == "Path_Uniprot_ID :" :
                    if len(lines[17:len(lines)]) < 4 :
                        print("Path to Uniprot file is empty")
                        exit()
                    else :
                        path = lines[17:len(lines)].strip().strip('\n')
                        path_dict["Path_Uniprot_ID"] = path

                elif lines[0:21] == "Path_AlphaFold_Data :" :
                    if len(lines[21:len(lines)]) < 4 :
                        print("Path to AlphaFold data is empty, set by default on ./alphadata")
                        path_dict["Path_AlphaFold_data"] = "./alphadata"
                    else :
                        path = lines[21:len(lines)].strip().strip('\n')
                        if path[len(path)-1] == "/" :
                            path = path[0:len(path)-1]
                        else :
                            pass
                        path_dict["Path_AlphaFold_data"] = path

                elif lines[0:24] == "Path_Singularity_Image :" :
                    if len(lines[24:len(lines)]) < 4 :
                        print("Path to Singularity image is empty")
                        exit()
                    else :
                        path = lines[24:len(lines)].strip().strip('\n')
                        path_dict["Path_Singularity_Image"] = path

                elif lines[0:21] == "Path_Pickle_Feature :" :
                    if len(lines[21:len(lines)]) < 4 :
                        print("Path to Feature data is empty, set by default on ./feature")
                        path_dict["Path_Pickle_Feature"] = "./feature"
                    else :
                        path = lines[21:len(lines)].strip().strip('\n')
                        if path[len(path)-1] == "/" :
                            path = path[0:len(path)-1]
                        else :
                            pass
                        path_dict["Path_Pickle_Feature"] = path
    return(path_dict)

def remove_SP (file, org) :
        """
        Creating a new fasta file without signal peptide.

        Parameters:
        ----------
        file : object of class File_proteins
        org : type of organism
        
        Returns:
        ----------

        """
        final_file = str()
        SP_signal = 0
        prot_SP = dict()
        fasta_file = file.get_fasta_file()
        cmd = "signalp -fasta " + fasta_file + " -org " + org
        os.system(cmd)
        file_signalp = fasta_file.replace(".fasta","_summary.signalp5")
        with open(file_signalp,"r") as fh :
            for line in fh :
                new_line = line.split("\t")
                if new_line[1] != "OTHER" and new_line[0][0] != "#" :
                    prot_SP[new_line[0]] = new_line[6][11:13]
        with open(fasta_file, "r") as fa_file :
           for line2 in fa_file :
              new_line2 = line2
              if int(SP_signal) > 0 :
                 new_line2 = line2[int(SP_signal)+1:len(line2)]
                 SP_signal = 0
              if line2[0] == ">" :
                 if str(line2[1:len(line2)-1]) in prot_SP.keys() :
                    SP_signal = prot_SP[line2[1:len(line2)-1]]
              final_file = final_file + new_line2
        cmd2 = "rm " + fasta_file
        os.system(cmd2)
        with open(fasta_file, "w") as new_file2 :
           new_file2.write(final_file)

def create_feature (file, env_feature, data_dir, Path_Pickle_Feature) :
        """
        Launch command to generate features.

        Parameters:
        ----------
        env_feature : string
        data_dir : string
        file : object of class File_proteins

        Returns:
        ----------

        """
        fasta_file = file.get_fasta_file()
        cmd = f"#!/bin/bash --login \n source ~/.bashrc \n conda activate {env_feature}\n create_individual_features.py --fasta_paths=./{fasta_file} \--data_dir={data_dir} \--save_msa_files=True \--output_dir={Path_Pickle_Feature} \--max_template_date=2024-05-02 \--skip_existing=False"
        cmd2 = f"create_individual_features.py --fasta_paths=./{fasta_file} \--data_dir={data_dir} \--save_msa_files=True \--output_dir={Path_Pickle_Feature} \--max_template_date=2024-05-02 \--skip_existing=False"
        cmd3 = "#!/bin/bash --login \n source ~/.bashrc \n conda deactivate"
        if env_feature != None :
            os.system(cmd)
            os.system(cmd3)
        else :
            os.system(cmd2)

def Make_all_MSA_coverage (file,Path_Pickle_Feature) :
        """
        Creating a new fasta file without signal peptide.

        Parameters:
        ----------
        file : object of class File_proteins

        Returns:
        ----------

        """
        bad_MSA = str()
        proteins = file.get_proteins()
        for prot in proteins :
            pre_feature_dict = pickle.load(open(f'{Path_Pickle_Feature}/{prot}.pkl','rb'))
            feature_dict = pre_feature_dict.feature_dict
            msa = feature_dict['msa']
            if len(msa) <= 100 :
               bad_MSA += prot + " : " + str(len(msa)) + " sequences\n"
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
        with open("Bad_MSA.txt", "w") as MSA_file :
            MSA_file.write(bad_MSA)

def generate_APD_script (file, max_aa) :
        """
        Write two scripts in local to use AlphaPulldown, this scripts are build in function of maximum amino acid.

        Parameters:
        ----------
        max_aa : integer
        file : object of class File_proteins

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
                if int_lenght <= max_aa :
                    all_vs_all_script = all_vs_all_script + proteins[index_protein] + ";" +  proteins[index2_protein]+ "\n"
                else :
                    OOM_int = OOM_int + proteins[index_protein] + ";" +  proteins[index2_protein]+ "\n"
            lenght_homo = lenght
            for nbr_homo in range(2,21) :
                lenght_homo += lenght
                if lenght_homo <= max_aa :
                    homo_oligo_script = homo_oligo_script + proteins[index_protein] + "," + str(nbr_homo) + "\n"
                else :
                    OOM_int = OOM_int + proteins[index_protein] + "," + str(nbr_homo) + "\n"
        with open("homo_oligo.txt", "w") as homo_file:
            homo_file.write(homo_oligo_script)
        with open("all_vs_all.txt", "w") as all_file:
            all_file.write(all_vs_all_script)
        with open("OOM_int.txt", "w") as OOM_file :
            OOM_file.write(OOM_int)

### Generating Multimers

def Make_all_vs_all (env_multimers, data_dir, Path_Pickle_Feature) :
        """
        Use Alphapulldown script to generate all versus all interactions.

        Parameters:
        ----------
        env_multimers : string
        data_dir : string

        Returns:
        ----------

        """
        cmd = f"#!/bin/bash --login \n source ~/.bashrc \n conda activate {env_multimers}\n run_multimer_jobs.py --mode=custom \--num_cycle=3 \--num_predictions_per_model=1 \--compress_result_pickles=True \--output_path=result_all_vs_all \--data_dir={data_dir} \--protein_lists=all_vs_all.txt \--monomer_objects_dir={Path_Pickle_Feature}"
        cmd2 =f"run_multimer_jobs.py --mode=all_vs_all \--num_cycle=3 \--num_predictions_per_model=1 \--compress_result_pickles=True \--output_path=./result_all_vs_all \--data_dir={data_dir} \--protein_lists=all_vs_all.txt \--monomer_objects_dir={Path_Pickle_Feature}"
        cmd3 = "#!/bin/bash --login \n source ~/.bashrc \n conda deactivate"
        if env_multimers != None :
            os.system(cmd)
            os.system(cmd3)
        else :
            os.system(cmd2)

def add_iQ_score (dir_alpha) :
        """
        Launch command to generate all_vs_all result.

        Parameters:
        ----------
        env_multimers : string
        data_dir : string

        Returns:
        ----------

        """
        cmd4 = f"singularity exec --no-home --bind result_all_vs_all:/mnt {dir_alpha}/alpha-analysis_jax_0.4.sif run_get_good_pae.sh --output_dir=/mnt --cutoff=10"
        os.system(cmd4)
        with open("result_all_vs_all/predictions_with_good_interpae.csv", "r") as file1 :
            reader = csv.DictReader(file1)
            all_lines = "jobs,interface,Num_intf_residues,Polar,Hydrophobhic,Charged,contact_pairs, sc, hb, sb, int_solv_en, int_area,pi_score,iptm_ptm,iptm,mpDockQ/pDockQ,iQ_score\n"
            for row in reader :
                job = row['jobs']
                if '_and_' in job and row['pi_score'] != 'No interface detected' :
                    iQ_score = ((float(row['pi_score'])+2.63)/5.26)*40+float(row['iptm_ptm'])*30+float(row['mpDockQ/pDockQ'])*30
                    line = row['jobs']+","+row['interface']+","+row['Num_intf_residues']+","+row['Polar']+","+row['Hydrophobhic']+","+row['Charged']+","+row['contact_pairs']+","+row[' sc']+","+row[' hb']+","+row[' sb']+","+row[' int_solv_en']+","+row[' int_area']+","+row['pi_score']+","+row['iptm_ptm']+","+row['iptm']+","+row['mpDockQ/pDockQ']+","+str(iQ_score)+"\n"
                    all_lines = all_lines + line
        with open("result_all_vs_all/predictions_with_good_interpae.csv", "w") as file2 :
            file2.write(all_lines)

def create_out_fig (file) :
        """
        Generate result figure for validate interaction.

        Parameters:
        ----------

        Returns:
        ----------

        """
        iQ_score_dict = file.get_iQ_score_dict()
        print(iQ_score_dict)
        for interaction in iQ_score_dict.keys() :
            if float(iQ_score_dict[interaction]) >= 35 : #Plot figure of interest just for interesting interactions
                job1 = interaction[0] + "_and_" + interaction[1]
                #plot_Distogram("./result_all_vs_all/" + job1)
                make_table_res_int("./result_all_vs_all/" + job1)
        hiQ_score_dict = file.get_hiQ_score_dict()
        #for homo_oligo in hiQ_score_dict.keys() :
         #   if float(hiQ_score_dict[homo_oligo][0]) >= 50 :
          #      job2 = homo_oligo + "_homo_" + hiQ_score_dict[homo_oligo][1] + "er"
           #     plot_Distogram("./result_homo_oligo/" + job2)
            #    make_table_res_int("./result_homo_oligo/" + job2)

def make_table_res_int (path_int) :
        """
        Generate a table of residue in interactions.

        Parameters:
        ----------
        path_int : string

        Returns:
        ----------

        """
        def make_table_res_int (file, path_int) :
        """
        Generate a table of residue in interactions with distance and PAE.

        Parameters:
        ----------
        file : object of class File_proteins
        path_int : string

        Returns:
        ----------

        """
        ranking_results = json.load(open(os.path.join(f'{path_int}/ranking_debug.json')))
        best_model = ranking_results["order"][0]
        lenght_prot = file.get_lenght_prot()
        seq_prot = file.get_proteins_sequence()
        names = path_int.split("/")[2].split("_and_")
        chains = path_int.split("/")[2]
        dict_interface = dict()
        with open(os.path.join(f'{path_int}/result_{best_model}.pkl.gz'), 'rb') as gz_file :
            pickle_dict = pickle.load(gzip.open(gz_file))
            pae_mtx = pickle_dict['predicted_aligned_error']
            bin_edges = pickle_dict["distogram"]["bin_edges"]
            bin_edges = np.insert(bin_edges, 0, 0)
            distogram_softmax = softmax(pickle_dict["distogram"]["logits"], axis=2)
            dist = np.sum(np.multiply(distogram_softmax, bin_edges), axis=2) #center of the residue
            dict_interface[chains] = [[names[0]," "+names[1]," Distance Ä"," PAE score"]]
            for line in range(lenght_prot[names[1]],len(dist)) :
                hori_index = 0
                for distance in dist[line] :
                    hori_index += 1
                    if hori_index <= lenght_prot[names[1]] :
                        if distance <= 7 :
                            if pae_mtx[line][hori_index] < 10 :
                                residue1 = seq_prot[names[0]][line-lenght_prot[names[1]]]
                                residue2 = seq_prot[names[1]][hori_index]
                                dict_interface[chains].append([residue1+" "+str(line-lenght_prot[names[1]]+1)," "+residue2+" "+str(hori_index+1)," "+str(distance), " "+str(pae_mtx[line][hori_index])]) #+1 to match with pdb model
        file.define_interface(dict_interface[chains],names)
        fileout = chains+"_res_int.csv"
        np_table = np.array(dict_interface[chains])
        with open(f"{path_int}/"+fileout, "w", newline="") as file :
            mywriter = csv.writer(file, delimiter=",")
            mywriter.writerows(np_table)

def plot_Distogram (job) :
        """
        Generate distogram for interactions of interest.

        Parameters:
        ----------
        job : string
        
        Returns:
        ----------

        """
        pickle_list = glob.glob(job + "/result_*.pkl")
        for i, pickle_output in enumerate(pickle_list):
            logging.warning(
                f"Processing pickle file {i+1}/{len(pickle_list)}: {pickle_output}")
            with open(pickle_output, "rb") as p:
                results = pickle.load(p)
                bin_edges = results["distogram"]["bin_edges"]
                bin_edges = np.insert(bin_edges, 0, 0)
                distogram_softmax = softmax(results["distogram"]["logits"], axis=2)
                dist = np.sum(np.multiply(distogram_softmax, bin_edges), axis=2)
                np.savetxt(f"{pickle_output}.dmap", dist)
                lenght_list = []
            with open(pickle_output,'rb') as handle :
                result = pickle.load(handle)
                for seq in result["seqs"] :
                    lenght_list.append(len(seq))
            print("make png")
            initial_lenght = 0
            fig, ax = plt.subplots()
            d = ax.imshow(dist)
            plt.colorbar(d, ax=ax, fraction=0.046, pad=0.04)
            ax.title.set_text("Distance map")
            for index in range(len(lenght_list)-1) :
                initial_lenght += lenght_list[index]
                ax.axhline(initial_lenght, color="black", linewidth=1.5)
                ax.axvline(initial_lenght, color="black", linewidth=1.5)
            plt.savefig(f"{pickle_output}.dmap.png", dpi=600)
            plt.close()

def Make_homo_oligo (env_multimers, data_dir, Path_Pickle_Feature) :
        """
        Use Alphapulldown script to generate all homo oligomer.

        Parameters:
        ----------
        env_multimers : string
        data_dir : string
        
        Returns:
        ----------

        """
        cmd = f"#!/bin/bash --login \n source ~/.bashrc \n conda activate {env_multimers}\n run_multimer_jobs.py --mode=homo-oligomer \--output_path=result_homo_oligo \--num_cycle=3 \--compress_result_pickles=True \--oligomer_state_file=homo_oligo.txt \--monomer_objects_dir={Path_Pickle_Feature} \--data_dir={data_dir} \--remove_result_pickles=False"
        cmd2 =f"run_multimer_jobs.py --mode=homo-oligomer \--output_path=result_homo_oligo \--num_cycle=3 \--compress_result_pickles=True \--oligomer_state_file=homo_oligo.txt \--monomer_objects_dir={Path_Pickle_Feature} \--data_dir={data_dir} \--remove_result_pickles=False"
        cmd3 = "#!/bin/bash --login \n source ~/.bashrc \n conda deactivate"
        if env_multimers != None :
            os.system(cmd)
            os.system(cmd3)
        else :
            os.system(cmd2)

def add_hiQ_score (dir_alpha) :
        """
        Generate a score table for all homo oligomer.

        Parameters:
        ----------
        dir_alpha : string
        
        Returns:
        ----------

        """
        #cmd4 = f"singularity exec --no-home --bind result_homo_oligo:/mnt {dir_alpha}/alpha-analysis_jax_0.4.sif run_get_good_pae.sh --output_dir=/mnt --cutoff=10"
        #os.system(cmd4)
        with open("result_homo_oligo/predictions_with_good_interpae.csv", "r") as file1 :
            reader = csv.DictReader(file1)
            all_lines = "jobs,pi_score,iptm_ptm,hiQ_score\n"
            all_homo = dict()
            for row in reader :
                job = row['jobs']
                if 'homo' in job and row['pi_score'] != 'No interface detected' :
                    if job not in all_homo.keys() :
                        all_homo[job] = (row['pi_score'],1,row)
                    else :
                        sum_pi_score = float(all_homo[job][0]) + float(row['pi_score'])
                        sum_int = all_homo[job][1] + 1
                        all_homo[job] = (sum_pi_score,sum_int,row)
            for key in all_homo.keys() :
                row = all_homo[key][2]
                hiQ_score = (((float(all_homo[key][0])/all_homo[key][1])+2.63)/5.26)*60+float(row['iptm_ptm'])*40 #cause iptm_ptm is always same for each homo of same protein
                line = key+","+str(all_homo[key][0])+","+row['iptm_ptm']+","+str(hiQ_score)+"\n"
                all_lines = all_lines + line
        with open("result_homo_oligo/predictions_with_good_interpae.csv", "w") as file2 :
            file2.write(all_lines)
            

def generate_interaction_network (file) :
    """
    Generate interaction network.
   
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
        if names not in [x[0] for x in valid_interactions] and float(iQ_score_dict[interactions]) >= 35 :
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
            inter0 = inter[0]+f"({dict_name[inter[0]]})" #set uniprotID with the name of protein
            inter1 = inter[1]+f"({dict_name[inter[1]]})"
            f.write(f'{inter0},{inter1},pp,{round(score,2)}\n')
            prots.add(inter0)
            prots.add(inter1)
            list_inter_score.append((inter0,inter1,round(score,2)))
    prots = list(prots)
    int_graph.add_nodes_from(prots)
    int_graph.add_weighted_edges_from(list_inter_score)
    fig, ax = plt.subplots(figsize=(8, 6))
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
    plt.savefig("network.png")
    plt.close()


#def generate_interaction_network2 (file) :
 #   interactions = list()
  #  with open("result_all_vs_all/predictions_with_good_interpae.csv", "r") as file1 :
   #     reader1 = csv.DictReader(file1)
    #    for row in reader1 :
     #       names = row['jobs'].split('_and_')
      #      if names not in [x[0] for x in interactions] and float(row['iQ_score']) >= 35 :
       #         interactions.append([names, float(row['iQ_score'])])
#    best_homo = dict()
 #   with open("result_homo_oligo/predictions_with_good_interpae.csv", "r") as file2 :
  #      reader2 = csv.DictReader(file2)
   #     for row in reader2 :
    #        if float(row["hiQ_score"]) >= 50 :
     #           prot_name = row['jobs'].split("_homo_")[0]
      #          if prot_name not in best_homo.keys() or float(row['hiQ_score']) >= best_homo[prot_name][0] :
       #             number_homo = int((row['jobs'].split("homo_")[1]).split("er")[0]) #to take the number of homo-oligomerisation of the protein and this score
        #            best_homo[prot_name] = (float(row['hiQ_score']),number_homo) #to take the number of homo-oligomerisation of the protein and this score
#    for key in best_homo :
 #       interactions.append([[key,key], best_homo[key][1]])
  #  H = xgi.Hypergraph()
   # list_inter = list()
    #prots = set()
#    dict_name = file.get_names()
 #   for inter, score in interactions :
  #      inter0 = inter[0]#+f"({dict_name[inter[0]]})" #set uniprotID with the name of protein
   #     inter1 = inter[1]#+f"({dict_name[inter[1]]})"
    #    prots.add(inter0)
     #   prots.add(inter1)
      #  list_inter.append([inter0,inter1])
#    prots = list(prots)
 #   H.add_nodes_from(prots)
  #  H.add_edges_from(list_inter)
   # pos = xgi.barycenter_spring_layout(H, seed=1)
    #print(H.edge_id())
#    xgi.draw(H, node_labels=True,hyperedge_labels=True, pos=pos)
 #   plt.savefig("network.png")
  #  plt.close()

def generate_heatmap (file):
    """
    Generate heatmap of interactions scores.
   
    Parameters:
    ----------
    file : object of File_proteins class

    Returns:
    ----------
    """
    iQ_data_matrix = list()
    iptm_ptm_data_matrix = list()
    iQ_score_dict = file.get_iQ_score_dict()
    proteins_list = file.get_proteins()
    proteins_name = file.get_names()
    index_prot = list()
    iptm_ptm_dict = dict()
    with open("result_all_vs_all/predictions_with_good_interpae.csv", "r") as file1 :
        reader = csv.DictReader(file1)
        for row in reader :
            iptm_ptm_dict[row['jobs']] = row["iptm_ptm"]
    for protein in proteins_list :
        index_prot.append(protein+"_"+proteins_name[protein])
    for protein1 in proteins_list :
        iQ_line = list()
        iptm_ptm_line = list()
        for protein2 in proteins_list :
            if protein1 == protein2 :
                iQ_line.append(0)
                iptm_ptm_line.append(0)
            else :
                if (protein1,protein2) in iQ_score_dict.keys() :
                    iQ_line.append(float(iQ_score_dict[(protein1,protein2)]))
                    iptm_ptm_line.append(float(iptm_ptm_dict[protein1+"_and_"+protein2]))
                elif (protein2,protein1) in iQ_score_dict.keys() :
                    iQ_line.append(float(iQ_score_dict[(protein2,protein1)]))
                    iptm_ptm_line.append(float(iptm_ptm_dict[protein2+"_and_"+protein1]))
                else :
                    iQ_line.append(0)
                    iptm_ptm_line.append(0)
                    print (protein1 + " and " + protein2 + " are not in score table or have bad PAE")
        iQ_data_matrix.append(iQ_line)
        iptm_ptm_data_matrix.append(iptm_ptm_line)
    iQ_complet_matrix = pd.DataFrame(iQ_data_matrix,index = index_prot, columns = index_prot)
    iptm_ptm_complet_matrix = pd.DataFrame(iptm_ptm_data_matrix,index = index_prot, columns = index_prot)
    ax2 = seaborn.heatmap(iQ_complet_matrix, cbar_kws = {'label' : 'iQ_score'}, xticklabels=True, yticklabels=True)
    ax2.figure.tight_layout()
    plt.savefig("iQ_score_heatmap.png")
    plt.close()
    ax2 = seaborn.heatmap(iptm_ptm_complet_matrix, cbar_kws = {'label' : 'iptm_ptm'}, xticklabels=True, yticklabels=True)
    ax2.figure.tight_layout()
    plt.savefig("iptm_ptm_heatmap.png")
    plt.close()
    #norm_data_matrix = list()
    #for protein in proteins_list :
    #    new_line = list()
    #    sum_iQ = 0
    #    max_iQ = 0
    #    for iQ_score in save_line[protein] :
    #        sum_iQ += iQ_score
    #        if iQ_score > max_iQ :
    #            max_iQ = iQ_score
    #    ave_iQ = sum_iQ / len(save_line)
    #    for iQ_score in save_line[protein] :
    #        new_line.append(((((iQ_score - ave_iQ) / max_iQ)+1)*0.5)*100) #all values are normalized to the mean and redistribute between 0 and 100.
    #    norm_data_matrix.append(new_line)
    #complet_matrix2 = pd.DataFrame(norm_data_matrix, index = index_prot, columns = index_prot)
    #ax2 = seaborn.heatmap(complet_matrix2, cbar_kws = {'label' : 'iQ_score*'}, xticklabels=True, yticklabels=True)
    #ax2.figure.tight_layout()
    #plt.savefig("Normalized_heatmap.png")
    #plt.close()

###Probably add to generate_interaction_network ???

    ###Probably add to generate_interaction_network ???

def redef_interface (file):
    alphabet = string.ascii_lowercase
    interface_dict = file.get_interface_dict()
    for proteins in interface_dict.keys() :
        interface_dict[proteins] = sorted(interface_dict[proteins]) #sorted all interface in function of number of resiudes
        for interface1 in range(len(interface_dict[proteins])) :
            if interface1 == 0 :
                interface_dict[proteins][interface1].insert(0,alphabet[0])
            for interface2 in range(interface1+1,len(interface_dict[proteins])) :
                list_inter = list(set(interface_dict[proteins][interface1]).intersection(set(interface_dict[proteins][interface2])))
                simi_inter = len(list_inter)/(len(set(interface_dict[proteins][interface1]).union(set(interface_dict[proteins][interface1])))-2) #indice jaccard *100 #-2 just to remove interface 'a' and uniprotID from .union()
                if simi_inter < 0.15 : #if interfaces got more than 0.15 of same residues, it's the same interface
                    interface_dict[proteins][interface2].insert(0,alphabet[interface2])
                else :
                    interface_dict[proteins][interface2].insert(0,interface_dict[proteins][interface1][0])

