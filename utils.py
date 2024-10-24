import numpy as np
import pickle
import csv
import matplotlib.pyplot as plt
from Bio import PDB
import copy
import glob
import logging
from scipy.special import softmax
from File_proteins import *
import os
import networkx as nx
from math import *
import seaborn
import pandas as pd

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

def create_feature (file, env_feature, data_dir) :
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
        cmd = f"#!/bin/bash --login \n source ~/.bashrc \n conda activate {env_feature}\n create_individual_features.py --fasta_paths=./{fasta_file} \--data_dir={data_dir} \--save_msa_files=True \--output_dir=./feature \--max_template_date=2024-05-02 \--skip_existing=False"
        cmd2 = f"create_individual_features.py --fasta_paths=./{fasta_file} \--data_dir={data_dir} \--save_msa_files=True \--output_dir=./feature \--max_template_date=2024-05-02 \--skip_existing=False"
        cmd3 = "#!/bin/bash --login \n source ~/.bashrc \n conda deactivate"
        if env_feature != None :
            os.system(cmd)
            os.system(cmd3)
        else :
            os.system(cmd2)

def Make_all_MSA_coverage (file) :
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
            pre_feature_dict = pickle.load(open(f'feature/{prot}.pkl','rb'))
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
            plt.savefig(f"feature/{prot+('_' if prot else '')}coverage.pdf")
            plt.close()
        with open("bad_MSA.txt", "w") as MSA_file :
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

def Make_all_vs_all (env_multimers, data_dir) :
        """
        Use Alphapulldown script to generate all versus all interactions.

        Parameters:
        ----------
        env_multimers : string
        data_dir : string
        Returns:
        ----------

        """
        cmd = f"#!/bin/bash --login \n source ~/.bashrc \n conda activate {env_multimers}\n run_multimer_jobs.py --mode=custom \--num_cycle=3 \--num_predictions_per_model=1 \--compress_result_pickles=True \--output_path=result_all_vs_all \--data_dir={data_dir} \--protein_lists=all_vs_all.txt \--monomer_objects_dir=./feature"
        cmd2 = "run_multimer_jobs.py --mode=all_vs_all \--num_cycle=3 \--num_predictions_per_model=1 \--compress_result_pickles=True \--output_path=./result_all_vs_all \--data_dir={data_dir} \--protein_lists=all_vs_all.txt \--monomer_objects_dir=./feature"
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
        for interaction in iQ_score_dict.keys() :
            if float(iQ_score_dict[interaction]) >= 35 : #Plot figure of interest just for interesting interactions
                job1 = interaction[0] + "_and_" + interaction[1]
                plot_Distogram("./result_all_vs_all/" + job1)
                make_table_res_int("./result_all_vs_all/" + job1)
        hiQ_score_dict = file.get_hiQ_score_dict()
        for homo_oligo in hiQ_score_dict.keys() :
            if float(hiQ_score_dict[homo_oligo][0]) >= 50 :
                job2 = homo_oligo + "_homo_" + hiQ_score_dict[homo_oligo][1] + "er"
                plot_Distogram("./result_homo_oligo/" + job2)
                make_table_res_int("./result_homo_oligo/" + job2)

def make_table_res_int (int) :
        """
        Generate a table of residue in interactions.

        Parameters:
        ----------
        int : string
        Returns:
        ----------

        """
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', int + "/ranked_0.pdb")
        dict_int = dict()
        dict_information = dict()
        atom_possible_contact = ["O","OH","NH2","NH1","OG","NE2","ND2","NZ","NE","N","OE1","OE2","OD2","OG1"]
        for model in structure:
            list_chain = model.get_list()
            print(str(len(list_chain))+" chains detected")
            for i in range(0,len(list_chain)) :
                chain1 = list_chain[i] #B and C
                for residue1 in chain1 : #number of residue (len of the chain)
                    for atom1 in residue1 : #type of the atom
                        for index in range(i+1,len(list_chain)) :
                            chain2 = list_chain[index]
                            for residue2 in chain2 :
                                for atom2 in residue2 :
                                    distance = atom1 - atom2
                                    if distance <= 3.64 :
                                        if chain1.get_id()+chain2.get_id() in dict_int.keys() : #to make different table for different interaction  
                                            dict_int[chain1.get_id()+chain2.get_id()].append([chain1.get_id()+":"+residue1.get_resname()+"  "+str(residue1.get_id()[1]),chain2.get_id()+":"+residue2.get_resname()+"   "+str(residue2.get_id()[1]),str(distance)])
                                            dict_information[chain1.get_id()+chain2.get_id()].append([residue1.get_resname()+str(residue1.get_id()[1]),atom1.get_id(),residue2.get_resname()+str(residue2.get_id()[1]),atom2.get_id(),str(distance)])
                                        else :
                                            dict_int[chain1.get_id()+chain2.get_id()] = [["Chain "+chain1.get_id(),"Chain "+chain2.get_id(),"Distance Ä"]]
                                            dict_information[chain1.get_id()+chain2.get_id()] = [[" "," ","Chain "+chain1.get_id(),"Chain "+chain2.get_id(),"Distance Ä"]]
                                            dict_int[chain1.get_id()+chain2.get_id()].append([chain1.get_id()+":"+residue1.get_resname()+"   "+str(residue1.get_id()[1]),chain2.get_id()+":"+residue2.get_resname()+"   "+str(residue2.get_id()[1]),str(distance)])
                                            dict_information[chain1.get_id()+chain2.get_id()].append([residue1.get_resname()+str(residue1.get_id()[1]),atom1.get_id(),residue2.get_resname()+str(residue2.get_id()[1]),atom2.get_id(),str(distance)])
                                    else :
                                        pass

        int_list = list()
        save_dict = copy.deepcopy(dict_int)
        for chains in dict_int.keys() :
            for line in range(len(save_dict[chains])) :
                if dict_information[chains][line][1] not in atom_possible_contact or dict_information[chains][line][3] not in atom_possible_contact and dict_information[chains][line][0] != " " : #sort in function of atom id
                    for index in range(len(dict_int[chains])) :
                        if dict_int[chains][index-1] == save_dict[chains][line] :
                            dict_int[chains].pop(index-1)
        save_dict2 = copy.deepcopy(dict_int)
        for chains2 in save_dict2.keys() : #sort in function of the distance
            for line2 in range(len(save_dict2[chains2])) :
                int_list.append(save_dict2[chains2][line2])
        for chains2 in save_dict2.keys() : #sort in function of the distance
            for line2 in range(len(save_dict2[chains2])) :
                for interaction in range(len(int_list)) :
                    if int_list[interaction][0] == save_dict2[chains2][line2][0] and int_list[interaction][1] == save_dict2[chains2][line2][1] and float(save_dict2[chains2][line2][2]) > float(int_list[interaction][2]) :
                        for index in range(len(dict_int[chains2])) :
                            if dict_int[chains2][index-1] == save_dict2[chains2][line2] :
                                dict_int[chains2].pop(index-1)
                    elif int_list[interaction][0] == save_dict2[chains2][line2][0] and int_list[interaction][1] == save_dict2[chains2][line2][1] and float(save_dict2[chains2][line2][2]) < float(int_list[interaction][2]) :
                        for index in range(1,len(dict_int[chains2])) :
                            if dict_int[chains2][index-1] == int_list[interaction] :
                                dict_int[chains2].pop(index-1)
                    else :
                        pass
            fileout = chains2+"_res_int.csv"
            np_table = np.array(dict_int[chains2])
            with open(f"{int}/"+fileout, "w", newline="") as file :
                mywriter = csv.writer(file, delimiter=",")
                mywriter.writerows(np_table)
            print("Write table")

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

def Make_homo_oligo (env_multimers, data_dir) :
        """
        Use Alphapulldown script to generate all homo oligomer.

        Parameters:
        ----------
        env_multimers : string
        data_dir : string
        
        Returns:
        ----------

        """
        cmd = f"#!/bin/bash --login \n source ~/.bashrc \n conda activate {env_multimers}\n run_multimer_jobs.py --mode=homo-oligomer \--output_path=result_homo_oligo \--num_cycle=3 \--compress_result_pickles=True \--oligomer_state_file=homo_oligo.txt \--monomer_objects_dir=feature \--data_dir={data_dir} \--remove_result_pickles=False"
        cmd2 = "run_multimer_jobs.py --mode=homo-oligomer \--output_path=result_homo_oligo \--num_cycle=3 \--compress_result_pickles=True \--oligomer_state_file=homo_oligo.txt \--monomer_objects_dir=feature \--data_dir={data_dir} \--remove_result_pickles=False"
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
        cmd4 = f"singularity exec --no-home --bind result_homo_oligo:/mnt {dir_alpha}/alpha-analysis_jax_0.4.sif run_get_good_pae.sh --output_dir=/mnt --cutoff=10"
        os.system(cmd4)
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
    #for homo_oligomer in hiQ_score_dict.keys() :
    #    if float(hiQ_score_dict[homo_oligomer][0]) >= 50 :
    #        valid_interactions.append([[homo_oligomer,homo_oligomer], hiQ_score_dict[homo_oligomer][1]])
    int_graph = nx.Graph()
    list_inter_score = list()
    prots = set()
    dict_name = file.get_names()
    for inter, score in valid_interactions :
        inter0 = inter[0]+f"({dict_name[inter[0]]})" #set uniprotID with the name of protein
        inter1 = inter[1]+f"({dict_name[inter[1]]})"
        prots.add(inter0)
        prots.add(inter1)
        list_inter_score.append((inter0,inter1,round(score,2)))
    prots = list(prots)
    int_graph.add_nodes_from(prots)
    int_graph.add_weighted_edges_from(list_inter_score)
    pos = nx.spring_layout(int_graph, k = len(prots)+5, scale = 3, seed=7)
    nx.draw_networkx_nodes(int_graph,pos)
    nx.draw_networkx_edges(int_graph,pos, node_size=100, edgelist=int_graph.edges, width=1, style="dashed")
    nx.draw_networkx_labels(int_graph, pos, font_size=10, font_family="sans-serif", hide_ticks = 'True')
    edge_labels = nx.get_edge_attributes(int_graph, "weight")
    nx.draw_networkx_edge_labels(int_graph, pos, edge_labels)
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
    data_matrix = list()
    iQ_score_dict = file.get_iQ_score_dict()
    proteins_list = file.get_proteins()
    proteins_name = file.get_names()
    index_prot = list()
    for protein in proteins_list :
        index_prot.append(protein+"_"+proteins_name[protein])
    save_line = dict()
    for protein1 in proteins_list :
        line = [] 
        for protein2 in proteins_list :
            if protein1 == protein2 :
                line.append(0)
            else :
                if (protein1,protein2) in iQ_score_dict.keys() :
                    line.append(float(iQ_score_dict[(protein1,protein2)]))
                elif (protein2,protein1) in iQ_score_dict.keys() :
                    line.append(float(iQ_score_dict[(protein2,protein1)]))
                else :
                    line.append(float(0))
                    print (protein1 + " and " + protein2 + " are not in score table or have bad PAE")
        save_line[protein1] = line
        data_matrix.append(line)
    complet_matrix = pd.DataFrame(data_matrix,index = index_prot, columns = index_prot )
    ax2 = seaborn.heatmap(complet_matrix, cbar_kws = {'label' : 'iQ_score'}, xticklabels=True, yticklabels=True)
    ax2.figure.tight_layout()
    plt.savefig("heatmap.png")
    plt.close()
    norm_data_matrix = list()
    for protein in proteins_list :
        new_line = list()
        sum_iQ = 0
        max_iQ = 0
        for iQ_score in save_line[protein] :
            sum_iQ += iQ_score
            if iQ_score > max_iQ :
                max_iQ = iQ_score
        ave_iQ = sum_iQ / len(save_line)
        for iQ_score in save_line[protein] :
            print(protein)
            print(max_iQ)
            new_line.append(((((iQ_score - ave_iQ) / max_iQ)+1)*0.5)*100) #all values are normalized to the mean and redistribute between 0 and 100.
        norm_data_matrix.append(new_line)
    complet_matrix2 = pd.DataFrame(norm_data_matrix, index = index_prot, columns = index_prot)
    ax2 = seaborn.heatmap(complet_matrix2, cbar_kws = {'label' : 'iQ_score*'}, xticklabels=True, yticklabels=True)
    ax2.figure.tight_layout()
    plt.savefig("Normalized_heatmap.png")
    plt.close()