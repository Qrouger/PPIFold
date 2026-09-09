#!/usr/bin/env python3
#Adapted from get_good_inter_pae.py (https://github.com/KosinskiLab/AlphaPulldown/blob/main/alphapulldown/analysis_pipeline/alpha_analysis_jax0.4.def)

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "script_pi_score"))
import pickle
import json
import shutil
import gemmi
import logging
import subprocess
import gzip
import numpy as np
import pandas as pd
from pathlib import Path
from calculate_mpdockq import *


def examine_inter_pae(pae_mtx, length, cutoff, type_int) :
    """Check inter-chain PAE only between the last chain and the others"""
    pae = pae_mtx.copy()
    if type_int == "PPI" or type_int == "Compounds" :
        start_last = sum(length[:-1])
        # mask all 
        pae[:] = 50

        # restore only inter-PAE with prey
        pae[start_last:, :start_last] = pae_mtx[start_last:, :start_last]
        pae[:start_last, start_last:] = pae_mtx[:start_last, start_last:]
        check = np.where(pae < cutoff)[0].size != 0
    if type_int == "homo" :
        start = 0

        for l in length:
            end = start + l
            pae[start:end, start:end] = 50  # masque intra-chaîne
            start = end

        check = np.any(pae < cutoff)
    return check


def obtain_chain_coord(work_dir,pkl_dict=None) :
    """Returns chain_coords,chain_CB_inds,plddt_per_chain,best_plddt,pdb_path"""
    interaction = work_dir.split("/")[-1]
    pdb_path = os.path.join(work_dir,f'{interaction}_ranked_0.pdb')
    pdb_chains, chain_coords, chain_CA_inds, chain_CB_inds = read_pdb(pdb_path)
    if pkl_dict ==  None :
        best_plddt = extract_plddt_from_pdb(pdb_path)
    else :
        best_plddt = pkl_dict['plddt'] 
    plddt_per_chain = read_plddt(best_plddt,chain_CA_inds)
    return chain_coords,chain_CB_inds,plddt_per_chain,best_plddt,pdb_path

def obtain_mpdockq2(chain_coords,chain_CB_inds,plddt_per_chain,best_plddt,pdb_path) :
    """Returns mpDockQ if more than two chains otherwise return pDockQ"""
    complex_score,num_chains = score_complex(chain_coords,chain_CB_inds,plddt_per_chain)
    if complex_score is not None and num_chains>2:
        mpDockq_or_pdockq = calculate_mpDockQ(complex_score)
    elif complex_score is not None and num_chains==2:
        chain_coords,plddt_per_chain = read_pdb_pdockq(pdb_path)
        mpDockq_or_pdockq = calc_pdockq(chain_coords,plddt_per_chain,t=8)
    else:
        mpDockq_or_pdockq = "None"
    return mpDockq_or_pdockq


def extract_plddt_from_pdb(pdb_file) :
    """
    Extract plddt from b-factor in pdb
    """
    plddt_values = []
    with open(pdb_file, "r") as f :
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA" :
                try :
                    b_factor = float(line[60:66].strip())
                    plddt_values.append(b_factor)
                except ValueError :
                    pass
    return np.array(plddt_values, dtype=float)

def run_and_summarise_pi_score(jobs, surface_thres, ccp4_setup) :
    """
    A function to calculate all predicted models' pi_scores and make a pandas df of the results.
    Instrumented to log timing per major step.
    """
    output_df = pd.DataFrame()
    for t in ["micromamba", "mamba", "conda"] :
        if shutil.which(t) :
            tool = t
    output_df = pd.DataFrame()
    for job in jobs :
        direc = os.path.dirname(job)
        file_pdb = job.split("/")[-1]
        name_job = f"{file_pdb.split('.pdb')[0]}"
        if os.path.isdir("/scratch") :
            tmp_dir = f"/scratch/tmp/{name_job}"
        else :
            tmp_dir = f"/tmp/{name_job}"

        logging.info(f"Creating temporary directory {tmp_dir} for pi_score outputs")
        subprocess.run(f"rm -rf {tmp_dir} && mkdir -p {tmp_dir}/pi_score_outputs", shell=True, executable="/bin/bash", check=True)
        pi_score_outputs = os.path.join(tmp_dir, "pi_score_outputs")
        
        cwd = os.path.dirname(os.path.abspath(__file__))

        if not os.path.isfile(os.path.join(direc, f"{file_pdb}")) :
            logging.error(f"{job} failed. Cannot find {file_pdb} in {direc}")
            sys.exit()


        output_dir = os.path.join(pi_score_outputs)

        cmd = (
            f"source {ccp4_setup}/bin/ccp4.setup-sh && "
            f"{tool} run -n pi_score python {cwd}/script_pi_score/run_piscore_wc.py "
            f"-p {job} -o {output_dir} -s {surface_thres} -ps 10"
        )

        proc = subprocess.Popen(cmd, shell=True, executable="/bin/bash", close_fds=True)
        proc.wait()

   

        subdir = os.path.join(pi_score_outputs)
        csv_files = [f for f in os.listdir(subdir) if 'filter_intf_features' in f]
        pi_score_files = [f for f in os.listdir(subdir) if 'pi_score_' in f]

        if not csv_files or not pi_score_files :
            logging.info(f"Warning: missing CSV or pi_score files for {name_job}")
            continue
        for csv_f in csv_files:
            filtered_df = pd.read_csv(os.path.join(subdir, csv_f))

            if filtered_df.shape[0] == 0 :
                for column in filtered_df.columns :
                    filtered_df[column] = ["None"]
                filtered_df['jobs'] = str(name_job)
                filtered_df['pi_score'] = "No interface detected"
                output_df = pd.concat([output_df, filtered_df])
                continue

            interface_id = csv_f.split("filter_intf_features_")[-1].replace(".csv", "")

            pi_f = [f for f in pi_score_files if interface_id in f]
            if not pi_f :
                logging.info(f"Warning: no pi_score file for interface {interface_id}")
                continue

            pi_score = pd.read_csv(os.path.join(subdir, pi_f[0]))
            pi_score['jobs'] = str(name_job)


            if 'chains' in pi_score.columns :
                pi_score['interface'] = pi_score['chains']
            else:
                pi_score['interface'] = interface_id

            filtered_df['jobs'] = str(name_job)
            last_chain = get_last_chain_from_pdb(os.path.join(job))

            filtered_df = filtered_df[filtered_df['interface'].apply(lambda x: last_chain in x)]
            pi_score = pi_score.drop(columns=["#PDB", "pdb", "pvalue", "chains", "predicted_class"],errors="ignore")

            filtered_df = pd.merge(
                filtered_df,
                pi_score,
                on=['jobs', 'interface'],
                how='left'
            )
            output_df = pd.concat([output_df, filtered_df])

    subprocess.run(f"rm -rf {tmp_dir}", shell=True, executable='/bin/bash')
    return output_df
    
def get_last_chain_from_pdb(pdb_file) :
    """
    Get the last chain identifier from a PDB file. This is used to filter pi_score results for the last chain only.
    """
    last_chain = None
    with open(pdb_file, 'r') as f :
        for line in f :
            if line.startswith(('ATOM', 'HETATM')):
                chain = line[21]  # col 22 (index 21)
                last_chain = chain
    return last_chain


def main(job, cutoff, surface_thres, save_file, AF_version, ccp4_setup) :
    """
    Main function to check inter-PAE and run pi_score for good models. Returns a pandas df with all scores.
    structure confidence (ipTM)
    geometry confidence (PAE)
    interface quality (mpDockQ)
    physico-chemical scoring (PI-score)

    Parameters :
    ----------
    job : str
    cutoff : int
    surface_thres : int
    save_file : File_proteins object
    AF_version : str
    ccp4_setup : str
    """
    seq_no_SP = save_file.get_proteins_sequence()
    prot_length = save_file.get_lenght_prot()
    good_jobs = []
    iptm_ptm = list()
    iptm = list()
    name_jobs = list()
    mpDockq_scores = list()
    logging.info(f"Scoring {job}")
    result_subdir = os.path.join(job)
    interaction = job.split("/")[-1]
    length = list()
    if "_and_" in interaction :
        type_int = "PPI"
        for prot in interaction.split("_and_") :
            if "-" in prot and prot.split("_")[0] in seq_no_SP.keys() :
                prot = prot.split("_")[0]
            length.append(prot_length[prot])
    if "_homo_" in interaction :
        type_int = "homo"
        prot = interaction.split("_homo_")[0] 
        nbr = int(interaction.split("_homo_")[1].replace("er",""))
        for i in range(0,nbr) :
            if "-" in prot and prot.split("_")[0] in seq_no_SP.keys() :
                prot = prot.split("_")[0]
            length.append(prot_length[prot])
    if "Compounds" in job :
        type_int = "Compounds"
        for prot in interaction.split("_and_") :
            if "-" in prot and prot.split("_")[0] in seq_no_SP.keys() :
                prot = prot.split("_")[0]
            if prot == interaction.split("_and_")[-1] :
                l = 1
            else :
                l = prot_length[prot]
            length.append(l)

            

    ### AlphaFold3 ###
    if AF_version == "3" :
        cif_files = list(Path(result_subdir).glob("*model.cif"))
        pdb = "ranked_0.pdb"
        if os.path.isfile(os.path.join(result_subdir,f'{interaction}_ranked_0.pdb')) == False  : #create ranked_0.pdb for AF3 
            doc = gemmi.cif.read_file(str(cif_files[0]))
            block = doc.sole_block()
            structure = gemmi.make_structure_from_block(block)
            structure.write_pdb(os.path.join(result_subdir, f'{interaction}_ranked_0.pdb'))

        int_AF3 = str(cif_files[0]).split("_model")[0]
        if os.path.isfile(int_AF3+'_summary_confidences.json') :
            with open(int_AF3+'_summary_confidences.json','rb') as json_sum_f :
                json_sum = json.load(json_sum_f)
            if "iptm" in json_sum.keys() and "ptm" in json_sum.keys() :
                iptm_score = json_sum['iptm']
                ptm_score = json_sum['ptm']
                iptm_ptm_score = 0.8 * iptm_score + 0.2 * ptm_score
            with open(int_AF3+'_confidences.json','rb') as json_f :
                json_data = json.load(json_f)
            pae_mtx = np.array(json_data['pae'])
            chain_coords,chain_CB_inds,plddt_per_chain,best_plddt,pdb_path = obtain_chain_coord(os.path.join(job))
            check = examine_inter_pae(pae_mtx,length,cutoff=cutoff,type_int=type_int)
            mpDockq_score = obtain_mpdockq2(chain_coords,chain_CB_inds,plddt_per_chain,best_plddt,pdb_path)
            if check:
                good_jobs.append(str(job)+"/"+f"{job.split('/')[-1]}_ranked_0.pdb")
                iptm_ptm.append(iptm_ptm_score)
                iptm.append(iptm_score)
                mpDockq_scores.append(mpDockq_score)
                name_jobs.append(f"{job.split('/')[-1]}_{pdb.split('.')[0]}")
        else :
            logging.info(f"Cannot find summary_confidences.json for {job}, skipping.")


    ### AlphaFold2 ###
    if os.path.isfile(os.path.join(job,'ranking_debug.json')) :
        mod_index = -1
        os.system(f"cp {job}/ranked_0.pdb {job}/{interaction}_ranked_0.pdb") #rename pdb file with explicit name
        all_models = ["ranked_0.pdb"]
        with open(os.path.join(result_subdir,'ranking_debug.json'),'rb') as json_f :
            data = json.load(json_f)
        best_model = data['order']
        for pdb in all_models :
            mod_index += 1
            rank_model = best_model[mod_index]
            if "iptm" in data.keys() or "iptm+ptm" in data.keys() :
                iptm_ptm_score = data['iptm+ptm'][rank_model]
                if os.path.exists(os.path.join(result_subdir, f"result_{rank_model}.pkl")) :
                    pkl_path = os.path.join(result_subdir, f"result_{rank_model}.pkl")
                    with open(pkl_path, 'rb') as pkl :
                        check_dict = pickle.load(pkl)
                elif os.path.exists(os.path.join(result_subdir, f"result_{rank_model}.pkl.gz")) :
                    logging.info("Result pickle for the best model not found. Now search for zipped pickle.")
                    pkl_path = os.path.join(result_subdir, f"result_{rank_model}.pkl.gz")
                    with gzip.open(pkl_path, 'rb') as pkl :
                        check_dict = pickle.load(pkl)
                else :
                    logging.info(f"Cannot find result pickle for {job}, skipping.")
                iptm_score = check_dict['iptm']
                pae_mtx = np.array(check_dict['predicted_aligned_error'])
                chain_coords,chain_CB_inds,plddt_per_chain,best_plddt,pdb_path = obtain_chain_coord(os.path.join(job),check_dict)
                if pdb == "ranked_0.pdb" :
                    check = examine_inter_pae(pae_mtx,length,cutoff=cutoff,type_int=type_int) #only check PAE for best model
                mpDockq_score = obtain_mpdockq2(chain_coords,chain_CB_inds,plddt_per_chain,best_plddt,pdb_path)
                if check :
                    good_jobs.append(str(f"{job}/{job.split('/')[-1]}_{pdb}"))
                    iptm_ptm.append(iptm_ptm_score)
                    iptm.append(iptm_score)
                    mpDockq_scores.append(mpDockq_score)
                    name_jobs.append(f"{job.split('/')[-1]}_{pdb.split('.')[0]}")
        del data
        
    other_measurements_df=pd.DataFrame.from_dict({
        "jobs":name_jobs,
        "iptm_ptm":iptm_ptm,
        "iptm":iptm,
        "mpDockQ/pDockQ":mpDockq_scores})

    if good_jobs!=[] :
        pi_score_df = run_and_summarise_pi_score(good_jobs, surface_thres, ccp4_setup)
        pi_score_df = pd.merge(pi_score_df, other_measurements_df, on="jobs")
        columns = list(pi_score_df.columns.values)
        columns.pop(columns.index('jobs'))
        pi_score_df = pi_score_df[['jobs'] + columns]
        pi_score_df = pi_score_df.sort_values(by='iptm',ascending=False)
        
        return pi_score_df
    else :
        return None



