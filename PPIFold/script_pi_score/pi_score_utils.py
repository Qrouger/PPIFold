#Adapted from utils.py (https://gitlab.com/topf-lab/pi_score/-/tree/master/score_scripts)

import os
import json
import pandas as pd
import numpy as np
import pickle
import sys

def write_csv_with_features_wc_2(indir,outfle):
    out_csv = outfle
    lst_features = ['Num_intf_residues', 'Polar', 'Hydrophobhic', 'Charged']
    outfle = open(out_csv,'w')
    #writing header in csv file
    outfle.write('pdb,interface,')
    for features in lst_features:
        outfle.write(features + ',')
    outfle.write('contact_pairs')
    outfle.write(',sc,hb,sb,int_solv_en,int_area,pvalue')
    outfle.write('\n')
    for fle in os.listdir(indir):
        if os.path.isdir(os.path.join(indir,fle)):
            intf_dict = ''
            cons = {}
            pisa_dict = {}
            dict_num_contacts = {}
            sc_dict = {}

            for j in os.listdir(os.path.join(indir,fle)):
                nw_path = os.path.join(os.path.join(indir,fle),j)
                if not os.path.isfile(nw_path):
                    continue
                if not nw_path.endswith('json'):
                    continue

                if nw_path.endswith('_interface_chain_contacts.json'):
                    contacts_dict = nw_path
                    dict_num_contacts = count_contacts_at_interface(contacts_dict)

                elif nw_path.endswith('interface_properties_dict.json'):
                    intf_dict = nw_path

                elif nw_path.endswith('sc_scores.json'):
                    sc_dict1 = nw_path
                    with open(sc_dict1,'rb') as handle:
                        sc_dict = json.load(handle)
                elif nw_path.endswith('dict_pisa.json'):
                    with open(nw_path) as handle:
                        pisa_dict = json.load(handle)
                if intf_dict:
                    with open(intf_dict) as handle:
                        dict_intf_prop = json.load(handle)
                        for intf in dict_intf_prop:
                            outfle.write(fle + ',' + intf + ',')
                            for feature in lst_features:
                                if feature == lst_features[0]:
                                    outfle.write(str(dict_intf_prop[intf][feature]) + ',')
                                else:
                                    try:
                                        outfle.write(str(round(dict_intf_prop[intf][feature],3)) + ',')
                                    except:
                                        outfle.write('NA' + ',')
                            intf1 = intf.strip().split('_')[-1] + '_' + intf.strip().split('_')[0]
                            if dict_num_contacts:
                                outfle.write(str(dict_num_contacts[intf]) + ',')
                            else:
                                outfle.write('NA,')
                            if sc_dict:
                                for pdbid_sc in sc_dict.keys():
                                    if intf in sc_dict[pdbid_sc].keys() and len(sc_dict[pdbid_sc][intf]) > 0:
                                        outfle.write(str(sc_dict[pdbid_sc][intf]) + ',')
                                    elif intf1 in sc_dict[pdbid_sc].keys() and len(sc_dict[pdbid_sc][intf1]) > 0:
                                        outfle.write(str(sc_dict[pdbid_sc][intf1]) + ',')
                                    else:
                                        outfle.write('NA,')
                            else:
                                outfle.write('NA,')

                            if pisa_dict:
                                k = ''
                                if intf in pisa_dict.keys():
                                    k = intf
                                elif intf1 in pisa_dict.keys():
                                    k = intf1
                                if k:
                                    try:
                                        outfle.write(str(pisa_dict[k]['hb']) + ',')
                                    except:
                                        outfle.write('NA,')
                                    try:
                                        outfle.write(str(pisa_dict[k]['sb']) + ',')
                                    except:
                                        outfle.write('NA,')
                                    try:
                                        outfle.write(str(pisa_dict[k]['int_solv_en']) + ',')
                                    except:
                                        outfle.write('NA,')
                                    try:
                                        outfle.write(str(pisa_dict[k]['int_area']) + ',')
                                    except:
                                        outfle.write('NA,')
                                    try:
                                        outfle.write(str(pisa_dict[k]['pvalue']))
                                    except:
                                        outfle.write('NA')
                                else:
                                    outfle.write('NA,NA,NA,NA,NA')
                            outfle.write('\n')
    outfle.close()

def write_csv_with_features_wc(indir,outfle):
    out_csv = outfle
    lst_json_fles = []
    lst_features = ['Num_intf_residues', 'Polar', 'Hydrophobhic', 'Charged']
    outfle = open(out_csv,'w')
    #writing header in csv file
    outfle.write('pdb,interface,')
    for features in lst_features:
        outfle.write(features + ',')
    outfle.write('contact_pairs')
    outfle.write(', sc, hb, sb, int_solv_en, int_area, pvalue')
    outfle.write('\n')
    for fle in os.listdir(indir):
        if os.path.isdir(os.path.join(indir,fle)):
            intf_dict = ''
            cons = {}
            pisa_dict = {}
            dict_num_contacts = {}
            complete_sc_dict = {}

            for j in os.listdir(os.path.join(indir,fle)):
                nw_path = os.path.join(os.path.join(indir,fle),j)
                if not os.path.isfile(nw_path):
                    continue
                if not nw_path.endswith('json'):
                    continue

                if nw_path.endswith('_interface_chain_contacts.json'):
                    contacts_dict = nw_path
                    dict_num_contacts = count_contacts_at_interface(contacts_dict)

                elif nw_path.endswith('interface_properties_dict.json'):
                    intf_dict = nw_path

                elif nw_path.endswith('sc_scores.json'):
                    with open(nw_path,'rb') as handle:
                        sc_dict = json.load(handle)
                        for pdbid, intf_dict_sc in sc_dict.items():
                            if pdbid not in complete_sc_dict:
                                complete_sc_dict[pdbid] = {}
                            complete_sc_dict[pdbid].update(intf_dict_sc)

                elif nw_path.endswith('dict_pisa.json'):
                    with open(nw_path) as handle:
                        pisa_dict = json.load(handle)

            if intf_dict:
                with open(intf_dict) as handle:
                    dict_intf_prop = json.load(handle)
                    for intf in dict_intf_prop:
                        outfle.write(fle + ',' + intf + ',')
                        for feature in lst_features:
                            if feature == lst_features[0]:
                                outfle.write(str(dict_intf_prop[intf][feature]) + ',')
                            else:
                                try:
                                    outfle.write(str(round(dict_intf_prop[intf][feature],3)) + ',')
                                except:
                                    outfle.write('NA' + ',')
                        #if cons:
                        #    try:
                        #        outfle.write(str(round(cons[intf],3)) + ',')
                        #    except:
                        #        outfle.write('NA'+ ',')
                        #else:
                        #   outfle.write('NA'+ ',')
                        intf1 = intf.strip().split('_')[-1] + '_' + intf.strip().split('_')[0]
                        if dict_num_contacts:
                            outfle.write(str(dict_num_contacts[intf]) + ',')
                        else:
                            outfle.write('NA,')
                        sc_val = 'NA'
                        for pdbid_sc, intf_sc_dict in complete_sc_dict.items():
                            if intf in intf_sc_dict:
                                sc_val = intf_sc_dict[intf]
                                break
                            elif intf1 in intf_sc_dict:
                                sc_val = intf_sc_dict[intf1]
                                break
                        outfle.write(sc_val+",")

                        if pisa_dict:
                            k = ''
                            if intf in pisa_dict.keys():
                                k = intf
                            elif intf1 in pisa_dict.keys():
                                k = intf1
                            if k:
                                try:
                                    outfle.write(str(pisa_dict[k]['hb']) + ',')
                                except:
                                    outfle.write('NA,')
                                try:
                                    outfle.write(str(pisa_dict[k]['sb']) + ',')
                                except:
                                    outfle.write('NA,')
                                try:
                                    outfle.write(str(pisa_dict[k]['int_solv_en']) + ',')
                                except:
                                    outfle.write('NA,')
                                try:
                                    outfle.write(str(pisa_dict[k]['int_area']) + ',')
                                except:
                                    outfle.write('NA,')
                                try:
                                    outfle.write(str(pisa_dict[k]['pvalue']))
                                except:
                                    outfle.write('NA')
                            else:
                                outfle.write('NA,NA,NA,NA,NA')
                        outfle.write('\n')


def count_contacts_at_interface(contacts_dict):
    dict_num_contacts = {}
    with open(contacts_dict,'rb') as handle:
        dict_contacts = json.load(handle)
    #print dict_contacts
    for intf_name in dict_contacts:
        ct = 0
        for ch1_res in dict_contacts[intf_name].keys():
            for ch2_res in dict_contacts[intf_name][ch1_res]:
                ct += 1
        dict_num_contacts[intf_name] = ct
    return dict_num_contacts

def filter_csv(csvfile,outfile):
    '''
    Takes the CSV file as input and does a sanity check that all features are computed
    '''
    if os.path.isfile(outfile):
        os.remove(outfile)
    out = open(outfile,'w')
    with open(csvfile,'r') as infile:
        for lne in infile:
            if 'NA' not in lne:
                out.write(lne)
    out.close() 

def make_predictions(saved_sc, csvfile, saved_model, outfle, session_id):
    with open(saved_model,'rb') as file:
        clf_m = pickle.load(file)

    topredict = pd.read_csv(csvfile)
    if topredict.empty:
        print 'All features could not be calculated successfully for atleast one of the interfaces!'
        print 'Exiting'
        sys.exit(1)
    topredict_1 = topredict.drop(['pdb','interface','Num_intf_residues'], axis=1)
    with open(saved_sc, 'rb') as file:
        sc = pickle.load(file)
        topredict_1 = np.nan_to_num(sc.transform(topredict_1))
    
    pred_em = clf_m.predict(topredict_1)
    dec = clf_m.decision_function(topredict_1)

    out = open(outfle, 'w')
    out.write('#PDB,chains,predicted_class,pi_score\n')


    pdb_em = topredict['pdb']
    pdb_em_int = topredict['interface']

    for i in range(len(topredict_1)):
        out.write(str(pdb_em[i]) + ',' +
                  str(pdb_em_int[i]) + ',' +
                  str(pred_em[i]) + ',' +
                  str(round(dec[i],2)) + '\n')
    out.close()
