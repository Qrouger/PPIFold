#Adapted from run_piscore_wc.py (https://gitlab.com/topf-lab/pi_score/-/tree/master/score_scripts)
#code in python 2.7
import argparse
import os
import shutil
from interface_assess import *
from sc_utils import write_sc_script, run_sc
from clean_pdb import write_sel_pdb, change_atloc
from pisa_utils import run_pisa, parse_pisa_xml_outfle
from pi_score_utils import write_csv_with_features_wc, make_predictions, filter_csv
import time
import sys
import uuid 

author = 'Author: Sony Malhotra '


parser = argparse.ArgumentParser(description="Calculate PI-score for a given complex structure")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-p", "--pdb", help="PDB file to calculate the PI-score", default=None)
group.add_argument("-d", "--dirc", help="directory to assess the PDB files",default=None)
parser.add_argument("-ch", "--chains", help="Chains forming the interface which needs to be assessed",default=None)
parser.add_argument("-dc", "--distance_cutoff", help="Distance cutoff for interface definition",default=7, type =float)
parser.add_argument("-s", "--intf_size", help="Number of interface residues from each of the interacting subunits",default=10,type=int)
parser.add_argument("-ps", "--prot_size", help="Length of each of the chains",default=30,type=int)
parser.add_argument("-np", "--num_processors", help="Number of processors to use",default=40,type=int)
parser.add_argument("-o", "--out", help="Name of the output directory to write results",default=os.path.join(os.getcwd(),'pi_output'))
parser.add_argument("-c","--csv",help="Name of the CSV files where features will be stored",default=None)
parser.add_argument("-m","--model",help="SVM model to make new predictions",choices=['A','B','WC'], default='WC')
parser.add_argument("-r","--results",help="Name of the file where model predictions will be written",default=None)

def main(): 
    args = parser.parse_args()
    pdblist = []
    
    if args.dirc != None:
        #print ('Setting chains to all chains')
        args.chains = None

    #print ('*****Selecting model for training*****')
    dirname, filename = os.path.split(os.path.abspath(__file__))
    if args.model == 'A':
        saved_model = os.path.join(os.path.dirname(sys.argv[0]),'svm_model/finalized_model_A.sav')
        saved_sc = os.path.join(dirname,'svm_model/scaler_model_A.sav')
    elif args.model == 'B':
        saved_model = os.path.join(os.path.dirname(sys.argv[0]),'svm_model/finalized_model_B.sav')
        saved_sc = os.path.join(dirname,'svm_model/scaler_model_B.sav')
    if args.model == 'WC':
        saved_model = os.path.join(os.path.dirname(sys.argv[0]),'svm_model/finalized_model_wc.sav')
        saved_sc = os.path.join(dirname,'svm_model/scaler_model_wc.sav')
    working_path = os.getcwd()

    #print ('*****Setting output directory for results*****')
    if not os.path.isdir(args.out):
        os.mkdir(args.out)

    session_id = uuid.uuid4().hex[:8]

    if args.csv is None:
        current_time = time.strftime("%m.%d.%y_%H%M", time.localtime())
        obj = session_id
        args.csv = os.path.join(args.out, "intf_features_%s_%s.csv" % (current_time, session_id))
        
    if args.results is None:
        current_time = time.strftime("%m.%d.%y_%H%M", time.localtime())
        args.results = os.path.join(args.out,"pi_score_%s_%s_id.txt" % (current_time, session_id))
    
    if os.path.isfile(args.results):
        os.remove(args.results)
    
    out = open(args.results,'w')
    out.write('#PDB,chains,predicted_class,pi_score\n')
    out.close()
    
    filter_csvfle = os.path.join(args.out,"filter_intf_features_%s_%s.csv" % (current_time, session_id))

    #print ('*****Setting input for assessment*****')
    #get a list of pdbfiles to score interfaces for
    if args.pdb:
        pdblist = [args.pdb.rsplit('/',1)[-1]]
    elif args.dirc:
        for fl in os.listdir(args.dirc):
            if fl.endswith('pdb'):
                pdblist.append(fl)
    else:
        print ('Not a valid input')


    #list of directories where to run conservation
    #dir_infiles = []
    #all_conslst = []
    #all_consdir = []
    #ch_lst = []

    #print ('*****Iterating over each structure to assess*****')
    for pdb in pdblist:
        chains_in = []
        os.chdir(working_path)
        pd = pdb.rsplit('.',1)[0]
        #make an output directory
        if os.path.isdir(os.path.join(args.out,pd)):
            shutil.rmtree(os.path.join(args.out,pd))
        #print ('*****Making results subdirectory specific to each structure in output directory*****')
        #print (os.path.join(args.out,pd))
        os.mkdir(os.path.join(args.out,pd))
        #copy input files
        if args.dirc:
            copy_cmd = 'cp '+ os.path.join(args.dirc,pdb) + ' ' + os.path.join(args.out,pd)
        else:
            copy_cmd = 'cp '+ os.path.join(args.pdb) + ' ' + os.path.join(args.out,pd)
        #print ('*****Doing PDB,',pdb,' *****')
        os.system(copy_cmd)
        os.chdir(os.path.join(args.out,pd))
        #print ('*****Calculating interface*****')
        #Interface calculation
        dict_intf,dict_contact = interface_residues(pdb,
                                                    len_chains=args.prot_size,
                                                    numres_cut=args.intf_size,
                                                    dist_cutoff=args.distance_cutoff)
        intf_outfle_prefix = pd
        if not dict_intf: 
            intf_outfle_name = intf_outfle_prefix + '_no_int.txt'
            a = open(intf_outfle_name,'w')
            a.write('No interfaces found')
            a.close()
            continue
        
        #print ('*****Calculating contacts*****')
        #Writing the contact dictionary, which residue from which chains are in contact
        contact_name = intf_outfle_prefix + '_interface_chain_contacts.json'
        with open(contact_name,'w') as outfle:
            json.dump(dict_contact,outfle)
        if args.chains:
            for c in args.chain:
                chains_in.append(c)
        dict_intf_prop = {}
        
        #Write clean pdbfile for shape complementarity
        clean_pdbfile = write_sel_pdb(pdb)
        atm_file = change_atloc(clean_pdbfile)
        os.system('rm -f ' + clean_pdbfile)
        os.system('mv ' + atm_file + ' ' + clean_pdbfile)
        
        #print ('*****Iterating over interfaces in structure*****')
        #Iterating over all the interfaces
        index = 0
        for k in dict_intf:
            index+=1
            ch1 = k.split('_')[0]
            ch2 = k.split('_')[-1]
            flag = False
            if args.chains:
                if ch1 in args.chains and ch2 in args.chains:
                    flag = True
                    break
                else:
                    print ('No interfaces between the input chains,', args.chains)
                    continue
            #print ('*****Calculating interface features*****')
            if flag or not args.chains: 
                ch1_intf_res = dict_intf[k][0]
                ch2_intf_res = dict_intf[k][1]
                dict_intf_prop[k] = {'Intfresidues_ch1':[ch1,len(ch1_intf_res)]}
                dict_intf_prop[k] = {'Intfresidues_ch2':[ch2,len(ch2_intf_res)]}
                total_intf_res = ch1_intf_res + ch2_intf_res
                total_num_intf_res = len(total_intf_res)
                only_res_name_lst = []
                for re in total_intf_res:
                    only_res_name_lst.append(re[0:3])
                polar = polar_residues(only_res_name_lst)
                hydro = hydrophobic_residues(only_res_name_lst)
                charged = charged_residues(only_res_name_lst)
                if dict_intf_prop[k]:
                    dict_intf_prop[k].update({'Polar':polar, 'Charged': charged, 'Hydrophobhic':hydro, 'Num_intf_residues': total_num_intf_res})
                else:
                    dict_intf_prop[k] = {'Polar':polar, 'Charged': charged, 'Hydrophobhic':hydro, 'Num_intf_residues': total_num_intf_res}
                intf_outfle_name = intf_outfle_prefix + '_interface_properties_dict.json'
                with open(intf_outfle_name,'w') as outfle:
                    json.dump(dict_intf_prop,outfle)     

                #print ('*****Calculating shape complementarity*****')
                #calculate shape complementarity
                scriptfile = write_sc_script(chain_lst=k.split('_'),
                                            pdbfile=clean_pdbfile)
                sc_dir = args.out+"/"+pd+"/"+scriptfile
                run_sc(sc_dir)
                #parse_sc_output('tmp_sc.out')
                
        #print ('*****Calculating pisa *****')
        #calculate features using pisa
        outfle_name = run_pisa(clean_pdbfile)
        pisa_out = parse_pisa_xml_outfle(outfle_name)

    os.chdir(working_path)

    #print ('*****Writing CSV file with features *****')
    # # Write the CSV file with interface features
    
    write_csv_with_features_wc(indir=args.out, outfle=args.csv)
    # 
    # 
    # #filter CSV and have only interfaces where all features were computed successfully
    filter_csv(args.csv,filter_csvfle)

    #print ('*****Calculating PI-score*****')
    # # #Make predictions
    make_predictions(saved_sc=saved_sc, 
                    csvfile=filter_csvfle, 
                    saved_model=saved_model,
                    outfle=args.results,
                    session_id=session_id)
if __name__ == '__main__':
    main()

