""" Main file of PPIFold

    Author: Quentin Rouger
"""
import argparse

from .Utils import *
from .File_proteins import *
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'AlphaPulldown')))

def add_arguments(parser) :
    parser.add_argument("--use_mmseq", help = "Use mmseq for feature generation", required = False, default = True)
    parser.add_argument("--make_multimers", help = "If you just want make feature set on False", required = False, default = True)
    parser.add_argument("--max_aa" , help = "Maximum amino acids can be generate by your cluster", required = False, default = 2500, type = int)
    parser.add_argument("--use_signalP" , help = "Used or not SignalP", required = False, default = True)
    parser.add_argument("--org" , help = "Organism of interest : arch, gram+, gram- or euk", required = False, default = "gram-", type = str)

def main() :
    sys.stdout = open('./PPI.log', 'w')
    path_dict = define_path()
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args()
    PPI_object = File_proteins(path_dict["Path_Uniprot_ID"])
    PPI_object.find_proteins_sequence()
    PPI_object.find_prot_lenght()
    if len(PPI_object.already_pickle(path_dict["Path_Pickle_Feature"])) > 0 : #if new feature pickle is need
        PPI_object.create_fasta_file()
        if args.use_signalP == True :
            remove_SP(PPI_object,args.org)
        create_feature(PPI_object,path_dict["Path_AlphaFold_Data"],path_dict["Path_Pickle_Feature"],args.use_mmseq)
    Make_all_MSA_coverage(PPI_object,path_dict["Path_Pickle_Feature"])
    generate_APD_script(PPI_object, args.max_aa)
    if args.make_multimers == True :
        Make_all_vs_all(path_dict["Path_AlphaFold_Data"],path_dict["Path_Pickle_Feature"])
        add_iQ_score(path_dict["Path_Singularity_Image"])
        Make_homo_oligo(path_dict["Path_AlphaFold_Data"],path_dict["Path_Pickle_Feature"])
        add_hiQ_score(path_dict["Path_Singularity_Image"])
        PPI_object.update_iQ_score_hiQ_score()
        generate_heatmap(PPI_object)
        create_out_fig(PPI_object)
        generate_interaction_network(PPI_object)
        plot_sequence_interface(PPI_object,cluster_interface(PPI_object))
    sys.stdout.close()
