import argparse

from utils import *
from File_proteins import *
from alphapulldown.alphapulldown.scripts import create_individual_features
import sys
sys.path.append("alphapulldown")

def add_arguments(parser) :
    parser.add_argument("--make_multimers", help = "If you just want make feature set on False", required = False, default = True)
    parser.add_argument("--env_feature" , help = "Conda environemment name to make feature", required = False, default = None)
    parser.add_argument("--env_multimer" , help = "Conda environemment name to make multimers", required = False, default = None)
    parser.add_argument("--max_aa" , help = "Maximum amino acids can be generate by your cluster", required = False, default = 2500, type = int)
    parser.add_argument("--use_signalP" , help = "Used or not SignalP", required = False, default = True)
    parser.add_argument("--org" , help = "Organism of interest : arch, gram+, gram- or euk", required = False, default = "gram-", type = str)

def main(A4) :
    #if args.use_signalP == True :
    #    remove_SP(A4,args.org)
    #create_feature(A4,args.env_feature,path_dict["Path_AlphaFold_data"],path_dict["Path_Pickle_Feature"])
    #Make_all_MSA_coverage(A4,path_dict["Path_Pickle_Feature"])
        txt_name = A4.get_fasta_file()
        create_individual_features.main(txt_name,path_dict["Path_AlphaFold_data"],True,path_dict["Path_Pickle_Feature"],"2024-05-02",False)
        generate_APD_script(A4, args.max_aa)
    #if args.make_multimers == True :
        #Make_all_vs_all(args.env_multimer,path_dict["Path_AlphaFold_data"],path_dict["Path_Pickle_Feature"])
        add_iQ_score(path_dict["Path_Singularity_Image"])
    #    Make_homo_oligo(args.env_multimer,path_dict["Path_AlphaFold_data"],path_dict["Path_Pickle_Feature"])
        add_hiQ_score(path_dict["Path_Singularity_Image"])
        A4.update_iQ_hiQ_score()
        generate_heatmap(A4)
        #create_out_fig(A4)
        generate_interaction_network(A4)

if __name__ == "__main__" :
    path_dict = define_path()
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args()
    A4 = File_proteins(path_dict["Path_Uniprot_ID"])
    A4.find_proteins_sequence()
    A4.find_prot_lenght()
    A4.already_pickle(path_dict["Path_Pickle_Feature"])
    A4.create_fasta_file(path_dict["Path_Pickle_Feature"])
    main(A4)


#conf it's a file who contains path for differents data or file