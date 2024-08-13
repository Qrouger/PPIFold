from utils import *
from File_proteins import *
import argparse


def add_arguments(parser) :
    parser.add_argument("--txt_name", help = "Txt file name",required = True)
    parser.add_argument("--data_dir", help = "PATH to alphafold data",required = True)
    parser.add_argument("--dir_alpha_analysis" , help = "Directory of the singularity image", required =  True)
    parser.add_argument("--make_multimers", help = "If you just want make feature set False", required = False, default = True)
    parser.add_argument("--env_feature" , help = "Conda environemment to make feature", required = False, default = None)
    parser.add_argument("--env_multimer" , help = "Conda environemment to make multimers", required = False, default = None)
    parser.add_argument("--max_aa" , help = "Maximum amino acids can be generate by your cluster", required = False, default = 2500, type = int)
    parser.add_argument("--use_signalP" , help = "Don't use SignalP", required = False, default = True)

    def __init__(self, args) :
        super().__init__(args)

def main(A4, args) :
    if args.use_signalP == True :
        A4.remove_SP()
    A4.create_feature(args.env_feature,args.data_dir)
    A4.Make_all_MSA_coverage()
    A4.generate_APD_script(args.max_aa)
    if args.make_multimers == True :
        Make_all_vs_all(args.env_multimer,args.data_dir)
        add_iQ_score(args.dir_alpha_analysis)
        create_out_fig()
        Make_homo_oligo(args.env_multimer,args.data_dir)
        add_hiQ_score(args.dir_alpha_analysis)
        A4.generate_interaction_network()

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args()
    self.find_proteins_sequence()
    self.find_prot_lenght()
    self.create_fasta_file()
    A4 = EZFold(args)
    main(A4, args)
