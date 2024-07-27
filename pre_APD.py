from utils import *
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

class EZFold (pre_APD):

    def __init__(self) :
        parser = argparse.ArgumentParser()
        add_arguments(parser)
        self.args = parser.parse_args()
        super().__init__(self.args)

    def main(self) :
        super().find_proteins_sequence()
        super().find_prot_lenght()
        #super().create_fasta_file()
        #if args.use_signalP == True :
           #super().remove_SP()
        #super().create_feature(args.env_feature,args.data_dir)
        #super().Make_all_MSA_coverage()
        #super().generate_APD_script(args.max_aa)
        if self.args.make_multimers == True :
            #super().Make_all_vs_all(args.env_multimer,args.data_dir)
            super().add_iQ_score_and_make_cyt(self.args.dir_alpha_analysis)
            super().create_out_fig()
            #super().Make_homo_oligo(args.env_multimer,args.data_dir)
            super().add_hiQ_score(self.args.dir_alpha_analysis)


if __name__ == "__main__" :
    A4 = EZFold()
    A4.main()