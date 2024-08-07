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

class EZFold (File_proteins) :

    def __init__(self) :
        parser = argparse.ArgumentParser()
        add_arguments(parser)
        self.args = parser.parse_args()
        super().__init__(self.args)

    def main(self) :
        parser = argparse.ArgumentParser()
        add_arguments(parser)
        self.args = parser.parse_args()
        super().__init__(self.args)
        super().find_proteins_sequence()
        super().find_prot_lenght()
        super().create_fasta_file()
        if self.args.use_signalP == True :
           remove_SP(self)
        #create_feature(self.args.env_feature,self.args.data_dir,self)
        #Make_all_MSA_coverage(self)
        generate_APD_script(self.args.max_aa, self)
        if self.args.make_multimers == True :
            #Make_all_vs_all(self.args.env_multimer,self.args.data_dir)
            add_iQ_score(self.args.dir_alpha_analysis)
            create_out_fig()
            #Make_homo_oligo(self.args.env_multimer,self.args.data_dir)
            add_hiQ_score(self.args.dir_alpha_analysis)
            generate_interaction_network(self)

if __name__ == "__main__" :
    A4 = EZFold()
    A4.main()
