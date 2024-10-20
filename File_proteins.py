""" Create File_proteins object

    Author: Quentin Rouger
"""
import urllib.request
import re
from utils import *
import csv

class File_proteins() :
    """
    Manipulate and save file who contains all proteins
    """
    def __init__(self, args) :
        """
        Constructor : 
        Set attribute for one entry file.

        Parameters:
    	-----------
        args : ?
        """
        self.set_all_att(args.txt_name)

    def set_proteins_sequence(self, new_protein_sequence) :
        """
        Sets a list of all sequences.
        
        Parameters:
        ----------
        new_protein_sequence = list
        
        Returns:
        ----------
        """
        self.protein_sequence = new_protein_sequence

    def set_proteins(self, new_protein) :
        """
        Sets a list of all proteins names.
        
        Parameters:
        ----------
        new_protein = list
        
        Returns:
        ----------
        """
        self.protein = new_protein

    def set_file_name(self, filename) :
        """
        Sets new filename for the txt file.
        
        Parameters:
        ----------
        filename = string
        
        Returns:
        ----------
        """
        self.file_name = filename

    def set_fasta_file(self, filename) :
        """
        Sets new filename for the fasta file.
        
        Parameters:
        ----------
        filename = string
        
        Returns:
        ----------
        """
        self.fasta_file = filename

    def set_lenght_prot(self, lenght_prot) :
        """
        Sets lenght of all proteins.
        
        Parameters:
        ----------
        lenght_prot = dict
        
        Returns:
        ----------
        """
        self.lenght_prot = lenght_prot

    def set_names(self, name) :
        """
        Sets names of all proteins.
        
        Parameters:
        ----------
        name = dict
        
        Returns:
        ----------
        """
        self.name = name

    def set_iQ_score_dict(self, iQ_score_dict) :
        self.iQ_score_dict = iQ_score_dict

    def set_hiQ_score_dict(self, hiQ_score_dict) :
        self.hiQ_score_dict = hiQ_score_dict

    def get_proteins_sequence(self) :
        """
        Return the new amino acid sequence list.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        proteins_sequence : list
        """
        return self.protein_sequence
    
    def get_proteins (self) :
        """
        Return the new proteins name list.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        protein : list
        """
        return self.protein
    
    def get_file_name(self) :
        """
        Return the name of the file.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        file_name : str
        """
        return self.file_name
    
    def get_fasta_file(self) :
        """
        Return the name of the fasta file.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        fasta_file : str
        """
        return self.fasta_file
    
    def get_lenght_prot(self) :
        """
        Return the lenght of proteins.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        lenght_prot : dict
        """
        return self.lenght_prot
    
    def get_names(self) :
        """
        Return names of proteins.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        name : dict
        """
        return self.name

    def get_iQ_score_dict(self) :
        """
        Return iQ_score for all interactions.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        iQ_score_dict : dict
        """
        return self.iQ_score_dict

    def get_hiQ_score_dict(self) :
        """
        Return hiQ_score for all homo-oligomer.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        hiQ_score_dict : dict
        """
        return self.hiQ_score_dict
    
### Generating of features and pre-file to run multimer

    def set_all_att(self, filename) :
        """
        Set all values for all attribut for one txt file.
        
        Parameters:
        ----------
        filename : string
        
        Returns:
        ----------
        """
        new_proteins = []
        with open(filename,"r") as in_file :
            for line in in_file :
                new_line = (line.split(","))
                for prot in new_line :
                    new_proteins.append(prot.upper().strip())
        self.set_file_name(filename)
        self.set_proteins(new_proteins)
 
    def find_proteins_sequence (self) :
        """
        Search in the site uniprot the aa sequence and clean it.
        
        Parameters:
        ----------

        Returns:
        ----------
    
        """
        new_sequences = list()
        sequences = list()
        names = dict()
        pattern = r"SQ   SEQUENCE   .*  .*\n([\s\S]*)"
        pattern2 = r"GN   Name=([\w]*)"
        print(self.get_proteins())
        for proteins in self.get_proteins() :
            print("Search sequence for " + proteins)
            urllib.request.urlretrieve("https://rest.uniprot.org/uniprotkb/"+proteins+".txt","temp_file.txt")
            with open("temp_file.txt","r") as in_file:
                for seq in re.finditer(pattern, in_file.read()):
                    sequences.append(seq.group(1))
            with open("temp_file.txt","r") as in_file:
                for name in re.finditer(pattern2, in_file.read()) :
                    names[proteins] = name.group(1)
        print (names)
        for sequence in sequences :
            del_car = ["\n"," ","//"]
            for car in del_car :
                sequence = sequence.replace(car,"")
            new_sequences.append(sequence)
        self.set_proteins_sequence(new_sequences)
        self.set_names(names)

    def find_prot_lenght(self) :
        """
        Set the lenght for all proteins.

        Parameters:
        ----------

        Returns:
        ----------

        """
        proteins = self.get_proteins()
        sequences = self.get_proteins_sequence()
        lenght_prot = dict()
        for nbr_prot in range(len(proteins)) :
            lenght_prot[proteins[nbr_prot]] = len(sequences[nbr_prot])
        self.set_lenght_prot(lenght_prot)

    def create_fasta_file (self) :
        """
        Generate a fasta file with a txt file.

        Parameters:
        ----------

        Returns:
        ----------

        """
        line = str()
        proteins = self.get_proteins()
        sequences = self.get_proteins_sequence()
        for protein in range(0,len(proteins)) :
            line = line + ">" + proteins[protein] + "\n" + sequences[protein] + "\n"        
        file_name = self.get_file_name()
        file_out = file_name.replace("txt","fasta")
        with open(file_out,"w") as fh :
            fh.write(line)
        self.set_fasta_file(file_out)

    def update_iQ_hiQ_score(self) :
        """
        Generate two dictionary, first where the key is a tuple of interaction proteins(Uniprot) and the value is the iQ_score, a second where the key is the protein (Uniprot) and the value is a tuple of the better hiQ_score and this homo-oligomerisation.

        Parameters:
        ----------

        Returns:
        ----------

        """
        iQ_score_dic = dict()
        with open("result_all_vs_all/predictions_with_good_interpae.csv", "r") as file1 :
            reader1 = csv.DictReader(file1)
            for row in reader1 :
                names = row['jobs'].split('_and_')
                iQ_score_dic[(names[0],names[1])] = row['iQ_score']
        self.set_iQ_score_dict(iQ_score_dic)
        hiQ_score_dic = dict()
        with open("result_homo_oligo/predictions_with_good_interpae.csv", "r") as file2 :
            reader2 = csv.DictReader(file2)
            for row in reader2 :
                prot_name = row['jobs'].split("_homo_")[0]
                if prot_name not in hiQ_score_dic.keys() or float(row['hiQ_score']) >= hiQ_score_dic[prot_name][0] :
                    number_homo = int((row['jobs'].split("homo_")[1]).split("er")[0]) #to take the number of homo-oligomerisation of the protein and this score
                    hiQ_score_dic[prot_name] = (float(row['hiQ_score']),number_homo)
        self.set_hiQ_score_dict(hiQ_score_dic)
