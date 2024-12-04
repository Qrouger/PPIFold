""" Create File_proteins object

    Author: Quentin Rouger
"""
import urllib.request
import re
from utils import *
import csv
import os
import copy

class File_proteins() :
    """
    Manipulate and save file who contains all proteins.
    """
    def __init__ (self, path_txt_file) :
        """
        Constructor : 
        Set attribute for one entry file.

        Parameters:
    	-----------
        path_txt_file : string
        """
        self.set_all_att(path_txt_file)

    def set_proteins_sequence (self, new_protein_sequence) :
        """
        Sets a dict of all sequences.
        
        Parameters:
        ----------
        new_protein_sequence = dictionary
        
        Returns:
        ----------
        """
        self.protein_sequence = new_protein_sequence

    def set_proteins (self, new_protein) :
        """
        Sets a list of all proteins UniprotID.
        
        Parameters:
        ----------
        new_protein = list
        
        Returns:
        ----------
        """
        self.protein = new_protein

    def set_file_name (self, filename) :
        """
        Sets new filename for the txt file.
        
        Parameters:
        ----------
        filename = string
        
        Returns:
        ----------
        """
        self.file_name = filename

    def set_fasta_file (self, filename) :
        """
        Sets new filename for the fasta file.
        
        Parameters:
        ----------
        filename = string
        
        Returns:
        ----------
        """
        self.fasta_file = filename

    def set_lenght_prot (self, lenght_prot) :
        """
        Sets lenght of all proteins.
        
        Parameters:
        ----------
        lenght_prot = dictionary
        
        Returns:
        ----------
        """
        self.lenght_prot = lenght_prot

    def set_names (self, name) :
        """
        Sets names of all proteins.
        
        Parameters:
        ----------
        name = dictionary
        
        Returns:
        ----------
        """
        self.name = name

    def set_iQ_score_dict (self, iQ_score_dict) : 
        """
        Sets iQ_score for all proteins.
        
        Parameters:
        ----------
        iQ_score_dict = dictionary
        
        Returns:
        ----------
        """
        self.iQ_score_dict = iQ_score_dict

    def set_hiQ_score_dict (self, hiQ_score_dict) :
        """
        Sets hiQ_score for all proteins.
        
        Parameters:
        ----------
        hiQ_score_dict = dictionary
        
        Returns:
        ----------
        """
        self.hiQ_score_dict = hiQ_score_dict

    def set_interface_dict (self, interface_dict) :
        """
        Sets interface_dict for new interactions.
        
        Parameters:
        ----------
        interface_dict = dictionary
        
        Returns:
        ----------
        """
        self.interface_dict = interface_dict

    def set_new_pickle (self, new_pickle) :
        """
        Sets list of new pickle feature.
        
        Parameters:
        ----------
        new_pickle = list
        
        Returns:
        ----------
        """
        self.new_pickle = new_pickle
        
    def get_proteins_sequence (self) :
        """
        Return the new amino acid sequence list.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        proteins_sequence : dictionary
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
    
    def get_file_name (self) :
        """
        Return the name of the file.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        file_name : string
        """
        return self.file_name
    
    def get_fasta_file (self) :
        """
        Return the name of the fasta file.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        fasta_file : string
        """
        return self.fasta_file
    
    def get_lenght_prot (self) :
        """
        Return the lenght of proteins.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        lenght_prot : dictionary
        """
        return self.lenght_prot
    
    def get_names (self) :
        """
        Return names of proteins.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        name : dictionary
        """
        return self.name

    def get_iQ_score_dict (self) :
        """
        Return iQ_score for all interactions.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        iQ_score_dict : dictionary
        """
        return self.iQ_score_dict

    def get_hiQ_score_dict (self) :
        """
        Return hiQ_score for all homo-oligomer.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        hiQ_score_dict : dictionary
        """
        return self.hiQ_score_dict
        
    def get_interface_dict (self) :
        """
        Return all interfaces for all proteins.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        interface_dict : dictionary
        """
        return self.interface_dict
    
    def get_new_pickle (self) :
        """
        Return new UniprotID list who did not have pickle feature.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        new_pickle : list
        """
        return self.new_pickle
    
### Generating of features and pre-file to run multimer

    def set_all_att (self, path_txt) :
        """
        Set all values for all attribut for one txt file.
        
        Parameters:
        ----------
        filename : string
        
        Returns:
        ----------
        """
        new_proteins = list()
        interface_dict = dict()
        with open(path_txt,"r") as in_file :
            for line in in_file :
                new_line = (line.split(","))
                for prot in new_line :
                    new_proteins.append(prot.upper().strip())
        self.set_file_name(path_txt)
        self.set_proteins(new_proteins)
        self.set_interface_dict(interface_dict)
 
    def find_proteins_sequence (self) :
        """
        Search in the site uniprot the amino acid sequence and clean it.
        
        Parameters:
        ----------

        Returns:
        ----------
        """
        sequences = dict()
        names = dict()
        pattern = r"SQ   SEQUENCE   .*  .*\n([\s\S]*)"
        pattern2 = r"GN   Name=([\w]*)"
        del_car = ["\n"," ","//"]
        for proteins in self.get_proteins() :
            print("Search sequence for " + proteins)
            urllib.request.urlretrieve("https://rest.uniprot.org/uniprotkb/"+proteins+".txt","temp_file.txt")
            with open("temp_file.txt","r") as in_file:
                for seq in re.finditer(pattern, in_file.read()):
                    sequences[proteins] = seq.group(1)
            with open("temp_file.txt","r") as in_file:
                for name in re.finditer(pattern2, in_file.read()) :
                    names[proteins] = name.group(1)
            for car in del_car :
                sequences[proteins] = sequences[proteins].replace(car,"")
        os.remove("temp_file.txt")
        self.set_proteins_sequence(sequences)
        self.set_names(names)

    def find_prot_lenght (self, prot_dict = None) :
        """
        Set the lenght for all proteins.

        Parameters:
        ----------
        prot_dict = dictionary
        
        Returns:
        ----------
        """
        if prot_dict == None :
            proteins = self.get_proteins()
        else :
            proteins = prot_dict
        sequences = self.get_proteins_sequence()
        lenght_prot = dict()
        for protein in proteins :
            lenght_prot[protein] = len(sequences[protein])
        self.set_lenght_prot(lenght_prot)

    def create_fasta_file (self) :
        """
        Generate a fasta file with the initial txt file.

        Parameters:
        ----------

        Returns:
        ----------
        """
        line = str()
        proteins = self.get_new_pickle()
        sequences = self.get_proteins_sequence()
        for protein in proteins :
            line = line + ">" + protein + "\n" + sequences[protein] + "\n"        
        file_name = self.get_file_name()
        file_out = file_name.replace("txt","fasta")
        with open(file_out,"w") as fh :
            fh.write(line)
        self.set_fasta_file(file_out)

    def update_iQ_hiQ_score (self) :
        """
        Generate two dictionaries, first where the key is a tuple of interaction proteins(Uniprot) and the value is the iQ_score, a second where the key is the protein (Uniprot) and the value is a tuple of the better hiQ_score and this homo-oligomerisation.

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

    def already_pickle (self, pickle_path) :
        """
        Check if protein have already a feature pickle file, and return a list of proteins who doesn't.

        Parameters:
        ----------
        pickle_path : string

        Returns:
        ----------
        """
        prot_need_pkl = list()
        proteins = self.get_proteins()
        for uniprotID in proteins :
            if os.path.isfile(pickle_path + "/" + uniprotID + ".pkl") :
                pass
            else :
                prot_need_pkl.append(uniprotID)
        self.set_new_pickle(prot_need_pkl)

    def define_interface (self, list_of_list_int, int) :
        """
        Set a dictionnary with all residues in interaction, with the UniprotID.

        Parameters:
        ----------
        list_of_list_int : list
        int : string

        Returns:
        ----------
        """
        all_residues_int = copy.deepcopy(list_of_list_int)
        old_interface_dict = self.get_interface_dict()
        protein1 = int[0]
        protein2 = int[1]
        list_int_protein1 = list()
        list_int_protein2 = list()
        if protein1 not in old_interface_dict.keys() :
            old_interface_dict[protein1] = []
        if protein2 not in old_interface_dict.keys() :
            old_interface_dict[protein2] = []
        for line in all_residues_int :
            line[1] = line[1].strip() #remove readability spaces
            if line[0] != int[0] :
                if line[0] not in list_int_protein1 :
                    list_int_protein1.append(line[0])
                if line[1] not in list_int_protein2 :
                    list_int_protein2.append(line[1])
        list_int_protein1.append(protein2) #last values of each list is the second proteins
        list_int_protein2.append(protein1)
        old_interface_dict[protein1].append(list_int_protein1)
        old_interface_dict[protein2].append(list_int_protein2)
        self.set_interface_dict(old_interface_dict)
