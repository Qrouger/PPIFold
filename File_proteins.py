import os
import numpy as np
import urllib.request
import re
import pickle
import csv
import matplotlib.pyplot as plt
from Bio import PDB
import numpy as np
from tabulate import tabulate
import pandas as pd
import copy
import glob
import pickle
import logging
from scipy.special import softmax



class File_proteins() :
    """
    Manipulate all of Proteins.
    """
    def __init__(self, args) :
        """
        Constructor :method for Interactome:

        Parameters:
	-----------
        filename : string
        file of proteins of interest
        """
        self.set_all_att(args.txt_name)

    def set_proteins_sequence(self, new_protein_sequence) :
        """
        Sets a list of all sequence. 
        Parameters:
        ----------
        new_sequence = list
        
        Returns:
        ----------
        """
        self.protein_sequence = new_protein_sequence

    def set_proteins(self, new_protein) :
        """
        Sets a list of all proteins names. 
        Parameters:
        ----------
        new_sequence = list
        
        Returns:
        ----------
        """
        self.protein = new_protein

    def set_file_name(self, filename) :
        """
        Sets new filename for the use file. 
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
mpDockQ/pDockQ        filename = string
        
        Returns:
        ----------
        """
        self.fasta_file = filename

    def set_lenght_prot(self, lenght_prot) :
        """
        Sets lenght of all proteins.
        Parameters:
        ----------
        filename = string
        
        Returns:
        ----------
        """
        self.lenght_prot = lenght_prot

    def get_proteins_sequence(self) :
        """
        Return the new aa sequence list. 
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
        proteins_sequence : list
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
        pattern =r"SQ   SEQUENCE   .*  .*\n([\s\S]*)"
        print(self.get_proteins())
        for proteins in self.get_proteins() :
            print("Search sequence for " + proteins)
            urllib.request.urlretrieve("https://rest.uniprot.org/uniprotkb/"+proteins+".txt","temp_file.txt")
            with open("temp_file.txt","r") as in_file:
                for m in re.finditer(pattern, in_file.read()):
                    sequences.append(m.group(1))
        for sequence in sequences :
            del_car = ["\n"," ","//"]
            for car in del_car :
                sequence = sequence.replace(car,"")
            new_sequences.append(sequence)
        self.set_proteins_sequence(new_sequences)

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
    
