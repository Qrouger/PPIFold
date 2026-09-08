#Adapted from interface_assess.py (https://gitlab.com/topf-lab/pi_score/-/tree/master/score_scripts)

from Bio.PDB.MMCIFParser import MMCIFParser
import sys
from Bio.PDB import Selection, NeighborSearch, PDBParser, Select
import json

dict_aa = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
           'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
           'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
           'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

class SelectChains(Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)

def str_object(pdbfile):
    '''
    Given a pdb file creates a structure object
    '''
    if pdbfile.rsplit('.')[-1] == 'pdb':
        parser = PDBParser(QUIET=True)
    elif pdbfile.rsplit('.')[-1] == 'ent':
        parser = PDBParser(QUIET=True)
    elif pdbfile.rsplit('.')[0] == 'cif':
        parser = MMCIFParser()
    else:
        print 'Not a valid PDB format'
        sys.exit()
    structure = parser.get_structure('id', pdbfile)
    return structure

def get_num_chains(structure):
    '''
    Returns a list of chain objects in a structure
    '''
    ch_list = Selection.unfold_entities(structure, 'C')
    return ch_list

def get_CAatoms(chain_object):
    '''
    for a given chain object, returns a list of CA atoms
    '''
    ca_list = []
    atm_list = Selection.unfold_entities(chain_object, 'A')
    for a in atm_list:
        if a.get_id() == 'CA':
            ca_list.append(a)
    return ca_list
    
def interface_residues(pdbfile,
                    len_chains,
                    numres_cut,
                    dist_cutoff):
    '''
    Given a PDB file calculates interface residues
    '''
#     numres_cut = int(numres_cut)
    structure = str_object(pdbfile)
    ch_list = get_num_chains(structure)
    int_reslist_j = []
    int_reslist_i = []
    dict_intf = {}
    dict_residue_contacts = {}
    for i in xrange(0,len(ch_list)-1):
        for j in xrange(i+1,len(ch_list)):
            int_reslist_j = []
            int_reslist_i = []
            intf_ch = ch_list[j].get_id() + '_' + ch_list[i].get_id()
            ca_i = get_CAatoms(ch_list[i])
            ca_j = get_CAatoms(ch_list[j])
            # print ca_i, intf_ch
            if len(ca_i) <len_chains or len(ca_j) <len_chains: 
                # dict_intf = {} 
                # dict_residue_contacts = {}
                continue
            ns_i = NeighborSearch(ca_i)
            for atms in ca_j:
                n = ns_i.search(atms.get_coord(), dist_cutoff,'R')
                if len(n) >0 :
                    res = (str(atms.get_parent().get_resname()) 
                             + str(atms.get_parent().get_id()[1]))
                    int_reslist_j.append(res)
                    t = []
                    for at in n:
                        res1 = (str(at.get_resname()) 
                             + str(at.get_id()[1]))
                        int_reslist_i.append(res1)
                        t.append(res1)
                    if intf_ch in dict_residue_contacts.keys():  
                        dict_residue_contacts[intf_ch].update({res:t})
                    else:
                        dict_residue_contacts[intf_ch] = {res:t}
            #Number of interface residues from each chain is >=10
            if len(set(int_reslist_j)) >= numres_cut and len(set(int_reslist_i)) >=numres_cut:
                dict_intf[intf_ch] = [list(set(int_reslist_j)),list(set(int_reslist_i))]   
    intf_outfle_prefix = pdbfile.split('.')[0]
    intf_outfle_name = intf_outfle_prefix + '_interface_residues_dict.json'
    
    if dict_intf: 
        with open(intf_outfle_name,'w') as outfle:
            json.dump(dict_intf,outfle)
    return dict_intf, dict_residue_contacts
    
    
def polar_residues(reslist):
    assert len(reslist) >0 
    polar_ref = ['SER', 
                 'THR',
                 'ASN',
                 'GLN',
                 'HIS',
                 'TYR']
    polar_reslist = []
    for r in reslist:
        if r in polar_ref:
            polar_reslist.append(r)
    frac_polar_residue = len(polar_reslist)/float(len(reslist))
    return frac_polar_residue

def hydrophobic_residues(reslist):
    assert len(reslist) >0 
    hydro_ref = ['ALA',
                 'LEU',
                 'ILE',
                 'VAL',
                 'PHE',
                 'TRP',
                 'CYS',
                 'MET'
                ]
    hydro_reslist = []
    for r in reslist:
        if r in hydro_ref:
            hydro_reslist.append(r)
    frac_hydro_residue = len(hydro_reslist)/float(len(reslist))
    return frac_hydro_residue
            
def charged_residues(reslist):
    assert len(reslist) >0 
    charged_ref = ['LYS',
                   'ARG',
                   'ASP',
                   'GLU']
    charged_reslist = []
    for r in reslist:
        if r in charged_ref:
            charged_reslist.append(r)
    frac_charged_residue = len(charged_reslist)/float(len(reslist))
    return frac_charged_residue