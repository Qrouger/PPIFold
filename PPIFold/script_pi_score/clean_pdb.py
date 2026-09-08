#Adapted from clean_pdb.py (https://gitlab.com/topf-lab/pi_score/-/tree/master/score_scripts)

from Bio.PDB import *

dict_aa = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
           'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
           'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
           'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def str_object(pdbfile):
    '''
    Given a pdb file creates a structure object
    '''
    parser = PDBParser(QUIET=True)
    s = parser.get_structure('id', pdbfile)
    return s

class SelectResidues(Select):
    """ Only accept the specified residues when saving. """
    def accept_residue(self, residue):
        #print residue.get_resname(),residue.get_full_id()[3][0]
        if residue.get_full_id()[3][0] != 'W' and not residue.get_full_id()[3][0].startswith('H_'):
            if residue.get_resname() in dict_aa.keys():
                return 1

    def accept_atom(self,atom):
        if not atom.is_disordered() or atom.get_altloc() == 'A':
                if atom.element.strip() != 'H':
                        return 1

def change_atloc(pdbfile,
                out=None):
    if out == None:
        out = (pdbfile.rsplit('.')[0]
                + '_atm_loc'
                + '.pdb')
    
    outfle = open(out,'w')
    with open(pdbfile) as pdb:
        for lne in pdb:
            try:
                if len(lne.split()[3]) >3:
                    lne = lne.replace(lne.split()[3],' ' + lne.split()[3][1:])
                elif len(lne.split()[2]) > 3:
                    for keys in dict_aa.keys():
                        if keys in lne.split()[2]:
                            lne = lne.replace(lne.split()[2], lne.split()[2][0:3] + ' ' + lne.split()[2][4:])
            except: pass
            outfle.write(lne)
    outfle.close()
    return out

def write_sel_pdb(pdbfile):
    writer = PDBIO()
    struct = str_object(pdbfile)
    writer.set_structure(struct)
    outfle = (pdbfile.rsplit('.',1)[0]
                + '_clean'
                + '.pdb')
    writer.save(outfle, select=SelectResidues())

    return outfle
