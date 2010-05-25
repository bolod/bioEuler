import Bio.PDB.PDBParser as pdbp
from numpy import *

pdb_path = 'C:\\Users\\ivanagloriosi\\Documents\\My Dropbox\\UNI\\BioInformatica\\Marino-Spini\\1fpv.pdb'

parser = pdbp()
structure = parser.get_structure("test", pdb_path)
model = structure[0]

acc = 0

for chain in model.get_list():
    for residue in chain.get_list():
        for atom in residue.get_list():
            print atom.get_name() + " coord = " + str(atom.get_coord()) 
            acc += 1

print acc
