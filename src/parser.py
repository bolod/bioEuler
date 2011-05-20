import os
import sys

from Bio.PDB.PDBParser import PDBParser
from numpy import *

def get_structure(name):
	parser = PDBParser()
	structure = parser.get_structure(name, pdb_path + name + '.pdb')
	return structure

def get_structure_backbone(structure):
	backbone = []
	for model in structure.get_list():
		for chain in model.get_list():
			for residue in chain.get_list(): 
				if residue.has_id('CA'):
					atom = residue['CA']
					backbone.append(atom)
	return backbone

def get_backbone_centroid(backbone):
	n = len(backbone)
	coord = [0,0,0]
	for atom in backbone:
		coord += atom.get_coord()
	coord /= n
	return coord

def get_atom_euler(atom):
	coords = atom.get_coord()
	euler = [ coord * coords for coord in coords ]
	return euler

def get_backbone_euler(backbone):
	euler = sum([get_atom_euler(atom) for atom in backbone], axis=0)
	return euler

def get_backbone_axis(backbone):
	euler = get_backbone_euler(backbone)
	eigen_vec = transpose(linalg.eig(euler)[1])
	eigen_val = linalg.eigvals(euler)
	return [ eigen_vec[i] for i in eigen_val.argsort() ]

def translate_backbone(backbone):
	centroid = get_backbone_centroid(backbone)
	for atom in backbone:
		atom.coord -= centroid
	return backbone

def orient_backbone(backbone):
	translation = -1 * get_backbone_centroid(backbone)
	rotation = get_backbone_axis(backbone)
	for atom in backbone:
		atom.transform(rotation, translation)
	return backbone

def print_backbone(backbone):
	for atom in backbone:
		print atom.get_coord()
	print get_backbone_centroid(backbone)

def split_backbone(backbone):
	list = [[]]
	n = 0
	previous = backbone[0]
	for atom in backbone:
		if previous.get_coord()[2] * atom.get_coord()[2] < 0:
			list.append([])
			n += 1
		list[n].append(atom)
		previous = atom
	return list



if __name__ == "__main__":

	current_path = os.path.dirname(sys.argv[0])
	pdb_path = current_path + '/../pdb/'
	pdb_file_name = '1alc'

	structure = get_structure(pdb_file_name)
	backbone = get_structure_backbone(structure)
	euler = get_backbone_euler(backbone)

	backbone = translate_backbone(backbone)
	backbone = orient_backbone(backbone)

	print '1alc oriented backbone'
	print_backbone(backbone)
	print '\n'

	print '1alc splitted backbone'
	backbone_list = split_backbone(backbone)
	for backbone in backbone_list:
		print '***'
		print_backbone(backbone)