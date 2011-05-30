import os
import sys

from Bio.PDB import PDBParser, PPBuilder, Polypeptide
from numpy import array, matrix, linalg

def get_structure(name, path):
    parser = PDBParser()
    if not(path.endswith('/')):
        path += '/'
    structure = parser.get_structure(name, path + name + '.pdb')
    return structure


def first_split(pp_list):
    """
    Return a list of protein's segment derived by first splitting stage.

    @param
        - pp_list: list of Polypeptide_s
    """
    ca_splitted_segments = []

    for pp in pp_list:
            temp_seg = []
            for i in range(len(pp)-1):
                temp_seg.append(pp[i])
                caA = pp[i]['CA']
                caB = pp[i+1]['CA']
                if caA.get_coord()[2] * caB.get_coord()[2] < 0:
                    ca_splitted_segments.append(temp_seg)
                    temp_seg = []
            ca_splitted_segments.append(temp_seg)
    return ca_splitted_segments

def orient(pp_list):
    all_atoms = [a for pp in pp_list for res in pp for a in res]
    centroid = get_centroid(all_atoms)

    for a in all_atoms:
        a.transform([[1,0,0],[0,1,0],[0,0,1]], -1 * centroid)


    euler = reduce(lambda x,y: x+y, map(get_atom_euler, all_atoms))

    print euler

    rot = get_principal_axis(euler)

    print rot
    
    for a in all_atoms:
        a.transform(rot, [0,0,0])

def get_centroid(atom_list):
    n = len(atom_list)
    coord = array([0,0,0])
    for a in atom_list:
      coord += a.get_coord()
    coord /= n
    return coord    


def get_atom_euler(atom):
    coords = atom.get_coord()
    mass = atom.mass
    a_euler =  array([ coord * (mass * coords) for coord in coords ])
    return a_euler

def get_principal_axis(euler):
    eigen_vec = (linalg.eig(euler)[1]).T
    eigen_val = linalg.eigvals(euler)
    rot = [ eigen_vec[i] for i in eigen_val.argsort() ]

#   if det(self.rot)<0:
#       vt[2]=-vt[2]
#       self.rot=transpose(dot(transpose(vt), transpose(u)))
#   self.tran=av2-dot(av1, self.rot)

    return rot

# def get_structure_backbone(structure):
#   backbone = []
#   for model in structure.get_list():
#       for chain in model.get_list():
#           for residue in chain.get_list(): 
#               if residue.has_id('CA'):
#                   atom = residue['CA']
#                   backbone.append(atom)
#   return backbone




# def get_backbone_centroid(backbone):
#   n = len(backbone)
#   coord = [0,0,0]
#   for atom in backbone:
#       coord += atom.get_coord()
#   coord /= n
#   return coord

# def get_atom_euler(atom):
#   coords = atom.get_coord()
#   euler = [ coord * coords for coord in coords ]
#   return euler

# def get_backbone_euler(backbone):
#   euler = sum([get_atom_euler(atom) for atom in backbone], axis=0)
#   return euler

# def get_backbone_axis(backbone):
#   euler = get_backbone_euler(backbone)
#   eigen_vec = transpose(linalg.eig(euler)[1])
#   eigen_val = linalg.eigvals(euler)
#   rot = [ eigen_vec[i] for i in eigen_val.argsort() ]
    



# # if det(self.rot)<0:
# #             vt[2]=-vt[2]
# #             self.rot=transpose(dot(transpose(vt), transpose(u)))
# # self.tran=av2-dot(av1, self.rot)    

# def translate_backbone(backbone):
#   centroid = get_backbone_centroid(backbone)
#   for atom in backbone:
#       atom.coord -= centroid
#   return backbone

# def orient_backbone(backbone):
#   translation = -1 * get_backbone_centroid(backbone)
#   rotation = get_backbone_axis(backbone)
#   for atom in backbone:
#       atom.transform(rotation, translation)
#   return backbone

# def print_backbone(backbone):
#   for atom in backbone:
#       print atom.get_coord()
#   print get_backbone_centroid(backbone)

# def split_backbone(backbone):
#   list = [[]]
#   n = 0
#   previous = backbone[0]
#   for atom in backbone:
#       if previous.get_coord()[2] * atom.get_coord()[2] < 0:
#           list.append([])
#           n += 1
#       list[n].append(atom)
#       previous = atom
#   return list

if __name__ == "__main__":

    current_path = os.path.dirname(sys.argv[0])
    pdb_path = current_path + '../pdb/'
    pdb_id = '2vb1'


    structure = get_structure(pdb_id, pdb_path)
    model = structure[0]

    ppb = PPBuilder()
    pp_list = ppb.build_peptides(model)

    # orient
    orient(pp_list)

    # first split stage
    fs = first_split(pp_list)





    for seg in fs:
        pp = Polypeptide.Polypeptide(seg)
        print pp.get_sequence()










    #--------------------------------------------
    # backbone = get_structure_backbone(structure)
    # euler = get_backbone_euler(backbone)

    # backbone = translate_backbone(backbone)
    # backbone = orient_backbone(backbone)

    # print '1alc oriented backbone'
    # print_backbone(backbone)
    # print '\n'

    # print '1alc splitted backbone'
    # backbone_list = split_backbone(backbone)
    # for backbone in backbone_list:
    #     print '***'
    #     print_backbone(backbone)