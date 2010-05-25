from elementtree import ElementTree
from numpy import array

# Parse an xml file from PDB database
tree = ElementTree.parse("1fpv.xml")

# Get the root of the tree for experimental purpouses
root = tree.getroot()

# Get all the atom_site nodes
namespace="http://pdbml.pdb.org/schema/pdbx-v32.xsd"
searchString = ".//{%s}atom_site" % namespace
atoms = tree.findall(searchString)

# Using list comprehension build a list from atoms:
# each element of the list will be another list,
# representing a 3-dimensional vector that contains
# the XYZ coordinates for the nodes (atoms) of type N or CA
NCA_atoms_xyz = [[float(node.find('.//{%s}Cartn_x' % namespace).text),
                  float(node.find('.//{%s}Cartn_y' % namespace).text),
                  float(node.find('.//{%s}Cartn_z' % namespace).text)]
                 for node in atoms
                 if node.find('.//{%s}label_atom_id' % namespace).text == 'N'
                 or node.find('.//{%s}label_atom_id' % namespace).text == 'C'
                 or node.find('.//{%s}label_atom_id' % namespace).text == 'CA']

# print NCA_atoms_xyz

# Calculate the center of mass:
# NB: seems like every atom weights 1 (???),
# so just sum up all the vectors component-wisely
# and then divide by the number of atoms.
mass_center = [0.0, 0.0, 0.0]
num_atoms = 0
for atom_xyz in NCA_atoms_xyz:
    num_atoms += 1
    mass_center[0] += atom_xyz[0]
    mass_center[1] += atom_xyz[1]
    mass_center[2] += atom_xyz[2]

mass_center[0] /= num_atoms
mass_center[1] /= num_atoms
mass_center[2] /= num_atoms

#print mass_center

atoms_xyz = [(float(node.find('.//{%s}Cartn_x' % namespace).text),
              float(node.find('.//{%s}Cartn_y' % namespace).text),
              float(node.find('.//{%s}Cartn_z' % namespace).text),
              node.find('.//{%s}type_symbol' % namespace).text
              )
             for node in atoms
             ]

#print atoms_xyz

atom_weights = { 'C' : 12.0107,
                 'N' : 14.0067, 
                 'O' : 15.99994,
                 'S' : 32.065
                 }

tensors = [[x*x, x*y, x*z, x, y*y, y*z, y, z*z, z, atom_weights[m]] for (x, y, z, m) in atoms_xyz]

#print lista kaoz84@gmail.com

sum_tensor = array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

for tensor in tensors:
    sum_tensor += array(tensor)

print sum_tensor
