from elementtree import ElementTree
from numpy import array
from pyplasm import *

atom_weights = { 'C' : 12.0107,
                 'N' : 14.0067, 
                 'O' : 15.99994,
                 'S' : 32.065
                 }

# Parse an xml file from PDB database
tree = ElementTree.parse("1fpv.xml")

# Get the root of the tree for experimental purpouses
root = tree.getroot()

# Get all the atom_site nodes
namespace="http://pdbml.pdb.org/schema/pdbx-v32.xsd"
searchString = ".//{%s}atom_site" % namespace
atoms = tree.findall(searchString)

atoms_xyzm = [(float(node.find('.//{%s}Cartn_x' % namespace).text),
              float(node.find('.//{%s}Cartn_y' % namespace).text),
              float(node.find('.//{%s}Cartn_z' % namespace).text),
              atom_weights[node.find('.//{%s}type_symbol' % namespace).text]
              )
             for node in atoms
             ]

atom_list = [(x, y, z) for (x, y, z, m) in atoms_xyzm]

atom_per_m_list = [(x * m, y * m, z * m) for (x, y, z, m) in atoms_xyzm]

sum_m = sum([m for (x, y, z, m) in atoms_xyzm])

bar_x = sum([xi_per_mi for (xi_per_mi, yi_per_mi, zi_per_mi) in atom_per_m_list]) / sum_m
bar_y = sum([yi_per_mi for (xi_per_mi, yi_per_mi, zi_per_mi) in atom_per_m_list]) / sum_m
bar_z = sum([zi_per_mi for (xi_per_mi, yi_per_mi, zi_per_mi) in atom_per_m_list]) / sum_m

bar = (bar_x, bar_y, bar_z)

print bar

atom_trans_list = [(x - bar_x, y - bar_y, z - bar_z) for (x, y, z, m) in atoms_xyzm]

def point(x,y,z):
	s = CUBOID([.5,.5,.5])
	s = T([1,2,3])([x,y,z])(s)
	return s
	
point_list = []

point_trans_list = []

for (x, y, z) in atom_trans_list:
	point_trans_list.append(point(x, y, z))
	
for (x, y, z) in atom_list:
	point_list.append(COLOR(RED)(point(x,y,z)))

big_point = COLOR(GREEN)(CUBOID([5,5,5]))

bar_point = T([1,2,3])([bar_x, bar_y, bar_z])(big_point)

VIEW(STRUCT(point_list + point_trans_list + [bar_point]))
