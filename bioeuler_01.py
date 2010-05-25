import xml.etree.ElementTree as ET
from numpy import *
from pyplasm import *

# Parse an xml file from PDB database
tree = ET.parse("1alc.xml")

# Get the root of the tree for experimental purpouses
root = tree.getroot()

# Get all the atom_site nodes
namespace="http://pdbml.pdb.org/schema/pdbx-v32.xsd"
searchString = ".//{%s}atom_site" % namespace
atoms = tree.findall(searchString)

atom_weights = { 'C' : 12.0107,
                 'CA': 12.0107,
                 'N' : 14.0067, 
                 'O' : 15.99994,
                 'S' : 32.065 }

atoms_xyzm = array([ [float(i.find(".//{%s}Cartn_x" % namespace).text),
                     float(i.find(".//{%s}Cartn_y" % namespace).text),
                     float(i.find(".//{%s}Cartn_z" % namespace).text),
                     atom_weights[
                     i.find(".//{%s}type_symbol" % namespace).text]]
                    for i in atoms ])

# Centro di massa = sum(m*r) / sum(m)
sum_mr = sum ( [  i[:-1] * i[-1]  
                for i in atoms_xyzm ], axis=0 )

sum_m = sum ( [ i[-1] 
                for i in atoms_xyzm ], axis=0 )

mass_center = sum_mr / sum_m
# fine calcolo centro di massa

# Baricentro = sum(r) / n
centroid = (sum ([ i[:-1] for i in atoms_xyzm ], axis=0)) / atoms_xyzm.shape[0]

def cubo(xyz):
    	s = T([1,2,3])(xyz)(CUBOID([.5,.5,.5]))
    	return s
    
nuvola = [cubo(x[:-1]) for x in atoms_xyzm ]
centro_massa_nuvola = [COLOR(RED)(cubo(mass_center))]
centroid_nuvola = [COLOR(GREEN)(cubo(centroid))]

VIEW ( STRUCT(nuvola + centro_massa_nuvola + centroid_nuvola) )
