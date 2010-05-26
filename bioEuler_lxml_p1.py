from lxml import etree
from numpy import array

masse_atomiche = {'C':12.0107, 'N':14.0067, 'O':15.99994, 'S':32.065}

# xml_file_path = '/home/ivanagloriosi/Dropbox/UNI/BioInformatica/tut-biomedica/1fpv.xml'
xml_file_path = 'C:\\Users\\ivanagloriosi\\Documents\\My Dropbox\\UNI\\BioInformatica\\tut-biomedica\\1fpv.xml'

tree = etree.parse(xml_file_path)
ns = etree.FunctionNamespace('http://pdbml.pdb.org/schema/pdbx-v32.xsd')
ns.prefix = 'PDBx'

nodes = tree.xpath("//PDBx:atom_site")

atoms = []

for node in nodes:
    sym = str(node.xpath("./PDBx:type_symbol/text()")[0])
    x = float(node.xpath('./PDBx:Cartn_x/text()')[0])
    y = float(node.xpath('./PDBx:Cartn_y/text()')[0])
    z = float(node.xpath('./PDBx:Cartn_z/text()')[0])
    atom = (sym, (x, y, z))

    atoms.append(atom)

##eul_mat_list = []

##for atom in atoms:
##    x = float(atom[1][0]) 
##    y = float(atom[1][1])
##    z = float(atom[1][2])
##
##    xx = x * x
##    xy = x * y
##    xz = x * z
##
##    yy = y * Y
##    yz = y * z
##
##    zz = z * z
##
##    m = masse_atomiche[atom[0]]
##
##    eul_mat = (xx, xy, xz, x, yy, yz, y, zz, z, m)
##    
##    eul_mat_list.append(eul_mat)

eul_mat_list = [(x*x, x*y, x*z, x, y*y, y*z, y, z*z, z, masse_atomiche[sym]) for (sym, (x,y,z)) in atoms]

sum_a = array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

for e in eul_mat_list:
    sum_a += array(e)

massa = sum_a[9]
x_sum = sum_a[3]
y_sum = sum_a[6]
z_sum = sum_a[8]

baricentro = (x_sum/massa, y_sum/massa, z_sum/massa)

atoms_new_coord = [(x-baricentro[0], y-baricentro[1], z-baricentro[2]) for (sym, (x,y,z)) in atoms]

print atoms_new_coord

