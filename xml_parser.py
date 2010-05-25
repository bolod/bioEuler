from lxml import etree

tree = etree.parse('/home/ivanagloriosi/Dropbox/UNI/BioInformatica/tut-biomedica/1fpv.xml')
ns = etree.FunctionNamespace('http://pdbml.pdb.org/schema/pdbx-v32.xsd')
ns.prefix = 'PDBx'

nodes = tree.xpath("//PDBx:atom_site[./PDBx:label_atom_id/text() = 'N' or ./PDBx:label_atom_id/text() = 'C' or ./PDBx:label_atom_id/text() = 'CA']")

atoms = []

for node in nodes:
    sym = node.xpath("./PDBx:label_atom_id/text()")[0]
    x = node.xpath('./PDBx:Cartn_x/text()')[0]
    y = node.xpath('./PDBx:Cartn_y/text()')[0]
    z = node.xpath('./PDBx:Cartn_z/text()')[0]
    atom = (sym, (x, y, z))

    atoms.append(atom)


centro_di_massa = [0.0, 0.0, 0.0]
for a in atoms:
    centro_di_massa[0] += float(a[1][0])
    centro_di_massa[1] += float(a[1][1])
    centro_di_massa[2] += float(a[1][2])

dim = len(atoms)

centro_di_massa[0] = centro_di_massa[0]/dim
centro_di_massa[1] = centro_di_massa[1]/dim
centro_di_massa[2] = centro_di_massa[2]/dim
    
print centro_di_massa




