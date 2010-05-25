from lxml import etree
from numpy import array

masse_atomiche = {'C':12.0107, 'N':14.0067, 'O':15.99994, 'S':32.065}

# xml_file_path = '/home/ivanagloriosi/Dropbox/UNI/BioInformatica/tut-biomedica/1fpv.xml'
xml_file_path = 'C:\\Users\\ivanagloriosi\\Documents\\My Dropbox\\UNI\\BioInformatica\\tut-biomedica\\1fpv.xml'

tree = etree.parse(xml_file_path)
ns = etree.FunctionNamespace('http://pdbml.pdb.org/schema/pdbx-v32.xsd')
ns.prefix = 'PDBx'

nodes = tree.xpath("//PDBx:atom_site")

n_atomi = len(nodes)

print n_atomi
