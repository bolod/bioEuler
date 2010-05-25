from elementtree import ElementTree

tree = ElementTree.parse("1fpv.xml")
root = tree.getroot()

print root

namespace="http://pdbml.pdb.org/schema/pdbx-v32.xsd"
searchString = ".//{%s}atom_site" % (namespace, namespace)
atoms = tree.findall(searchString)

for node in atoms:
    print node.find('.//{%s}Cartn_x' % namespace).text,
    print node.find('.//{%s}Cartn_y' % namespace).text,
    print node.find('.//{%s}Cartn_z' % namespace).text




