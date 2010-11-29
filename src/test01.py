import os
from PDBReader import PDBReader
from UtilsForView import *
from AtomList import *
from Imprint import *
from Cluster import *

repo = os.path.join("..", "repo")
imprints = os.path.join("..", "imprints")

if __name__ == "__main__":
    prot1 = '1alc'
    prot2 = '2vb1'

    #prot1
    pdbReader = PDBReader(prot1, repo)
    protein1 = pdbReader.get_protein_from_PDB()

    #prot2
    pdbReader = PDBReader(prot2, repo)
    protein2 = pdbReader.get_protein_from_PDB()

    al_backbone1 = protein1.get_backbone().align()
    al_backbone2 = protein2.get_backbone().align()


    al_backbone1.scale([1,-1,-1])

    s = STRUCT(
        [al_backbone1.plasm_polyline(RED)]
        + [(al_backbone2.plasm_polyline(GREEN))]
    )

    VIEW(s)

    impt1 = al_backbone1.get_imprint(5)
    name_file = os.path.join(imprints, prot1) + ".Bim"
    impt1.write_on_file(name_file)

    impt2 = al_backbone2.get_imprint(5)
    name_file = os.path.join(imprints, prot2) + ".Bim"
    impt2.write_on_file(name_file)


    imp_vects = create_vectors_from_file("imp_mirrored.txt")


    input = [x for (n,x) in imp_vects]


    dic = dictionary_clustering(imp_vects, 10**10)
    for k in dic.keys():
        print "---------------"
        print k
        tuples = sorted(dic[k], key=lambda x: x[1])
        for t in tuples:
            print t


        

    print '\n\ncomplete successfully!'

    
    
