from UtilsForView import ProViewer
import os
from PDBReader import PDBReader
from AtomList import *
from numpy import *
from pyplasm import *

repo = os.path.join("..", "repo")
imprints = os.path.join("..", "imprints")

def orient(al):
    al.align()
    scale_values = [1, 1, 1]

    for i in [0, 1, 2]:
        [up, down] = al.split(i)
        up_val = up.get_turning_radius_weight()
        down_val = down.get_turning_radius_weight()
        print al.get_name(), 'up_val:', up_val, 'down_val:', down_val
        if  up_val < down_val:
            scale_values[i] = -1

    print al.get_name(), scale_values
    return al.scale(scale_values)


if __name__ == "__main__":
    #prot1 = '1xwc'   #the first protein
    #prot2 = '3trx'   #the second protein

    prot1 = '2vb1'   #the first protein
    prot2 = '1alc'   #the second protein

    #prot1
    pdbReader = PDBReader(prot1, repo)
    protein1 = pdbReader.get_protein_from_PDB()

    #prot2
    pdbReader = PDBReader(prot2, repo)
    protein2 = pdbReader.get_protein_from_PDB()

    al_backbone1 = protein1.get_backbone().align()
    al_backbone2 = protein2.get_splitted_backbone()[0].align()

    all_atoms1 = protein1.get_all_atoms().align()
    all_atoms2 = protein2.get_all_atoms().align()



    ss = [[-1, -1, -1],
    	  [1, -1, -1],
    	  [-1, 1, -1],
    	  [1, 1, -1],
    	  [-1, -1, 1],
    	  [1, -1, 1],
    	  [-1, 1, 1],
    	  [1, 1, 1]]


    pv = ProViewer()

    for s in ss:
        print '##########', s, '##########'
        al_backbone2.scale(s)
        pv.view(all_atoms1.get_proviewer_batches(GREEN, "van der Waals") +
        all_atoms2.get_proviewer_batches(RED, "van der Waals"))
        al_backbone2.scale(s)

#    pv.view(al_backbone1.get_proviewer_batches(GREEN, "stick and ball") +
#        al_backbone2.get_proviewer_batches(RED, "stick and ball"))

    



#    s = STRUCT(
#        [al_backbone1.plasm_polyline(RED)]
#
#        + [ al_backbone2.plasm_polyline(GREEN)]
#    )
#
#    VIEW(s)

#    ss = [[-1, -1, -1],
#    	  [1, -1, -1],
#    	  [-1, 1, -1],
#    	  [1, 1, -1],
#    	  [-1, -1, 1],
#    	  [1, -1, 1],
#    	  [-1, 1, 1],
#    	  [1, 1, 1]]
#
#
#
#
#    for s in ss:
#        s = STRUCT(
#            [al_backbone1.plasm_polyline(RED)]
#
#            + [ S([1,2,3])(s)(al_backbone2.plasm_polyline(GREEN))]
#        )
#
#        VIEW(s)
