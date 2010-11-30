from UtilsForView import ProViewer, Batcher
import os
from PDBReader import PDBReader
from AtomList import *
from numpy import *
from pyplasm import *

repo = os.path.join("..", "repo")
imprints = os.path.join("..", "imprints")

if __name__ == "__main__":
    prot1 = '1xwc'   #the first protein
    prot2 = '3trx'   #the second protein

    prot3 = '2vb1'   #the third protein
    prot4 = '1alc'   #the fourth protein

    prot5 = '1bb0'   #the fifth protein
    prot6 = '2drp'   #the sixth protein


    #prot1 single chain
    pdbReader = PDBReader(prot1, repo)
    protein1 = pdbReader.get_protein_from_PDB()
    p1cA = protein1.get_chain('A').orient()
    bb1 = protein1.get_backbone().orient()
    batcher_p1cA = Batcher(p1cA, GREEN)
    batcher_bb1 = Batcher(bb1, GREEN)

    #prot2 single chain
    pdbReader = PDBReader(prot2, repo)
    protein2 = pdbReader.get_protein_from_PDB()
    p2cA = protein2.get_chain('A').orient()
    bb2 = protein2.get_backbone().orient()
    batcher_p2cA = Batcher(p2cA, RED)
    batcher_bb2 = Batcher(bb2, RED)

    #prot3 single chain
    pdbReader = PDBReader(prot3, repo)
    protein3 = pdbReader.get_protein_from_PDB()
    p3cA = protein3.get_chain('A').orient()
    bb3 = protein3.get_backbone().orient()
    batcher_p3cA = Batcher(p3cA, YELLOW)
    batcher_bb3 = Batcher(bb3, YELLOW)

    #prot4 single chain
    pdbReader = PDBReader(prot4, repo)
    protein4 = pdbReader.get_protein_from_PDB()
    p4cA = protein4.get_chain('A').orient()
    bb4 = protein4.get_backbone().orient()
    batcher_p4cA = Batcher(p4cA, BLUE)
    batcher_bb4 = Batcher(bb4, BLUE)

    #prot5 - has multiple chains
    pdbReader = PDBReader(prot5, repo)
    protein5 = pdbReader.get_protein_from_PDB()
    p5cA = protein5.get_chain('A').orient()
    bb5 = protein5.get_splitted_backbone()[0].orient()
    batcher_p5cA = Batcher(p5cA, PURPLE)
    batcher_bb5 = Batcher(bb5, PURPLE)

    #prot6 - has multiple chains
    pdbReader = PDBReader(prot6, repo)
    protein6 = pdbReader.get_protein_from_PDB()
    p6cA = protein6.get_chain('A').orient()
    bb6 = protein6.get_splitted_backbone()[0].orient()
    batcher_p6cA = Batcher(p6cA, GRAY)
    batcher_bb6 = Batcher(bb6, GRAY)

    
    #ProViewer
    pv = ProViewer()

    pv.view(batcher_p1cA.vanDerWaals() + batcher_p2cA.vanDerWaals())
    pv.view(batcher_bb1.stickAndBall() + batcher_bb2.stickAndBall())

    pv.view(batcher_p3cA.vanDerWaals() + batcher_p4cA.vanDerWaals())
    pv.view(batcher_bb3.stickAndBall() + batcher_bb4.stickAndBall())

    # AtomList.orinet() works on protein5 and protein6 A-chain (p5cA - p6cA)
    # but fails with protein5 and protein6 backbone only (bb5 - bb6)
    pv.view(batcher_p5cA.vanDerWaals() + batcher_p6cA.vanDerWaals())
    pv.view(batcher_bb5.stickAndBall() + batcher_bb6.stickAndBall())



    # ss = [[-1, -1, -1],
    #       [1, -1, -1],
    #       [-1, 1, -1],
    #       [1, 1, -1],
    #       [-1, -1, 1],
    #       [1, -1, 1],
    #       [-1, 1, 1],
    #       [1, 1, 1]]

    # for s in ss:
    # print '##########', s, '##########'
    # al_backbone2.scale(s)
    # pv.view(b_bb_1.stickAndBall() + b_bb_2.stickAndBall())
    # al_backbone2.scale(s)


