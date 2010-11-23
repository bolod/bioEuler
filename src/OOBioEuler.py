"""
OOBioEuler Project
Object Oriented BioEuler Project
"""
from Clusterizer import *

class BioEuler:
    def __init__(self, max_level=5, max_size_ratio=10, min_size=10, max_euler_distance=200000):
        self.max_level = max_level
        self.max_size_rate = max_size_ratio
        self.min_size = min_size
        self.max_euler_distance = max_euler_distance   
        
    def create_imprints_from_file(self, file, type, max_level=5):

        list = []
        list_impt = []
        a = []
        FILE = open(file,"r")
        iter(FILE)

        for line in FILE:
            pdbReader = PDBReader(line[:4])
            protein = pdbReader.get_protein_from_PDB()
            splitted_line = line.split(" ")
            name = str(splitted_line[0].strip('\n\t'))
            chain_num = (len(splitted_line)-1)/2
            atom_list = AtomList()
           
            if (chain_num >= 1):
                start = 0
                end = 0
                i = 1
                while(chain_num > 0):
                    start = int(splitted_line[i].strip(('\n\t')))
                    i += 1
                    end = int(splitted_line[i].strip(('\n\t')))
                    i += 1
                    atom_list.add_atom_list( protein.get_typed_atom_list(type).filter(start,end))
                    chain_num -= 1
            else:
                atom_list.add_atom_list(protein.get_typed_atom_list(type))

            list.append(atom_list)
            list[-1].set_name(name)
        FILE.close()

        for atom_list in list:
            impt = atom_list.get_imprint(max_level)
            list_impt.append(impt)
            name_file = atom_list.get_name() + "." + type + "im"
            impt.write_on_file(name_file)

        return list

#    def create_imprints_from_file(self, file, type, max_level=5):
#
#        list = []
#        list_impt = []
#        a = []
#        FILE = open(file,"r")
#        iter(FILE)
#
#        for line in FILE:
#            pdbReader = PDBReader(line[:4])
#            protein = pdbReader.get_protein_from_PDB()
#            a = line.split(" ")
#            name = str(a[0].strip('\n\t'))
#            start = 0
#            end = 0
#            lg = (len(a)-1)/2
#            i=0
#            atom_list = AtomList()
#
#            while(lg>0):
#                print 'Sono nel ciclo'
#                temp = a[1:][i*2:i*2+2]
#                start = int(temp[0].strip(('\n\t')))
#                end = int(temp[1].strip(('\n\t')))
#                i=+1
#                lg=-2
#                atom_list.add_atom_list( protein.get_typed_atom_list(type).filter(start,end))
#
#            list.append(atom_list)
#            list[-1].set_name(name)
#            print line, atom_list
#        FILE.close()
#
#        for atom_list in list:
#            impt = atom_list.get_imprint(max_level)
#            print file, impt
#            list_impt.append(impt)
#            name_file = atom_list.get_name() + "." + type + "im"
#            impt.write_on_file(name_file)
#
#        return list




#    def create_imprints_from_file(self, filename, type, max_level=5):
#        
#        list = []
#        list_impt = [] 
#        a = []
#        file = open(filename, "r")
#        iter(file)
#        
#        for line in file:
#            pdbReader = PDBReader(line[:4])
#            protein = pdbReader.get_protein_from_PDB()
#            a = line.split(" ")
#            name = str(a[0].strip('\n\t'))
#            start = 0
#            end = 0
#            if (len(a)>2):
#                start = int(a[1].strip('\n\t'))
#                end = int(a[2].strip('\n\t')) 
#            prot = protein.get_typed_atom_list(type).filter(start,end)    
#            list.append(prot)
#            list[-1].set_name(name)
#        file.close()
#            
#        for atom_list in list:
#            impt = atom_list.get_imprint(max_level)
#            list_impt.append(impt)
#            name_file = atom_list.get_name()+"." + type + "im"
#            impt.write_on_file(name_file)
#            
#        return list

#    def imprint(self, p1):
#        return p1.get_imprint()





    def compare(self, p1, p2):
        """
        Compares protein p1 to protein p2.

        Verify if the two proteins may be compared and do it if:
        - each protein has sufficient atoms to compare
        - the level of comparison is under a max level

        To compare two proteins it's necessary
        to make a structural alignment by principal axis.

        Returns
        -------
        [out_1, out_2] : [AtomList list, AtomList list]


        """
        def size_rate_test(a, b):
            return max(a.get_size(), b.get_size()) / min(a.get_size(), b.get_size()) < self.max_size_rate

        def size_test(a, b):
            return a.get_size() > self.min_size and b.get_size() > self.min_size

        def euler_distance_test(a, b):
            return a.get_euler_distance(b) < self.max_euler_distance

        def level_test(level):
            return level <= self.max_level

        if not(size_rate_test(p1, p2)):
            return 0

        out_1 = []
        out_2 = []
        p1 = p1.clone()
        p2 = p2.clone()
        q = [[p1, p2, 0]]
        out_1.append(p1)
        out_2.append(p2)
        while q:
            a, b, level = q.pop(0)

            a.align()
            b.align()

            if level_test(level) and size_test(a, b) and euler_distance_test(a, b) :
                a_split = a.split_z()
                b_split = b.split_z()
                q.extend([[a_split[0], b_split[0], level+1], [a_split[1], b_split[1], level+1]])
                out_1.extend([a_split[0], a_split[1]])
                out_2.extend([b_split[0], b_split[1]])
            else:
                q = []

        return [out_1, out_2]


# MAIN
if __name__ == "__main__":
    from PDBReader import PDBReader
    from UtilsForView import *
    from AtomList import *
    from Imprint import *
    from Cluster import *
#
#
#
#    name_01 = "2vb1"
#    name_02 = "3HF4"
#
#
#    p_01 = PDBReader(name_01).get_protein_from_XML()
#    p_02 = PDBReader(name_02).get_protein_from_XML()
#
#    b_01 = p_01.get_backbone()
#    b_02 = p_02.get_backbone()
#
    bioeuler = BioEuler()

#    bioeuler.compare(b_01, b_02)
#
#    b_01.align()
#    b_02.align()

#    plasm = Plasm()

#   plasm.view([CUBOID([0,0,0])] + [b_01.plasm_polyline(WHITE)]
#       + b_01.plasm_cube_list(0.25, RED)
#       + b_02.plasm_cube_list(0.25, GREEN)
#   )

#    imp01 = p_01.get_backbone_imprint(5)
#    imp01.write_file("C:\\Users\\ivanagloriosi\\Documents\\My Dropbox\\UNI\\BioInformatica\\Marino-Spini")
#    print p_01.get_backbone().get_mass()

    
    bioeuler.create_imprints_from_file("init_test.txt", "B", 6)
    imp_vects = create_vectors_from_file("impronte_test.txt")


    input = [x for (n,x) in imp_vects]

#    D = 10**5
    D = linalg.norm(sum(input)/len(input))/3


#    while D > 10**4:
#        print "\n*************\n"
#        print "Distance ", D
#        clusters = []
#        for c in qt_clustering(input, D):
#            clusters.append(QTCluster(c, imp_vects))
#            print clusters[-1]
#        D -= 10000
#        print "\n*************\n"

    dic = dictionary_clustering(imp_vects, 10**10)
    for k in dic.keys():
        print "---------------"
        print k
        tuples = sorted(dic[k], key=lambda x: x[1])
        for t in tuples:
            print t


        

    print '\n\ncomplete successfully!'
