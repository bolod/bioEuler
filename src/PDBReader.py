from Bio import PDB
from numpy import *
import os
import urllib
from lxml import etree

from Atom import Atom
from AtomList import AtomList
from Protein import Protein

class PDBReader:

    def __init__(self, pdb_id, path="."):
        self.name = pdb_id
        self.path = path

    def _get_file(self, file_name, file_path):
        """
            Downloads the corresponding file of protein's code.

            If the corresponding file is still on disk,
            doesn't' proceed and uses the file on disk.

            Parameters
            ----------
            name : string
                   the name of the file to download
            ext :  string
                   the file extension

            Returns
            -------
            filename : file on disk

            """
        file = os.path.join(file_path, file_name)
        if os.path.exists(file):
            print "File %s already exists" % file
        else:
            url = "http://www.rcsb.org/pdb/files/"
            print 'retrieving %s' % url + file_name
            lines = urllib.urlopen(url + file_name).read()
            open(file, 'wb').write(lines)
            print "File %s is now here %s" % (file_name, file)

        return file_name

    def _create_Atom(self, pdb_atom, res_number):
        label = pdb_atom.get_name()
        symbol = label[0]
        coords = array([float(coord) for coord in pdb_atom.get_coord()])
        atom = Atom(symbol, label, coords, res_number)
        return atom

    def get_protein_from_PDB(self):
        file_name = self.name + ".pdb"
        self._get_file(file_name, self.path)

        p = PDB.PDBParser(PERMISSIVE=1)

        s = p.get_structure(self.name + "_struct", os.path.join(self.path,file_name))

        allAtoms = AtomList([], file_name)
        backbone = AtomList([], file_name + '_backbone')
        splitted_backbone = []
        CA_backbone = AtomList([], file_name + '_CA_backbone')
        CA_splitted_backbone = []

        ppb = PDB.Polypeptide.PPBuilder()
        #consideriamo solamente il primo modello tra quelli presenti nel file
        for pp in ppb.build_peptides(s[0]):
            CA_splitted_backbone.append(AtomList([]))
            splitted_backbone.append(AtomList([]))
            for res in pp:
                for pdb_atom in res.get_unpacked_list():
                    if not(pdb_atom.is_disordered()) or (pdb_atom.is_disordered() and pdb_atom.altloc == 'A'):
                        res_number = res.get_id()[1]
                        allAtoms.add_atom(self._create_Atom(pdb_atom, res_number))
                        if (pdb_atom.get_name() == 'CA'):
                            CA_backbone.add_atom(self._create_Atom(pdb_atom, res_number))
                            CA_splitted_backbone[-1].add_atom(self._create_Atom(pdb_atom, res_number))
                            res = pdb_atom.get_parent()
                            backbone.add_atom(self._create_Atom(res["N"], res_number))
                            backbone.add_atom(self._create_Atom(res["CA"], res_number))
                            backbone.add_atom(self._create_Atom(res["C"], res_number))
                            splitted_backbone[-1].add_atom(self._create_Atom(res["N"], res_number))
                            splitted_backbone[-1].add_atom(self._create_Atom(res["CA"], res_number))
                            splitted_backbone[-1].add_atom(self._create_Atom(res["C"], res_number))

        return Protein(self.name, allAtoms, backbone, splitted_backbone, CA_backbone, CA_splitted_backbone)

    def _get_Atom_from_lxml_element(self, elem):
        label = elem.xpath("./PDBx:label_atom_id/text()")[0]
        symbol = label[0]
        coords = [float(elem.xpath("./PDBx:Cartn_x/text()")[0]), float(elem.xpath("./PDBx:Cartn_y/text()")[0]), float(elem.xpath("./PDBx:Cartn_z/text()")[0])]
        #print "x = " + elem.xpath("./PDBx:Cartn_x/text()")[0] + " y = " + elem.xpath("./PDBx:Cartn_y/text()")[0] + " z = " + elem.xpath("./PDBx:Cartn_z/text()")[0]
        return Atom(symbol, label, coords)

    def get_protein_from_XML(self):
        file_name = self.name + ".xml"
        self._get_file(file_name)

        tree = etree.parse(file_name)
        ns = etree.FunctionNamespace('http://pdbml.pdb.org/schema/pdbx-v32.xsd')
        ns.prefix = 'PDBx'

        ATOM_nodes = tree.xpath("//PDBx:atom_site[./PDBx:group_PDB/text() = 'ATOM' and (./PDBx:label_alt_id/text() = 'A' or ./PDBx:label_alt_id/text() = not(string(.)))]")
        backboneCA_nodes = tree.xpath("//PDBx:atom_site[./PDBx:group_PDB/text() = 'ATOM' \
                        and (./PDBx:label_alt_id/text() = 'A' or ./PDBx:label_alt_id/text() = not(string(.))) \
                        and PDBx:label_atom_id/text() = 'CA']")

        CA_Atoms = [self._get_Atom_from_lxml_element(CA_node) for CA_node in backboneCA_nodes]
        N_Atoms = [self._get_Atom_from_lxml_element(CA_node.xpath("preceding-sibling::PDBx:atom_site[PDBx:label_atom_id/text() = 'N'][1]")[0]) for CA_node in backboneCA_nodes]
        C_Atoms = [self._get_Atom_from_lxml_element(CA_node.xpath("following-sibling::PDBx:atom_site[PDBx:label_atom_id/text() = 'C'][1]")[0]) for CA_node in backboneCA_nodes]

        z = zip(N_Atoms, CA_Atoms, C_Atoms)

        backbone = []
        for trip in z:
            backbone.extend(list(trip))

        seqs = [int(CA_node.xpath("./PDBx:label_seq_id/text()")[0]) for CA_node in backboneCA_nodes]

            #clono la lista tranne il primo elemento...
        seq_cpyd = seqs[1:]
        #un uno di sicuro...
        indexes = []
        #ma oltre al primo?
        #conto gli uno nella sequenza originale
        unos = seqs.count(1)
        #e prendo le lunghezze delle varie sottosequenze
        for i in range(unos):
            try:
                uno_indx = seq_cpyd.index(1) + 1
                indexes.append(uno_indx)
                del seq_cpyd[:uno_indx]
            except ValueError:
                indexes.append(seq_cpyd[-1])

        p_allAtoms = AtomList([self._get_Atom_from_lxml_element(a) for a in ATOM_nodes], file_name)

        p_CA_backbone = AtomList([CA_Atoms], file_name + '_CA_backbone')
        prefix_sum = [int(sum(indexes[:v])) for v in range(0, len(indexes))]
        zipped = zip(prefix_sum, indexes)
        p_CA_splitted_backbone = [AtomList(backbone[i:i + f])  for (i, f) in zipped]

        p_backbone = AtomList(backbone, file_name + '_backbone')
        indexes3 = 3 * array(indexes)
        prefix_sum3 = 3 * array(prefix_sum)
        zipped3 = zip(prefix_sum3, indexes3)
        p_splitted_backbone = [AtomList(backbone[i:i + f])  for (i, f) in zipped3]

        return Protein(self.name, p_allAtoms, p_backbone, p_splitted_backbone, p_CA_backbone, p_CA_splitted_backbone)
