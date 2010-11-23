class Protein:

    def __init__(self, name, all_atoms, backbone, splitted_backbone, CA_backbone, CA_splitted_backbone):
        self.name = name
        self.all_atoms = all_atoms
        self.backbone = backbone
        self.splitted_backbone = splitted_backbone
        self.CA_backbone = CA_backbone
        self.CA_splitted_backbone = CA_splitted_backbone

    def __repr__(self):
        return self.name + ' ' + str(self.all_atoms)

    def get_name(self):
        return self.name

    def get_typed_atom_list(self,type):
        
        types = { 'A' : self.get_all_atoms,
                  'B' : self.get_backbone,
                  'C' : self.get_CA_backbone
                }
        
        return types[type]()
        
    def get_all_atoms(self):
        return self.all_atoms

    def get_backbone(self):
        return self.backbone

    def get_splitted_backbone(self):
        return self.splitted_backbone

    def get_CA_backbone(self):
        return self.CA_backbone

    def get_CA_splitted_backbone(self):
        return self.CA_splitted_backbone

    def get_all_atoms_imprint(self, max_level):
        return self.all_atoms.get_imprint(self.name, "A", max_level)

    def get_backbone_imprint(self, max_level):
        return self.backbone.get_imprint(self.name, "B", max_level)

    def get_CA_backbone_imprint(self, max_level):
        return self.CA_backbone.get_imprint(self.name, "C", max_level)

