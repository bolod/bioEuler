"""BioEuler Project

 This program does the following things:
 - asks two protein's code and two colors as input
 - downloads the corresponding files in XML format,
 - parses the files to build an atoms list,
 - puts the correspondig proteins on principal axis,
 - splits in two parts the proteins until they are similar
 - draws the two protein on principal axis, the hellipsoids and convex hull

"""

import xml.etree.ElementTree as ET
from numpy import *
from pyplasm import *
import urllib, os

symbol_to_mass = { 
    'C': 12.0107,
    'N': 14.0067,
    'O': 15.99994,
    'S': 32.065,
    'H': 1.00794 
}

symbol_to_radius = {
    'C': 0.70,
    'N': 0.65,
    'O': 0.60,
    'S': 1.00,
    'H': 0.25
}

colors = {
    '1': WHITE,
    '2': RED,
    '3': GREEN,
    '4': BLUE 
}

max_level = 5

threshold = 200000

class Atom():
    def __init__(self, symbol, label, coords):
        self.symbol = symbol
        self.label = label
        self.coords = coords
        self.m = symbol_to_mass[symbol]
        self.radius = symbol_to_radius[symbol]
        
    def get_coords(self):
        """
        Gets the coordinates of this atom.
        
        Returns
        -------
        coords : ndarray, shape (3, )
            coords of this atom
        
        """
        return self.coords
        
    def get_x(self):
        """
        Gets the x coordinate of this atom.
        
        Returns
        -------
        x : float
            x coord of this atom
        
        """
        return self.coords[0]
    
    def get_y(self):
        """
        Gets the y coordinate of this atom.
        
        Returns
        -------
        y : float
            y coord of this atom
        
        """ 
        return self.coords[1]

    def get_z(self):
        """
        Gets the z coordinate of this atom.
        
        Returns
        -------
        z : float
            z coord of this atom
        
        """
        return self.coords[2]

    def get_mass(self):
        """
        Gets the mass of this atom.
        
        Returns
        -------
        mass : float
            mass of this atom
        
        """
        return self.m

    def get_radius(self):
        """
        Gets the x coordinate of this atom.
        
        Returns
        -------
        x : float
            x coord of this atom
        
        """
        return self.radius
    
    def get_euler(self):
        """
        Gets the euler tensor of this atom.
        
        Returns
        -------
        euler : ndarray, shape (3, 3)
            euler tensor of this atom
        
        """
        return array([ coord * self.coords for coord in self.coords ]) * self.m

    def rotate(self, rotation):
        """
        Rotates this atom by the given rotation matrix.
        
        Parameters
        ----------
        rotation: ndarray, shape(3, 3)
            the rotation matrix
        
        Returns
        -------
        atom : Atom
            this atom rotated
            
        """
        self.coords = dot(rotation, self.coords)
        return self
    
    def translate(self, translation):
        """
        Translates this atom by the given translation vector.
        
        Parameters
        ----------
        translation: ndarray, shape(3, )
            the translation vector
        
        Returns
        -------
        atom : Atom
            this atom translated

        """
        self.coords = self.coords - translation
        return self
            
    def get_symbol(self):
        """
        Gets the symbol of this atom.
    
        Returns
        -------
        symbol : String
            symbol of this atom

        """
        return self.symbol

    def cube(self, size, color):
        """
        Create a PLaSM cuboid with a color an put it on this atom coords.

        Parameters
        ----------
        size: float
            size of the cube        
        color : integer
            number of the color in the dictionary "colors"
        
        Returns
        -------
        cuboid : colored PLaSM cuboid traslated to this atom coords

        """
        return COLOR(color)(T([1,2,3])(self.coords)(CUBOID([size, size, size])))

class AtomList():
    def __init__(self, atom_list=[]):
        self.list = atom_list
    
    def add_atom(self, atom):
        """
        Adds the given atom to this atom list.

        Parameters
        ----------
        atom: Atom
            the atom to add
    
        Returns
        -------
        self : AtomList
            this atom list with the given atom added

        """
        self.list.append(atom)
        return self
    
    def rotate(self, rotation):
        """
        Rotates this atom list by the given rotation matrix.

        Parameters
        ----------
        translation: ndarray, shape(3, 3)
            the rotation matrix
    
        Returns
        -------
        self : AtomList
            this atom list rotated

        """
        for atom in self.list: 
            atom.rotate(rotation)
        return self
    
    def translate(self, translation):
        """
        Translates this atom list by the given translation vector.
        
        Parameters
        ----------
        translation: ndarray, shape(3, )
            the translation vector
        
        Returns
        -------
        self : AtomList
            this atom list translated
        
        """
        for atom in self.list: 
            atom.translate(translation)
        return self
    
    def get_mass(self):
        """
        Gets the mass of this atom list.
        
        Returns
        -------
        mass_center : ndarray, shape(3, )
            mass of this atom list
        
        """
        return sum([atom.get_mass() for atom in self.list])

    def get_masscenter(self):
        """
        Gets the mass center of this atom list.
        
        Returns
        -------
        mass_center : ndarray, shape(3, )
            coords of mass center of this atom list
        
        """
        return sum([ atom.get_coords() * atom.get_mass() for atom in self.list], axis=0) / self.get_mass()
    
    def get_centroid(self):
        """
        Gets the centroid of this atom list.
        
        Returns
        -------
        centroid : ndarray, shape(3, )
            coords of centroid of this atom list
        
        """

        return sum([ atom.get_coords() for atom in self.list]) / self.get_n_atoms()
    
    def get_euler(self):
        """
        Gets the Euler tensor of this atom list.
        
        Returns
        -------
        euler_tensor : ndarray, shape(3, 3)
            Euler tensor of this atom list
        
        """
        return array(sum([ atom.get_euler() for atom in self.list ], axis=0)) / self.get_mass()

    def get_n_atoms(self):
        """
        Gets the number of atoms of this atom list.
        
        Returns
        -------
        n_atoms : ndarray, shape(3, 3)
            number of atoms of this atom list
        
        """
        return len(self.list)

    def get_principal_axis(self):
        """
        Gets the principal axis of this atom list.
        
        Returns
        -------
        princ_axis : ndarray, shape(3, 3)
            principal axis of this atom list
        
        """
        def eigenvectors(m):
            return linalg.eig(m)[1]

        def eigenvalues(m):
            return linalg.eigvals(m)

        def get_ordered_matrix(mat):
            return transpose(sorted(transpose(mat), key=lambda x : x[0]))

        euler = self.get_euler()
        eigen_vectors = eigenvectors(euler)
        eigen_values = eigenvalues(euler)
        matrix = vstack((eigen_values, eigen_vectors))
        ordered_matrix = get_ordered_matrix(matrix)
        ord_eigenvalues = vsplit(ordered_matrix,[1])[0]
        ord_eigenvectors = vsplit(ordered_matrix,[1])[1]
        princ_axis = transpose(ord_eigenvectors)

        return princ_axis

    def ellipsoid(self, color):
        dim = linalg.eigvals(self.get_euler())
        return COLOR(color)(S([1,2,3])(dim)(SPHERE(1)([12,24])))

    def align(self):
        """
        Align this atom list by its principal axis.
        
        Do the following steps:
        1. translate the protein with the center of gravity in the origin
        2. rotate the protein by its principal axis
        
        Returns
        -------
        self : AtomList
            this atom list aligned
        
        """
        self.translate(self.get_masscenter())
        self.rotate(self.get_principal_axis())
        return self

    def cube_list(self, size, color):
        """
        Create a PLaSM cuboid's list 
        
        Parameters
        ----------
        size : float
            size of each cube
        color : integer
            number of the color in the dictionary "colors"
             
        Returns
        -------
        cuboid_list : colored PLaSM cuboid's list
        
        """
        return [ atom.cube(size, color) for atom in self.list ]

    def ellipsoid(self, color):
        """
        Create a PLaSM ellipsoid.
        
        Parameters
        ----------
        color : integer
            number of the color in the dictionary "colors"
             
        Returns
        -------
        ellipsoid : colored PLaSM ellipsoid
        
        """
        max_xyz = (abs(array([ atom.get_coords() for atom in self.list ]))).max(axis=0)

        return COLOR(color)(S([1,2,3])(max_xyz * 1.4)(SPHERE(1)([16,16])))

    def split_z(self):
        """
        Splits this atom list in two atom lists by omogeneous plane XY.

        Returns
        -------
        (first, second) : (AtomList, AtomList)
            first AtomList has all atoms which coordinate Z is >= 0
            second AtomList has all atoms which coordinate Z is < 0

        """
        return (
            AtomList([ atom for atom in self.list if atom.get_z() >= 0 ]),
            AtomList([ atom for atom in self.list if atom.get_z() < 0 ]))

    def get_euler_distance(self, atom_list):
        """
        Computes the distance from this atom list to the given atom list
        as the distance of the Euler tensor's diagonals of the atom lists.
        
        Parameters
        ----------
        atom_list : AtomList
            the atom list on which the distance is computed
          
        Returns
        -------
        distance : float
            the distance of the Euler tensor's diagonals of the atom lists
         
        """
        return linalg.norm(diag(self.get_euler()) - diag(atom_list.get_euler()))

    def compare(self, atom_list):
        """
        Compares this atom list with the given atom list.
        
        Verify if the two proteins may be compared and do it if:
        - each protein has sufficient atoms to compare
        - the level of comparison is under a max level
        
        To compare two proteins it's necessary
        to make a structural alignment by principal axis.
        
        """
    
        q = [[self, atom_list, 0]]

        while q:
            a, b, level = q.pop(0)

            if level < max_level and a.get_n_atoms() > 0 and b.get_n_atoms() > 0:
                print 'level ' + str(level)
                a.align()
                b.align()

                if a.get_euler_distance(b) < threshold:
                    a1, a2 = a.split_z()
                    b1, b2 = b.split_z()
                    q.extend([[a1,b1,level+1], [a2,b2,level+1]])
                else:
                    return level

        return level

    def filter_by_symbol(self, symbol_list):
        """
        Filter this atom list by atom symbol.

        Returns
        -------
        self : AtomList
            this atom list filtered

        """
        self.list = [ atom for atom in self.list if atom.get_symbol() in symbol_list ]
        return self

class PDB:
    def get_file(self, name, ext):
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
        file_name = name + "." + ext
        path = os.getcwd()
        file = os.path.join(path, file_name)
        if os.path.exists(file):
            print "File %s already exists" % file
        else:
            url = "http://www.rcsb.org/pdb/files/"
            print 'retrieving %s' % url + file_name
            lines = urllib.urlopen(url + file_name).read()
            open(file_name,'wb').write(lines)
            print "File %s is now here %s" % (file_name, file)

        return file_name

    def get_atom_list(self, file_name):
        """
        Parses a XML file of protein and create an atomList object.
        
        Parameters
        ----------
        filename : string
           the name of the file to parse
        
        Returns
        -------
        atom_list : AtomList
        
        """
        tree = ET.parse(file_name)
        root = tree.getroot()
        namespace = root.tag[:-9]
        search_string = ".//%satom_site" % namespace
        node_list = tree.findall(search_string)

        atom_list = AtomList()

        for node in node_list:
            x = float(node.find(".//%sCartn_x" % namespace).text)
            y = float(node.find(".//%sCartn_y" % namespace).text)
            z = float(node.find(".//%sCartn_z" % namespace).text)
            label = node.find(".//%stype_symbol" % namespace).text
            symbol = label[0]
            coords = [x, y, z]
            atom = Atom(symbol, label, array(coords))
            atom_list.add_atom(atom)

        return atom_list

class Viewer():

    def view(self, plasm_object_list):
        VIEW(STRUCT(plasm_object_list))

# MAIN
if __name__ == "__main__":
    name_01 = '1alc' #'2e1m'
    name_02 = '2vb1' #'2w9i'

    file_01 = PDB().get_file(name_01, "xml")
    file_02 = PDB().get_file(name_02, "xml")

    p1 = PDB().get_atom_list(file_01)
    p2 = PDB().get_atom_list(file_02)

    p1.align()

    Viewer().view(
          p1.cube_list(0.1, BLUE)
        + [p1.ellipsoid(RED)]
    )
