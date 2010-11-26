from numpy import *
from pyplasm import *
from Imprint import *
from BioTree import *
from BioNode import *

class AtomList():
    def __init__(self, atom_list=[], name=''):
        self.list = atom_list
        self.name = name

    def __repr__(self):
        """
        Prints info of this atom list.

        Returns
        -------
        info : String
            info of this atom list

        """
        s = '{\n'
        s += 'name: \'' + self.name + '\',\n'
        s += 'size: ' + str(self.get_size()) + ',\n'
        s += 'masscenter: ' + str(self.get_masscenter()) + ',\n'
        s += 'centroid: ' + str(self.get_centroid()) + ',\n'
        s += 'euler: \n'
        s += str(self.get_euler()) + '\n'
        s += '}\n'
        return s

    def get_name(self):
        """
        Gets the name of this atom list.

        Returns
        -------
        name : String
            name of this atom list

        """
        return self.name

    def get_atom_list(self):
        """
        Gets the atom list.

        Returns
        -------
        list : list of Atom
            the atom list

        """
        return self.list

    def get_coords_list(self):
        """
        Gets the atom coords list.

        Returns
        -------
        list : list of Atom coords
            the atom coords list

        """
        return [atom.get_coords() for atom in self.list]

    def set_name(self, name):
        """
        Sets the name of this atom list.

        Parameters
        ----------
        name : String
            name of this atom list

        Returns
        -------
        self : AtomList
            this atom list

        """
        self.name = name
        return self

    def set_atom_list(self, atom_list):
        """
        Sets this atom list.

        Parameters
        ----------
        list : list of Atom
            the atom list

        Returns
        -------
        self : AtomList
            this atom list

        """
        self.list = atom_list
        return self

    def clone(self):
        """
        Clones this atom list.

        Returns
        -------
        clone : AtomList
            the clone of this atom list

        """
        clone = AtomList([], self.name)
        clone.list = [atom.clone() for atom in self.list ]
        return clone

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
        if atom != None:
            self.list.append(atom)
        return self

    def add_atom_list(self, atom_list):
        for atom in atom_list.list:
            self.list.append(atom.clone())
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

    def scale(self, values):
        """
        Scales this atom list by the given scale value vector.

        Parameters
        ----------
        values: ndarray, shape(3, )
            the scale value vector

        Returns
        -------
        self : AtomList
            this atom list scaled

        """
        for atom in self.list:
            atom.scale(values)
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

        return sum([ atom.get_coords() for atom in self.list]) / self.get_size()

    def get_euler(self):
        """
        Gets the Euler tensor of this atom list.

        Returns
        -------
        euler_tensor : ndarray, shape(3, 3)
            Euler tensor of this atom list

        """
        return array(sum([ atom.get_euler() * atom.get_mass() for atom in self.list ], axis=0))

    def get_turning_radius(self):
        """
        Gets the turning radius values of this atom list.
        The turning radius of the three axis are:
        [ sqrt(I_xx / M) , sqrt(I_yy / M) , sqrt(I_zz / M) ]
        where [ I_xx, I_yy, I_zz ] is the diagonal of the euler matrix
        and M is the mass of this atom list.

        Returns
        -------
        turning_radius : ndarray, shape(1, 3)
            turning_radius of this atom list
        """
        return array(sqrt(diag(self.get_euler()) / self.get_mass()))

    def get_turning_radius_weight(self):
        """
        Gets the turning radius weight of this atom list.
        The weight is calculates as follow:
        w = sqrt( reduce(mul, turning_radius) / sum(turning_radius))
        
        Returns
        -------
        w : Real
        """
        turning_radius = self.get_turning_radius()
        return sqrt(reduce(mul, turning_radius) / sum(turning_radius))

    def get_size(self):
        """
        Gets the number of atoms of this atom list.

        Returns
        -------
        size : Real
            number of atoms of this atom list

        """
        return len(self.list)

    def get_index(self, atom):
        """
        Gets the index of the given atom in this atom list.

        Returns
        -------
        index : int
            index of the given atom in this atom list

        """
        return self.list.index(atom)

    def get_prev(self, atom):
        """
        Gets the previous atom of the given atom in this atom list if one exists.

        Returns
        -------
        prev : Atom
            previous atom of the given atom in this atom list if one exists.

        """
        try:
            return self.list[self.list.index(atom) - 1]
        except Exception:
            return None

    def get_next(self, atom):
        """
        Gets the next atom of the given atom in this atom list if one exists.

        Returns
        -------
        next : Atom
            next atom of the given atom in this atom list if one exists.

        """
        try:
            return self.list[self.list.index(atom) + 1]
        except Exception:
            return None

    def get_principal_axis(self):
        """
        Gets the principal axis of this atom list.

        Returns
        -------
        ndarray, shape(3, 3)
            principal axis of this atom list

        """
        euler = self.get_euler()
        eigen_vec = transpose(linalg.eig(euler)[1])
        eigen_val = linalg.eigvals(euler)

        return [ eigen_vec[i] for i in eigen_val.argsort() ]

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
        return self.translate(self.get_masscenter()).rotate(self.get_principal_axis())

    def get_ellipsoid_axis(self):
        """
        Get ellipsoid axis.

        Returns
        -------
        axis : ndarray, shape(3, )
            ellipsoid axis

        """
        return (abs(array([ atom.get_coords() for atom in self.list ]))).max(axis=0)

    def split(self, split_coord):
        """
        Splits this atom list in two atom lists by omogeneous plane.

        Parameters
        ----------
        split_coord: int {0, 1, 2}
            the coordinate to split

        Returns
        -------
        [first, second] : [AtomList, AtomList]
            first AtomList has all atoms which coordinate split_coord is >= 0
            second AtomList has all atoms which coordinate split_coord is < 0

        """
        split_list = [AtomList([], self.name + '-split_' + str(split_coord)), AtomList([], self.name + '-split_' + str(split_coord))]
        for atom in self.list:
            split_list[atom.get_coords()[split_coord] < 0].add_atom(atom)
        return split_list

    def split_x(self):
        """
        Splits this atom list in two atom lists by omogeneous plane YZ.
        
        Returns
        -------
        [first, second] : [AtomList, AtomList]
            first AtomList has all atoms which coordinate X is >= 0
            second AtomList has all atoms which coordinate X is < 0
        
        """
        split_list = [AtomList([], self.name + '-split_0'), AtomList([], self.name + '-split_0')]
        for atom in self.list:
            split_list[atom.get_x() < 0].add_atom(atom)
        return split_list

    def split_y(self):
        """
        Splits this atom list in two atom lists by omogeneous plane XZ.

        Returns
        -------
        [first, second] : [AtomList, AtomList]
            first AtomList has all atoms which coordinate Y is >= 0
            second AtomList has all atoms which coordinate Y is < 0

        """
        split_list = [AtomList([], self.name + '-split_1'), AtomList([], self.name + '-split_1')]
        for atom in self.list:
            split_list[atom.get_y() < 0].add_atom(atom)
        return split_list

    def split_z(self):
        """
        Splits this atom list in two atom lists by omogeneous plane XY.

        Returns
        -------
        [first, second] : [AtomList, AtomList]
            first AtomList has all atoms which coordinate Z is >= 0
            second AtomList has all atoms which coordinate Z is < 0

        """
        split_list = [AtomList([], self.name + '-split_2'), AtomList([], self.name + '-split_2')]
        for atom in self.list:
            split_list[atom.get_z() < 0].add_atom(atom)
        return split_list

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
        float
            the distance of the Euler tensor's diagonals of the atom lists

        """
        return linalg.norm(diag(self.get_euler()) - diag(atom_list.get_euler()))

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

    def filter_by_label(self, label_list):
        """
        Filter this atom list by atom label.

        Returns
        -------
        self : AtomList
            this atom list filtered

        """
        self.list = [ atom for atom in self.list if atom.get_label() in label_list ]
        return self

    def get_minmax_ca_distance(self):
        """
        Return a max and min value of distance between two consecutive Ca atom.

        Returns
        -------
        manmax: (min, max)
            min, min distance between two consecutive Ca atoms
            max, max distance between two consecutive Ca atoms
        """
        ca_list = self.clone().filter_by_label("CA")
        indices_list = zip(range(ca_list.get_size() - 1), range(1, ca_list.get_size()))
        diff = [ca_list.list[i1].get_euclid_distance_to(ca_list.list[i2])  for (i1,i2) in indices_list]
        return (min(diff), max(diff))

    def plasm_polyline(self, color=WHITE):
        """
        Create a PLaSM cuboid's list

        Parameters
        ----------
        color : COLOR
            plasm color

        Returns
        -------
        cuboid_list : colored PLaSM cuboid's list

        """
        return COLOR(color)(POLYLINE([ (atom.get_coords()).tolist() for atom in self.list ]))

    def plasm_cube_list(self, size=0.1, color=WHITE):
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
        return [ atom.plasm_cube(size, color) for atom in self.list ]

    def plasm_ellipsoid(self, color=WHITE):
        """
        Create a PLaSM ellipsoid of this atom list.

        Parameters
        ----------
        color : PLaSM color

        Returns
        -------
        ellipsoid : colored PLaSM ellipsoid

        """
        return COLOR(color)(S([1,2,3])(self.get_ellipsoid_axis() * 1.4)(SPHERE(1)([16,16])))

    def plasm_convex_hull(self, color=WHITE):
        """
        Create a PLaSM convex hull of this atom list.

        Parameters
        ----------
        color : PLaSM color

        Returns
        -------
        convex_hull : colored PLaSM ellipsoid

        """
        return COLOR(color)(JOIN(self.plasm_cube_list(0.1, color)))


    def get_proviewer_batches(self, color=WHITE, view_type="stick"):
        """
        Create a batches of this atom list.

        Parameters
        ----------
        color : PLaSM color

        viewType: string
            ["stick" | "stick and ball"]?

        Returns
        -------
        batches : list of batch for ProViewer viewer

        """

        def cylinder(p1, p2, color=WHITE, sx=0.15):
            vect = VECTDIFF([p2,p1])
            qz = UNITVECT(vect)
            qx = UNITVECT(VECTPROD([ vect,[0,0,1] ]))
            qy = VECTPROD([ qz,qx ])
            Rot = TRANS([qx,qy,qz])
            Rot = CAT([ Rot[0]+[0.], Rot[1]+[0.], Rot[2]+[0.], [0.,0.,0.,1.] ])
            h = VECTNORM(vect)

            batchCylinder = Batch(unitCylinder)
            batchCylinder.matrix = Mat4f.translate(*p1) * Mat4f(*Rot) * Mat4f.scale(sx,sx,h)
            batchCylinder.diffuse = color
            return batchCylinder

        def sphere(c, color=WHITE, sx=0.15):
            batchSphere = Batch(unitSphere)
            batchSphere.matrix = Mat4f.translate(*c)*Mat4f.scale(sx,sx,sx)
            batchSphere.diffuse = color
            return batchSphere

        batches = []
        unitCylinder = Batch.openObj("../resources/cylinder4x27.obj")[0]
        unitSphere = Batch.openObj("../resources/sphere18x27.obj")[0]
        ball_sx = 0.15

        if (view_type == "stick and ball"):
            ball_sx = 0.3

        for i in range(self.get_size() - 1):
            coords1 = self.list[i].get_coords()
            coords2 = self.list[i+1].get_coords()
            batches.append(sphere(coords1, color, ball_sx))
            batches.append(cylinder(coords1, coords2, color))

        batches.append(sphere(self.list[self.get_size() - 1].get_coords(), color, ball_sx))

        return batches
    

    def get_bio_tree(self, max_level=5):

        if(max_level <= 0 or self.get_size() <= 2):
            return BioTree(None)

        [atom_list_inf, atom_list_sup] = self.align().split_z()

        bio_tree = BioTree(BioNode(self.get_euler(), self.get_mass()))

        bio_tree_inf = atom_list_inf.get_bio_tree(max_level - 1)
        bio_tree_sup = atom_list_sup.get_bio_tree(max_level - 1)

        if(bio_tree_inf.is_empty() and bio_tree_sup.is_empty()):
            return bio_tree

        if(bio_tree_inf.get_node().get_weight() > bio_tree_sup.get_node().get_weight()):
            bio_tree.set_left(bio_tree_inf)
            bio_tree.set_right(bio_tree_sup)
        else:
            bio_tree.set_left(bio_tree_sup)
            bio_tree.set_right(bio_tree_inf)

        return bio_tree

    def get_imprint(self, max_level):
        return Imprint(self.get_bio_tree(max_level))

    def filter(self, start_res_number=0, end_res_number=0):
        if start_res_number == 0 and end_res_number == 0:
            return self.clone()

        filtered = AtomList([], self.name + '_filter')
        for atom in self.list:
            if atom.get_res_number() > start_res_number and atom.get_res_number() < end_res_number:
                filtered.add_atom(atom.clone())
        return filtered
