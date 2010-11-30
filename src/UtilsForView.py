from pyplasm import *
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

class ProViewer():
   def __init__(self):
      pass 
    
   def view(self, batches):
        """
        Render the given Batches

        Parameters
        ----------
        batches: list of Batch
            list of all of the batches to view

        Returns
        -------
        cuboid : PLaSM cuboid
            PLaSM cuboid with given size and coords.

        """
         # organize the batch in a loose octree
        octree = Octree(batches)
        # create the viewer and run it
        viewer = Viewer(octree)
        viewer.Run()

class Batcher():
    
    def __init__(self, atom_list, color=WHITE):
       self.atom_list = atom_list
       self.color = color
       self.unitCylinder = Batch.openObj("../resources/cylinder4x27.obj")[0]
       self.unitSphere = Batch.openObj("../resources/sphere18x27.obj")[0]

    def setColor(self, color):
        self.color = color

    def __sphere(self, c, color=WHITE, sx=0.15, ):
        batchSphere = Batch(self.unitSphere)
        batchSphere.matrix = Mat4f.translate(*c)*Mat4f.scale(sx,sx,sx)
        batchSphere.diffuse = color
        return batchSphere

    def __cylinder(self, p1, p2, color=WHITE, sx=0.15):
        vect = VECTDIFF([p2,p1])
        qz = UNITVECT(vect)
        qx = UNITVECT(VECTPROD([ vect,[0,0,1] ]))
        qy = VECTPROD([ qz,qx ])
        Rot = TRANS([qx,qy,qz])
        Rot = CAT([ Rot[0]+[0.], Rot[1]+[0.], Rot[2]+[0.], [0.,0.,0.,1.] ])
        h = VECTNORM(vect)

        batchCylinder = Batch(self.unitCylinder)
        batchCylinder.matrix = Mat4f.translate(*p1) * Mat4f(*Rot) * Mat4f.scale(sx,sx,h)
        batchCylinder.diffuse = color
        return batchCylinder

    def __stickAndBall(self, atom_list, ball_sx, color):
        batches = []
        for i in range(atom_list.get_size() - 1):
            coords1 = atom_list.list[i].get_coords()
            coords2 = atom_list.list[i+1].get_coords()
            batches.append(self.__sphere(coords1, color, ball_sx))
            batches.append(self.__cylinder(coords1, coords2, color))
        batches.append(self.__sphere(atom_list.list[self.atom_list.get_size() - 1].get_coords(), color, ball_sx))
        return batches

    def vanDerWaals(self):
        """
        Create batches in the van Der Waals radii type visualization.

        Parameters
        ----------
        None

        Returns
        -------
        batches : list of batch for ProViewer viewer

        """

        batches = []
        ball_sx = 0.15
        for atom in self.atom_list.list:
            ball_sx = atom.get_van_Der_Waals_radius()/100.
            coords = atom.get_coords()
            batches.append(self.__sphere(coords, self.color, ball_sx))
        return batches

    def stick(self):
        """
        Create batches in the van stick type visualization.

        Parameters
        ----------
        None

        Returns
        -------
        batches : list of batch for ProViewer viewer

        """

        ball_sx = 0.15
        return self.__stickAndBall(self.atom_list, ball_sx, self.color)

    def stickAndBall(self):
        """
        Create batches in the van stick and ball type visualization.

        Parameters
        ----------
        None

        Returns
        -------
        batches : list of batch for ProViewer viewer

        """

        ball_sx = 0.3
        return self.__stickAndBall(self.atom_list, ball_sx, self.color)

class Plasm():

    def create_cube(self, coords, size=0.1, color=WHITE):
        """
        Create a PLaSM cuboid with given size and coords.

        Parameters
        ----------
        size: float
            size of the cube

        Returns
        -------
        cuboid : PLaSM cuboid
            PLaSM cuboid with given size and coords.

        """
        return (T([1,2,3])(coords)(CUBOID([size, size, size])))

    def create_cube_list(self, coords_list, size=0.1, color=WHITE):
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
        return [ self.create_cube(coords, size, color) for coords in coords_list ]

    def create_ellipsoid(self, coords, size=0.1, color=WHITE):
        """
        Create a PLaSM ellipsoid of this atom list.

        Parameters
        ----------
        color : integer
            number of the color in the dictionary "colors"

        Returns
        -------
        ellipsoid : colored PLaSM ellipsoid

        """
        return COLOR(color)(T([1,2,3])(coords)((S([1,2,3])(size)(SPHERE(1)([16,16])))))

    def create_convex_hull(self, coords_list, size=0.1, color=WHITE):
        """
        Create a PLaSM convex hull of this atom list.

        Parameters
        ----------
        color : integer
            number of the color in the dictionary "colors"

        Returns
        -------
        convex_hull : colored PLaSM ellipsoid

        """
        return COLOR(color)(JOIN(self.create_cube_list(coords_list, size)))

    def create_polyline(self, coords_list):
        """
        Create a PLaSM polyline of this atom list.

        Parameters
        ----------
        coords_list : array_like, shape(dim, 3)
            coords list
        Returns
        -------
            polyline
        """
        return POLYLINE(coords_list)

    def color(self, plasm_object, color=WHITE):
        """
        Add color to a PLaSM objects.

        Parameters
        ----------
        plasm_object : PLaSM object
            PLaSM object to color
        color : PLaSM color
            PLaSM color
        """
        colors = {
            '1': WHITE,
            '2': RED,
            '3': GREEN,
            '4': BLUE }

        return COLOR(color)(plasm_object)

    def view(self, plasm_object_list):
        """
        View a list of PLaSM objects.

        Parameters
        ----------
        plasm_object_list : PLaSM object list
            list of PLaSM objects
        """
        VIEW(STRUCT(plasm_object_list))