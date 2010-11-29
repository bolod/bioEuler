class BioTree():
    def __init__(self, bio_node, parent=None):
        self.bio_node = bio_node
        self.parent = parent
        if (bio_node != None):
            self.left = BioTree(None)
            self.right = BioTree(None)

    def __repr__(self):
        return str([str(value) for value in self.vectorialize()])

    def is_empty(self):
        return self.bio_node == None

    def has_childs(self):
        return not self.is_empty() and self.left != None and self.right != None
    
    def set_node(self, bio_node):
        self.bio_node = bio_node
        return self

    def set_left(self, bio_tree):
        self.left = bio_tree
        bio_tree.parent = self
        return self.left

    def set_right(self, bio_tree):
        self.right = bio_tree
        bio_tree.parent = self
        return self.right

    def set_parent(self, parent):
        self.parent = parent
        return self.parent

    def get_level(self):
        if(self.parent == None):
            return 0
        else:
            return 1 + self.parent.get_level()

    def get_node(self):
        return self.bio_node

    def get_left(self):
        return self.left

    def get_right(self):
        return self.right

    def get_parent(self):
        return self.parent

    def vectorialize(self):
        vector = []

        tree = [self]
        while tree:
            node = tree.pop(0)
            if(node.get_node() != None):
                vector.append(node.get_node())
                if (node.has_childs()):
                    tree.extend([node.get_left(), node.get_right()])

        return vector

    def switch(self, bio_tree, distance):
        if (not(self.has_childs() and bio_tree.has_childs())):
            return

        (left1, right1) = (self.get_left(), self.get_right())
        (left2, right2) = (bio_tree.get_left(), bio_tree.get_right())

        dl1l2 = distance(left1.get_node().get_value(), left2.get_node().get_value())
        dr1r2 = distance(right1.get_node().get_value(), right2.get_node().get_value())
        dr1l2 = distance(right1.get_node().get_value(), left2.get_node().get_value())
        dl1r2 = distance(left1.get_node().get_value(), right2.get_node().get_value())

        if (dl1l2 + dr1r2 > dr1l2 + dl1r2):
            bio_tree.set_left(right1)
            bio_tree.set_right(left1)

        left1.switch(left2, distance)
        right1.switch(right2, distance)

