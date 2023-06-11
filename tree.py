# class to represent a node in a Tree
# a Tree(node) has an id, parent, distance to parent, and list of children 
class Tree:
    def __init__(self, id, parent=None,distance=0.0,childrenLst=None,):
        self.id = id
        self.parent = parent
        self.distance = distance
        self.childrenLst = []
        if childrenLst is not None:
            for child in childrenLst:
                self.add_child(child)

    def add_child(self, node, distance):
        assert isinstance(node, Tree)
        self.childrenLst.append(node)
