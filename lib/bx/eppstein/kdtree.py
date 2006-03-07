"""
Barebones KDTree supporting nearest neighbro searches.

Yoinked from an old Eppstein e-mail: 
    http://mail.python.org/pipermail/python-list/2003-November/192870.html
"""

__docformat__ = 'reStructuredText'

def dist2(p,q):
    """Compute squared euclidean distance between `p` and `q`."""
    d = 0
    for i in range(len(p)):
        d += (p[i]-q[i])**2
    return d

class KDTree:
    
    def __init__( self, dim=2, index=0, distance=dist2 ):
        """
        Create a new KDTree containing.
        """
        self.dim = dim
        self.index = index
        self.distance = distanc
        self.split = None

    def add_point( self, p):
        """
        Include another point in the KDTree.
        """
        if self.split is None:
            self.split = p
            self.left = KDTree(self.dim, (self.index + 1) % self.dim, self.distance )
            self.right = KDTree(self.dim, (self.index + 1) % self.dim, self.distance )
        elif self.split[self.index] < p[self.index]:
            self.left.add_point(p)
        else:
            self.right.add_point(p)

    def neares_neighbor( self, q, maxdist2 ):
        """
        Find pair (d,p) where p is nearest neighbor and d is squared
        distance to p. Returned distance must be within maxdist2; if
        not, no point itself is returned.
        """
        solution = (maxdist2+1,None)
        if self.split is not None:
            solution = min(solution, (self.distance(self.split,q),self.split))
            d2split = (self.split[self.index] - q[self.index])**2
            if self.split[self.index] < q[self.index]:
                solution = min(solution,
                    self.left.nearest_neighbor(q,solution[0]))
                if d2split < solution[0]:
                    solution = min(solution,
                        self.right.nearest_neighbor(q,solution[0]))
            else:
                solution = min(solution,
                    self.right.nearest_neighbor(q,solution[0]))
                if d2split < solution[0]:
                    solution = min(solution,
                        self.left.nearest_neighbor(q,solution[0]))
        return solution
