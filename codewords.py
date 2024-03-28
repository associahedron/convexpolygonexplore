import numpy as np
import matplotlib.pyplot as plt

class Node:
    def __init__(self):
        self.left = None
        self.right = None
    
    def inorder(self, x):
        if self.left:
            self.left.inorder(x)
        self.x = x[0]
        x[0] += 1
        if self.right:
            self.right.inorder(x)
        
    def draw(self, ax, opts, depth=0):
        from utils import draw_curve
        x_offset, dx, y_offset, dy = opts["x_offset"], opts["dx"], opts["y_offset"], opts["dy"]
        x1 = x_offset + self.x*dx
        y1 = y_offset - dy*depth
        y2 = y_offset - dy*(depth+1)
        ax.scatter(x1, y1, c='k', zorder=101)
        if self.left:
            x2 = x_offset + self.left.x*dx
            draw_curve([x1, y1], [x2, y2])
            self.left.draw(ax, opts, depth+1)
        if self.right:
            x2 = x_offset + self.right.x*dx
            draw_curve([x1, y1], [x2, y2])
            self.right.draw(ax, opts, depth+1)

class HEdge:
    def __init__(self):
        self.pair = None
        self.next = None
    
    def pair_other(self, other):
        self.pair = other
        other.pair = self

    def get_tree(self):
        root = Node()
        if self.next.pair: # Left node
            root.left = self.next.pair.get_tree()
        else:
            root.left = Node() # Leaf node
        if self.next.next.pair: # Right node
            root.right = self.next.next.pair.get_tree()
        else:
            root.right = Node() # Leaf node
        return root
    
def draw_tree(ax, root, opts=None):
    if not opts:
        opts = {}
    if not "x_offset" in opts:
        opts["x_offset"] = 0
    if not "dx" in opts:
        opts["dx"] = 1
    if not "y_offset" in opts:
        opts["y_offset"] = 0
    if not "dy" in opts:
        opts["dy"] = 1
    x = [0]
    root.inorder(x)
    root.draw(ax, opts, 0)


def are_rotation_neighbors(w1, w2):
    """
    Figure out if two codewords are rotation neighbors, as per Lemma 2 in the paper

    Parameters
    ----------
    w1: ndarray(N)
        First codeword
    w2: ndarray(N)
        Second codeword

    Returns
    -------
    True if w1 and w2 are rotation neighbors, and False otherwise
    """
    diff = w1 - w2
    is_neighbor = True
    # Make sure one entry is 1 and the other is -1
    if np.max(np.abs(diff) == 1) and np.sum(diff) == 0 and np.sum(np.abs(diff)) == 2:
        a = np.argmax(diff == -1)
        b = np.argmax(diff == 1)
        if a < b:
            # Case 1
            for k in range(a+1, b):
                is_neighbor = is_neighbor and np.sum(w1[k:b]) < b-k
            if a == 0:
                is_neighbor = is_neighbor and np.sum(w1[0:b]) >= b-a-1
            else:
                is_neighbor = is_neighbor and np.sum(w1[a:b]) >= b-a
        else:
            # Case 2
            a, b = b, a
            for k in range(a+1, b):
                is_neighbor = is_neighbor and np.sum(w1[k:b]) < b-k
            if a == 0:
                is_neighbor = is_neighbor and np.sum(w1[0:b]) >= b-a
            else:
                is_neighbor = is_neighbor and np.sum(w1[a:b]) > b-a
    else:
        is_neighbor = False
    return is_neighbor

class Codeword:
    def __init__(self, w):
        if type(w) is int:
            N = w
            self.w = np.zeros(N, dtype=int)
            self.w[0] = N-1
        else:
            self.w = np.array(w, dtype=int)
    
    def get_num_extreme(self):
        """
        Return the number of "length i suffixes" (indices 0:i) of this codeword 
        that are in extreme form, i <= 2 <= n (def pg. 15)
        """
        w = self.w
        n = w.size
        is_extreme = w[0] == n-1-np.sum(w[1:])
        num = 0
        k = 1
        while k < n and is_extreme:
            is_extreme = is_extreme and ((w[k] == 0) or (w[k] == n-k-np.sum(w[k+1:])))
            if is_extreme:
                num += 1
            k += 1
        return num

    def is_rotation_neighbor(self, other):
        """
        Figure out if two codewords are rotation neighbors, as per Lemma 2 in the paper

        Parameters
        ----------
        Codeword: other
            Other codeword
        """
        return are_rotation_neighbors(self.w, other.w)

    def get_edges(self, min_idx=0):
        """
        Extract the internal edges from the codeword
        """
        w = self.w
        N = len(w)
        visible = np.ones(N+2)
        i = N-1
        edges = []
        while i >= min_idx:
            wi = w[i]
            j = i+2
            while wi > 0:
                # Find closest visible vertex
                while visible[j] == 0:
                    j += 1
                edges.append((i, j))
                visible[i+1:j] = 0
                wi -= 1
                j += 1
            i -= 1
        return set(edges)
    
    def get_tris(self):
        """
        A quick n' dirty O(N^3) method to extract triangles from an edge list
        """
        N = len(self.w)+2
        edges = self.get_edges()
        edges.add((0, N-1))
        for i in range(N-1):
            edges.add((i, i+1))
        tris = []
        for i in range(N):
            for j in range(N+1):
                if (i, j) in edges:
                    for k in range(N+2):
                        if (j, k) in edges and (i, k) in edges:
                            tris.append((i, j, k))
        return tris

    def get_tree(self):
        """
        Create the tree from this codeword, using the half-edge as an
        intermediate data structure for the triangulation
        """
        N = len(self.w)
        tris = self.get_tris()
        hedges = {}
        for ijk in tris:
            hs = []
            for x in range(3):
                e = (ijk[x], ijk[(x+1)%3])
                h = HEdge()
                hedges[e] = h
                hs.append(h)
            for i in range(3):
                hs[i].next = hs[(i+1)%3]
        for ijk in tris:
            for x in range(3):
                e1 = (ijk[x], ijk[(x+1)%3])
                e2 = (ijk[(x+1)%3], ijk[x])
                if e2 in hedges:
                    hedges[e1].pair_other(hedges[e2])
        return hedges[(N+1, 0)].get_tree()

    def draw(self, ax, d, c, options=None):
        """
        Parameters
        ----------
        d: float
            Diameter of circle in which the polygon is inscribed
        c: ndarray(2)
            Center of circle in which the polygon is inscribed
        options: {
            show_codeword: bool
                If True, show codeword values at vertices (default True)
            min_idx: int
                Minimum index from which to draw edges or numbers
            bold_idxs: set
                Bold the indices in this set (useful for stacks)
            bold_color
                Color to draw bolded items
            color:
                Color to draw polygon (default 'k')
            circled_vertices: list of int
                Indices of vertices to circle
            dotted_edges: list of [int, int]
                Edges to draw dotted (not necessarily in the triangulation)
            dotted_color: 
                Color to draw dotted edges if they exist
            draw_index: bool
                If True, draw indices instead of values
        }
        """
        if not options:
            options = {}
        if not "show_codeword" in options:
            options["show_codeword"] = True
        if not "min_idx" in options:
            options["min_idx"] = 0
        if not "bold_idxs" in options:
            options["bold_idxs"] = set([])
        if not "bold_color" in options:
            options["bold_color"] = 'C1'
        if not "color" in options:
            options["color"] = 'k'
        if not "circled_vertices" in options:
            options["circled_vertices"] = []
        if not "dotted_edges" in options:
            options["dotted_edges"] = []
        if not "dotted_color" in options:
            options["dotted_color"] = 'k'
        if not "draw_index" in options:
            options["draw_index"] = False

        min_idx = options["min_idx"]
        ## Step 1: Draw polygon boundary
        r = d/2
        w = self.w
        N = len(w)
        theta = np.linspace(0, 2*np.pi, N+3)[0:N+2] + np.pi/2 + np.pi/(N+2)
        X = np.zeros((N+2, 2))
        X[:, 0] = r*np.cos(theta) + c[0]
        X[:, 1] = r*np.sin(theta) + c[1]
        self.X = X
        XTheta = np.array(X)
        XTheta[:, 0] += r*0.2*np.cos(theta)
        XTheta[:, 1] += r*0.2*np.sin(theta)
        ax.scatter(X[:, 0], X[:, 1], c=options["color"], zorder=100) 
        rg = N
        if options["draw_index"]:
            rg = N+2
        if options["show_codeword"]:
            for i in range(rg):
                if i >= min_idx:
                    c = options["color"]
                    weight = "regular"
                    if i in options["bold_idxs"]:
                        c = options["bold_color"]
                        weight = "black"
                    val = i
                    if not options["draw_index"]:
                        val = w[i]
                    ax.text(XTheta[i, 0], XTheta[i, 1], "{}".format(val), c=c, ha="center", va="center", weight=weight)
        X = np.concatenate((X, X[0, :][None, :]), axis=0)
        ax.plot(X[:, 0], X[:, 1], c=options["color"])

        ## Step 2: Draw polygon edges
        w = np.concatenate((w, np.array([0, 0])))
        i = N-1
        for (i, j) in self.get_edges(min_idx):
            ax.plot(X[[i, j], 0], X[[i, j], 1], c=options["color"])

        ## Step 3: Circle any vertices
        for idx in options["circled_vertices"]:
            ax.scatter(XTheta[idx, 0], XTheta[idx, 1], s=200, facecolors='none', edgecolors='C3', zorder=101)

        ## Step 4: Draw any dotted edges
        for [i, j] in options["dotted_edges"]:
            ax.plot(X[[i, j], 0], X[[i, j], 1], c=options["dotted_color"], linestyle='--')


class Associahedron:
    def __init__(self, n, opts=None, ax=None):
        if not opts:
            opts = {}
        if not "diameter" in opts:
            opts["diameter"] = 1
        if not "g_x_offset" in opts:
            opts["g_x_offset"] = 0
        if not "g_y_offset" in opts:
            opts["g_y_offset"] = 0
        if not "draw_tree" in opts:
            opts["draw_tree"] = False
        w = np.zeros(n, dtype=int)
        self.stack_index = np.zeros(n, dtype=int)
        self.codewords = []
        self.last_codeword = None
        self.make_stack_rec(w, n-1, opts["g_x_offset"], opts["g_y_offset"], opts, ax)


    def make_stack_rec(self, w, d, g_x_offset, y_offset, opts, ax=None):
        n = w.size
        diam = opts["diameter"]
        dy = -0.1*diam
        h = w.size-d-np.sum(w[d+1:])+1
        x_offset = g_x_offset - 1.5*diam*(n-d-1)
        if ax:
            ax.text(x_offset, y_offset+1.2*diam/2, "d = {}, h = {}".format(d, h), ha="center", va="center")
        y1 = y_offset+1.4*diam/2
        n_items = 0
        vals = list(range(h))
        stackorder = np.array(self.stack_index, dtype=int)
        if self.stack_index[d]%2 == 1:
            vals = reversed(vals)
        self.stack_index[d] += 1
        for val in vals:
            wi = np.array(w)
            wi[d] = val
            if d == 1:
                ## Base case
                wi[0] = n-1-np.sum(wi[1:])
                ni = 1
                c = Codeword(wi)
                self.codewords.append(dict(
                    c=c,
                    s=stackorder,
                    x=x_offset,
                    y=y_offset
                ))
                if ax:
                    dotted_edges = []
                    circled_vertices = []
                    if len(self.codewords) > 1:
                        # Indicate quad where flip happened
                        e1 = c.get_edges()
                        c2 = self.last_codeword
                        e2 = c2.get_edges()
                        dotted_edges = np.array(list(e2.difference(e1)), dtype=int)
                        circled_vertices = np.where(c.w != c2.w)[0]
                    c.draw(ax, diam, np.array([x_offset, y_offset+dy]), {
                        "bold_idxs":set([1]),
                        "circled_vertices": circled_vertices,
                        "dotted_edges":dotted_edges
                    })
                    s = "".join([str(u) for u in wi])
                    ax.text(x_offset-2*diam-n*0.07, y_offset+dy, s, va="center")
                    if opts["draw_tree"]:
                        root = c.get_tree()
                        draw_tree(ax, root, {"x_offset":x_offset-diam*1.7, "y_offset":y_offset+0.35*diam, "dx":diam/(2*n+1), "dy":diam/(n+1)})
                self.last_codeword = c
            else:
                if ax:
                    c = Codeword(wi)
                    c.draw(ax, diam, np.array([x_offset, y_offset+dy]), {
                        "min_idx":d,
                        "bold_idxs":set([d])
                    })
                    s = "".join([str(u) for u in wi[d:]])
                    ax.text(x_offset, y_offset+diam/3.2+dy, s, ha="center", va="center")
                ni = self.make_stack_rec(wi, d-1, g_x_offset, y_offset, opts, ax)
            n_items += ni
            y_offset -= diam*1.5*ni
        y2 = y_offset + 1.4*diam/2
        x1 = x_offset - 1.5*diam/2
        x2 = x1 + 1.5*diam
        if d == 1:
            x1 -= 0.1*n*diam
            if opts["draw_tree"]:
                x1 -= diam*1.2 # Make room for tree
        if ax:
            xbox = [x1, x1, x2, x2, x1]
            ybox = [y1, y2, y2, y1, y1]
            ax.plot(xbox, ybox, c='k')
            if self.stack_index[d]%2 == 0:
                ax.fill(xbox, ybox, c=[0.9]*3)
        return n_items


def make_octagon_stack():
    """
    As an example, make a stack of octagons
    """
    plt.figure(figsize=(20, 400))
    ax = plt.subplot(111)
    a6 = Associahedron(6, {"draw_tree":True}, ax=ax)
    plt.axis("equal")
    plt.savefig("octagonstack.svg", bbox_inches='tight')    

def make_n_stack():
    fac = 1
    plt.figure(figsize=(fac*20, fac*32)) #*42/14
    ax = plt.subplot(111)
    a2 = Associahedron(2, {"draw_tree":True}, ax=ax)
    a3 = Associahedron(3, {"g_y_offset":-14, "draw_tree":True}, ax=ax)
    a4 = Associahedron(4, {"g_x_offset":5, "draw_tree":True}, ax=ax)
    
    plt.axis("equal")
    plt.savefig("stacks.svg", bbox_inches='tight')    


def non_rotation_example():
    ## Passes check 1 but fails check 2
    #c1 = Codeword([1, 1, 0, 2, 1, 0])
    #c2 = Codeword([1, 0, 1, 2, 1, 0])
    
    #c1 = Codeword([2, 0, 0, 1, 1, 1])
    #c2 = Codeword([2, 1, 0, 0, 1, 1])
    
    #c1 = Codeword([1, 3, 0, 0, 1, 0])
    #c2 = Codeword([1, 2, 0, 0, 2, 0])
    
    #c1 = Codeword([1, 1, 3, 0, 0, 0])
    #c2 = Codeword([1, 1, 2, 0, 0, 1])

    ## Fails check 1 but passes check 2
    #c1 = Codeword([0, 1, 1, 1, 1, 1])
    #c2 = Codeword([1, 1, 1, 1, 0, 1])

    ## Fails both checks
    #c1 = Codeword([3, 1, 0, 0, 1, 0])
    #c2 = Codeword([3, 0, 0, 0, 1, 1])

    ## Passes both checks
    c1 = Codeword([0, 5, 0, 0, 0, 0])
    c2 = Codeword([0, 4, 0, 0, 0, 1])

    circled_vertices = np.where(c1.w != c2.w)[0]
    e1 = c1.get_edges()
    e2 = c2.get_edges()
    dotted_edges = np.array(list(e1.difference(e2)), dtype=int)
    ax = plt.figure(figsize=(6, 6))
    c1.draw(ax, 1, np.array([0, 0]), {"color":"C0", "circled_vertices":circled_vertices})
    c2.draw(ax, 1, np.array([1.5, 0]), {"color":"C1", "circled_vertices":circled_vertices, "dotted_edges":dotted_edges, "dotted_color":"C0"})
    plt.axis("equal")
    plt.savefig("RotExample.svg", bbox_inches='tight')

def quantify_extremes():
    ax = plt.figure(figsize=(12, 3))
    c1 = Codeword([0, 4, 0, 0, 1, 0])
    n1 = c1.get_num_extreme()
    c1b = Codeword([0, 3, 0, 0, 2, 0])
    n1b = c1b.get_num_extreme()

    c2 = Codeword([0, 4, 0, 1, 0, 0])
    n2 = c2.get_num_extreme()
    c2b = Codeword([0, 5, 0, 0, 0, 0])
    n2b = c2b.get_num_extreme()

    dx = -0.25
    x = 0
    c1.draw(ax, 1, np.array([x, 0]), {"bold_idxs":np.arange(1, n1+1), "bold_color":"C0", "circled_vertices":[4]})
    s = ",".join([str(i) for i in range(1, n1+1)])
    plt.text(x+dx-0.025*n1, 0.7, "{}-extreme".format(s))
    
    x += 1.7
    c1b.draw(ax, 1, np.array([x, 0]), {"bold_idxs":np.arange(1, n1b+1), "bold_color":"C0"})
    s = ",".join([str(i) for i in range(1, n1b+1)])
    plt.text(x+dx-0.025*n1b, 0.7, "{}-extreme".format(s))
    
    x += 1.7
    c2.draw(ax, 1, np.array([x, 0]), {"bold_idxs":np.arange(1, n2+1), "bold_color":"C0", "circled_vertices":[2]})
    s = ",".join([str(i) for i in range(1, n2+1)])
    plt.text(x+dx-0.025*n2, 0.7, "{}-extreme".format(s))

    x += 1.7
    c2b.draw(ax, 1, np.array([x, 0]), {"bold_idxs":np.arange(1, n2b+1), "bold_color":"C0"})
    s = ",".join([str(i) for i in range(1, n2b+1)])
    plt.text(x+dx-0.025*n2b, 0.7, "{}-extreme".format(s))
    plt.axis("equal")
    plt.savefig("Extreme.svg", bbox_inches='tight')

def tree_example():
    c = Codeword([1, 0, 4, 0, 0, 0])
    T = c.get_tris()
    print(T)

    ax = plt.subplot(121)
    c.draw(ax, 1, np.array([0, 0]), {"draw_index":True})
    for ijk in T:
        x = np.mean(c.X[ijk, :], axis=0)
        plt.scatter(x[0], x[1])
    plt.subplot(122)
    root = c.get_tree()
    draw_tree(root)
    plt.show()

if __name__ == '__main__':
    make_octagon_stack()
    #make_n_stack()
        
    #non_rotation_example()
    #quantify_extremes()

    #tree_example()
