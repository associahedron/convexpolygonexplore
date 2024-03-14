import numpy as np
import matplotlib.pyplot as plt

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
            is_extreme = is_extreme and (w[k] == 0 or w[k] == n-k-np.sum(w[k+1:]))
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
        w1 = self.w
        w2 = other.w
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

    def get_edges(self, min_idx=0):
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
    
    def get_triangles(self):
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


    def draw(self, d, c, options=None):
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
        self.XTheta = XTheta
        plt.scatter(X[:, 0], X[:, 1], c='k', zorder=100) 
        if options["show_codeword"]:
            for i in range(N):
                if i >= min_idx:
                    c = 'k'
                    if i in options["bold_idxs"]:
                        c = 'C1'
                    plt.text(XTheta[i, 0], XTheta[i, 1], "{}".format(w[i]), c=c, ha="center", va="center")
        X = np.concatenate((X, X[0, :][None, :]), axis=0)
        plt.plot(X[:, 0], X[:, 1], c='k')

        ## Step 2: Draw polygon edges
        w = np.concatenate((w, np.array([0, 0])))
        i = N-1
        for (i, j) in self.get_edges(min_idx):
            plt.plot(X[[i, j], 0], X[[i, j], 1], c='k')
        plt.axis("off")


class Associahedron:
    def __init__(self, n, opts=None, draw=True):
        if not opts:
            opts = {}
        if not "diameter" in opts:
            opts["diameter"] = 1
        if not "g_x_offset" in opts:
            opts["g_x_offset"] = 0
        if not "g_y_offset" in opts:
            opts["g_y_offset"] = 0
        w = np.zeros(n, dtype=int)
        self.stack_index = np.zeros(n, dtype=int)
        self.codewords = []
        self.make_stack_rec(w, n-1, opts["g_x_offset"], opts["g_y_offset"], opts, draw)


    def make_stack_rec(self, w, d, g_x_offset, y_offset, opts, draw=True):
        n = w.size
        diam = opts["diameter"]
        dy = -0.1*diam
        h = w.size-d-np.sum(w[d+1:])+1
        x_offset = g_x_offset + 1.5*diam*(n-d-1)
        if draw:
            plt.text(x_offset, y_offset+1.2*diam/2, "d = {}, h = {}".format(d, h), ha="center", va="center")
        y1 = y_offset+1.4*diam/2
        n_items = 0
        vals = list(range(h))
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
                self.codewords.append(c)
                if draw:
                    c.draw(diam, np.array([x_offset, y_offset+dy]), {
                        "bold_idxs":set([1])
                    })
                    s = "".join([str(u) for u in wi])
                    plt.text(x_offset+1.2*diam/2, y_offset+dy, s, va="center")
                    # Indicate quad where flip happened
                    if len(self.codewords) > 1:
                        e1 = c.get_edges()
                        c2 = self.codewords[-2]
                        e2 = c2.get_edges()
                        e = np.array(list(e2.difference(e1)), dtype=int).flatten()
                        plt.plot(c.X[e, 0], c.X[e, 1], c='k', linewidth=1, linestyle='--')
                        for idx in np.where(c.w != c2.w)[0]:
                            plt.text(c.XTheta[idx, 0], c.XTheta[idx, 1], "{}".format(c.w[idx]), ha="center", va="center", c='C3')

            else:
                if draw:
                    c = Codeword(wi)
                    c.draw(diam, np.array([x_offset, y_offset+dy]), {
                        "min_idx":d,
                        "bold_idxs":set([d])
                    })
                    s = "".join([str(u) for u in wi[d:]])
                    plt.text(x_offset, y_offset+diam/3.2+dy, s, ha="center", va="center")
                ni = self.make_stack_rec(wi, d-1, g_x_offset, y_offset, opts, draw)
            n_items += ni
            y_offset -= diam*1.5*ni
        y2 = y_offset + 1.4*diam/2
        x1 = x_offset - 1.4*diam/2
        x2 = x1 + 1.5*diam
        if d == 1:
            x2 += 0.1*n*diam
        if draw:
            xbox = [x1, x1, x2, x2, x1]
            ybox = [y1, y2, y2, y1, y1]
            plt.plot(xbox, ybox, c='k')
            if self.stack_index[d]%2 == 0:
                plt.fill(xbox, ybox, c=[0.9]*3)
        return n_items


def make_octagon_stack():
    """
    As an example, make a stack of octagons
    """
    plt.figure(figsize=(20, 400))
    a6 = Associahedron(6, draw=True)
    plt.axis("equal")
    plt.savefig("octagonstack.svg", bbox_inches='tight')    

def make_n_stack():
    fac = 1
    plt.figure(figsize=(fac*20, fac*32)) #*42/14
    a2 = Associahedron(2, draw=True)
    a3 = Associahedron(3, {"g_y_offset":-14}, draw=True)
    a4 = Associahedron(4, {"g_x_offset":5}, draw=True)
    
    plt.axis("equal")
    plt.savefig("stacks.svg", bbox_inches='tight')    


make_octagon_stack()
#make_n_stack()
    
"""
c = Codeword([0, 1, 2, 2, 0, 0])
c.draw(1, np.array([0, 0]))
T = c.get_triangles()
print(T)
for ijk in T:
    idx = np.array(list(ijk), dtype=int)
    x = np.mean(c.X[ijk, :], axis=0)
    plt.scatter(x[0], x[1])
plt.show()
"""