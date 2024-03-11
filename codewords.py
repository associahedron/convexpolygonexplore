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
            bold_first: bool
                If true, draw the first index in red (useful for stacks).
                Default False
        }
        """
        if not options:
            options = {}
        if not "show_codeword" in options:
            options["show_codeword"] = True
        if not "min_idx" in options:
            options["min_idx"] = 0
        if not "bold_first" in options:
            options["bold_first"] = False

        min_idx = options["min_idx"]
        ## Step 1: Draw polygon boundary
        r = d/2
        w = self.w
        N = len(w)
        theta = np.linspace(0, 2*np.pi, N+3)[0:N+2] + np.pi/2 + np.pi/(N+2)
        X = np.zeros((N+2, 2))
        X[:, 0] = r*np.cos(theta) + c[0]
        X[:, 1] = r*np.sin(theta) + c[1]
        plt.scatter(X[:, 0], X[:, 1], c='k') 
        if options["show_codeword"]:
            for i in range(N):
                if i >= min_idx:
                    c = 'k'
                    if i == min_idx and options["bold_first"]:
                        c = 'r'
                    xi = X[i, 0] + r*0.2*np.cos(theta[i])
                    yi = X[i, 1] + r*0.2*np.sin(theta[i])
                    plt.text(xi, yi, "{}".format(w[i]), c=c)
        X = np.concatenate((X, X[0, :][None, :]), axis=0)
        plt.plot(X[:, 0], X[:, 1], c='k')

        ## Step 2: Draw polygon edges
        w = np.concatenate((w, np.array([0, 0])))
        visible = np.ones(N+2)
        i = N-1
        while i >= min_idx:
            wi = w[i]
            j = i+2
            while wi > 0:
                # Find closest visible vertex
                while visible[j] == 0:
                    j += 1
                plt.plot(X[[i, j], 0], X[[i, j], 1], c='k')
                visible[i+1:j] = 0
                wi -= 1
                j += 1
            i -= 1
        plt.axis("off")

def make_stack_rec(w, d, y_offset, opts):
    n = w.size
    diam = opts["diameter"]
    h = w.size-d-np.sum(w[d+1:])+1
    title = "d = {}, h = {}".format(d, h)
    x_offset = 1.5*diam*(n-d-1)

    plt.text(x_offset-diam/2, y_offset+1.2*diam/2, title)
    y1 = y_offset+1.4*diam/2
    n_items = 0
    for i in range(h):
        wi = np.array(w)
        wi[d] = i
        if d == 1:
            ## Base case
            wi[0] = n-1-np.sum(wi[1:])
            ni = 1
            c = Codeword(wi)
            c.draw(diam, np.array([x_offset, y_offset]))
            plt.text(x_offset+diam*1.2/2, y_offset, "{}".format(wi))
        else:
            c = Codeword(wi)
            c.draw(diam, np.array([x_offset, y_offset]), {
                "min_idx":d,
                "bold_first":True
            })
            ni = make_stack_rec(wi, d-1, y_offset, opts)
        n_items += ni
        y_offset -= diam*1.5*ni
    y2 = y_offset + 1.4*diam/2
    x1 = x_offset - 1.2*diam/2
    x2 = x1 + 1.5*diam
    if d == 1:
        x2 += 0.7*diam
    plt.plot([x1, x1, x2, x2, x1], [y1, y2, y2, y1, y1], c='k')
    return n_items


def make_stack(n, opts=None):
    if not opts:
        opts = {}
    if not "diameter" in opts:
        opts["diameter"] = 1
    w = np.zeros(n, dtype=int)
    make_stack_rec(w, n-1, 0, opts)

def make_octagon_stack():
    """
    As an example, make a stack of octagons
    """
    plt.figure(figsize=(20, 400))
    make_stack(6)
    plt.axis("equal")
    plt.savefig("octagonstack.svg", bbox_inches='tight')    

make_octagon_stack()