class Codeword {
    constructor(w) {
        this.w = w;
    }

    /**
     * Extract the internal edges from the codeword
     * @param {int} min_idx Minimum index at which to start
     * throwing down edges
     * @returns 
     */
    get_edges(min_idx) {
        const w = this.w;
        const N = w.length;
        const visible = new Int32Array(N+2);
        for (let k = 0; k < N+2; k++) {
            visible[k] = 1;
        }
        let i = N-1;
        let edges = [];
        while (i >= min_idx) {
            console.log("i", i);
            let wi = w[i];
            let j = i+2;
            while (wi > 0) {
                // Find closest visible vertex
                while (j < N+2 && visible[j] == 0) {
                    j += 1;
                }
                edges.push([i, j]);
                for (let k = i+1; k < j; k++) {
                    visible[k] = 0;
                }
                wi -= 1;
                j += 1;
            }
            i -= 1;
        }
        return edges;
    }

    /**
     * @param {svg Element} g SVG element on which to draw this
     * @param {float} d Diameter of circle in which the polygon is inscribed
     * @param {list of [float, float]} c  Center of circle in which the polygon is inscribed
     * @param {object} options 
     *  {
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
     */
    draw(g, d, c, options) {
        if (options == undefined) {
            options = {};
        }
        if (!("show_codeword" in options)) {
            options["show_codeword"] = true;
        }
        if (!("min_idx" in options)) {
            options["min_idx"] = 0;
        }
        if (!("bold_idxs" in options)) {
            options["bold_idxs"] = [];
        }
        if (!("bold_color" in options)) {
            options["bold_color"] = d3.rgb(100, 50, 6);
        }
        if (!("color" in options)) {
            options["color"] = d3.rgb(0, 0, 0);
        }
        if (!("circled_vertices" in options)) {
            options["circled_vertices"] = [];
        }
        if (!("dotted_edges" in options)) {
            options["dotted_edges"] = [];
        }
        if (!("dotted_color" in options)) {
            options["dotted_color"] = d3.rgb(0, 0, 0);
        }
        if (!("draw_index" in options)) {
            options["draw_index"] = true;
        }
        
        // Step 1: Draw polygon boundary
        const min_idx = options["min_idx"];
        const r = d/2;
        const w = this.w;
        const N = w.length;
        const dTheta = (2*Math.PI/(N+2));
        let theta = Math.PI/2 + Math.PI/(N+2);
        let Xx = [];
        let Tx = [];
        let Xy = [];
        let Ty = [];
        for (let i = 0; i < N+3; i++) {
            const x = r*Math.cos(theta) + c[0];
            const y = -r*Math.sin(theta) + c[1];
            Xx.push(x);
            Xy.push(y);
            Tx.push(r*1.2*Math.cos(theta) + c[0]);
            Ty.push(-r*1.2*Math.sin(theta) + c[1]);
			g.append("circle")
				.attr("r", 5)
				.attr("fill", options["color"])
				.attr("cx", x).attr("cy", y);
            theta += dTheta;
        }
        let rg = N;
        if (options["draw_index"]) {
            rg = N+2;
        }
        if (options["show_codeword"]) {
            for (let i = 0; i < rg; i++) {
                if (i >= min_idx) {
                    let c = options["color"];
                    if (i in options["bold_idxs"]) {
                        c = options["bold_color"];
                    }
                    let val = i;
                    if (!options["draw_index"]) {
                        val = w[i];
                    }  
                    g.append("text")
                        .attr("x", Tx[i])
                        .attr("y", Ty[i])
                        .attr("text-anchor", "middle")
                        .attr("fill", c)
                        .text(""+val)
                }

            }
        }
            
        

        // Step 2: Draw polygon edges
        // Step 2a: Boundary edges
        for (let i = 0; i < N+2; i++) {
            g.append("line")
            .attr("x1", Xx[i])
            .attr("y1", Xy[i])
            .attr("x2", Xx[i+1])
            .attr("y2", Xy[i+1])
            .attr("stroke", options["color"])
            .attr("stroke-width", 1);
        }
        // Step 2b: Internal edges
        const edges = this.get_edges(min_idx);
        for (let k = 0; k < edges.length; k++) {
            let e = edges[k];
            console.log(e);
            const i = e[0];
            const j = e[1];
            console.log(Xx[i], Xy[i])
            g.append("line")
            .attr("x1", Xx[i])
            .attr("y1", Xy[i])
            .attr("x2", Xx[j])
            .attr("y2", Xy[j])
            .attr("stroke", options["color"])
            .attr("stroke-width", 1);
        }

        // Step 3: Circle any vertices
        /*for idx in options["circled_vertices"]:
            ax.scatter(XTheta[idx, 0], XTheta[idx, 1], s=200, facecolors='none', edgecolors='C3', zorder=101)

        ## Step 4: Draw any dotted edges
        for [i, j] in options["dotted_edges"]:
            ax.plot(X[[i, j], 0], X[[i, j], 1], c=options["dotted_color"], linestyle='--')*/

    }
}


/*
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
        self.codeword_obj = {}
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
        if self.stack_index[d]%2 == 1:
            vals = reversed(vals)
        self.stack_index[d] += 1
        stackorder = np.array(self.stack_index, dtype=int)
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
                self.codeword_obj[tuple(wi)] = self.codewords[-1]
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
*/