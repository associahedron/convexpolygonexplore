/**
 * Sum up all of the elements of the array starting
 * at index i
 * @param {list} arr Array
 * @param {int} idx Start index
 */
function arrsum(arr, idx) {
    let res = 0;
    for (let i = idx; i < arr.length; i++) {
        res += arr[i];
    }
    return res;
}

/**
 * Make a comma separated string out of an array
 * 
 * @param {list} arr Array
 * @param {str} sep Separator
 */
function arrstr(arr, sep) {
    let s = "";
    if (sep == undefined) {
        sep = ""
    }
    for (let i = 0; i < arr.length; i++) {
        s += arr[i];
        if (i < arr.length-1) {
            s += sep;
        }
    }
    return s;
}

/**
 * Search for an array in an array of arrays
 * @param {array of arrays} arr Array of arrays
 * @param {array} inarr Array we're searching for
 */
function arrInArr(arr, inarr) {
    let found = false;
    for (let i = 0; i < arr.length; i++) {
        let equals = (arr[i].length == inarr.length);
        if (equals) {
            for (let k = 0; k < arr[i].length; k++) {
                equals = equals && (arr[i][k] == inarr[k]);
            }
            found = found || equals;
        }
    }
    return found;
}

class Codeword {
    constructor(w) {
        this.w = new Int32Array(w);
    }

    /**
     * Extract the internal edges from the codeword
     * @param {int} min_idx Minimum index at which to start
     * throwing down edges
     * @returns 
     */
    get_edges(min_idx) {
        if (min_idx === undefined) {
            min_idx = 0;
        }
        const w = this.w;
        const N = w.length;
        const visible = new Int32Array(N+2);
        for (let k = 0; k < N+2; k++) {
            visible[k] = 1;
        }
        let i = N-1;
        let edges = [];
        while (i >= min_idx) {
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
            options["bold_color"] = d3.rgb(255, 128, 12);
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
            options["draw_index"] = false;
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
            Tx.push(r*1.25*Math.cos(theta) + c[0]);
            Ty.push(-r*1.25*Math.sin(theta) + c[1]);
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
                        .attr("dy", 5)
                        .attr("fill", c)
                        .text(""+val);
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
            const i = e[0];
            const j = e[1];
            g.append("line")
            .attr("x1", Xx[i])
            .attr("y1", Xy[i])
            .attr("x2", Xx[j])
            .attr("y2", Xy[j])
            .attr("stroke", options["color"])
            .attr("stroke-width", 1);
        }

        // Step 3: Circle any vertices
        for (let i = 0; i < options["circled_vertices"].length; i++) {
            let idx = options["circled_vertices"][i];
			g.append("circle")
            .attr("r", 10)
            .attr("fill", "none")
            .attr("stroke", "red")
            .attr("cx", Tx[idx]).attr("cy", Ty[idx]);
        }

        // Step 4: Draw any dotted edges
        for (let idx = 0; idx < options["dotted_edges"].length; idx++) {
            const e = options["dotted_edges"][idx];
            const i = e[0];
            const j = e[1];
            g.append("line")
            .attr("x1", Xx[i])
            .attr("y1", Xy[i])
            .attr("x2", Xx[j])
            .attr("y2", Xy[j])
            .attr("stroke", options["dotted_color"])
            .attr("stroke-width", 1)
            .style("stroke-dasharray", "10, 5");
        }

    }
}



class Associahedron {
    constructor(n, opts, g) {
        if (opts == undefined) {
            opts = {};
        }
        if (!("diameter" in opts)) {
            opts["diameter"] = 1;
        }
        if (!("g_x_offset" in opts)) {
            opts["g_x_offset"] = opts["diameter"]*2*n;
        }
        if (!("g_y_offset" in opts)) {
            opts["g_y_offset"] = 0;
        }
        this.w = new Int32Array(n);
        this.stack_index = new Int32Array(n);
        this.codewords = [];
        this.last_codeword = null;
        this.codeword_obj = {};
        this.make_stack_rec(this.w, n-1, opts["g_x_offset"], opts["g_y_offset"], opts, g);
    }

    make_stack_rec(w, d, g_x_offset, y_offset, opts, g) {
        const n = w.length;
        const diam = opts["diameter"];
        let dy = -0.1*diam;
        let h = w.length-d-arrsum(w, d+1)+1;
        let x_offset = g_x_offset - 1.5*diam*(n-d-1);
        g.append("text")
        .attr("x", x_offset)
        .attr("y", y_offset)
        .attr("text-anchor", "middle")
        .text("d = " + d + ", h = " + h, dy=0.65*diam);
        let y1 = y_offset+1.4*diam/2;
        let n_items = 0;
        let vals = [];
        if (this.stack_index[d]%2 == 0) {
            vals = (new Int32Array(h)).map((_, idx) => idx);
        }
        else {
            vals = (new Int32Array(h)).map((_, idx) => h-idx-1);            
        }
        this.stack_index[d] += 1;
        const stackorder = new Int32Array(this.stack_index);
        for (let ival = 0; ival < vals.length; ival++) {
            let val = vals[ival];
            let wi = new Int32Array(w);
            wi[d] = val;
            let ni = 0;
            if (d == 1) {
                // Base case
                wi[0] = n-1-arrsum(wi, 1);
                ni = 1;
                let c = new Codeword(wi);
                this.codewords.push({
                    "c":c, 
                    "s":stackorder,
                    "x":x_offset,
                    "y":y_offset
                });
                this.codeword_obj[arrstr(wi, ",")] = this.codewords[this.codewords.length-1];
                
                let dotted_edges = [];
                let circled_vertices = [];
                if (this.codewords.length > 1) {
                    // Indicate quad where flip happened
                    let e1 = c.get_edges();
                    let c2 = this.last_codeword;
                    let e2 = c2.get_edges();
                    for (let k = 0; k < e2.length; k++) {
                        if (!arrInArr(e1, e2[k])) {
                            dotted_edges.push(e2[k]);
                        }
                    }
                    for (let k = 0; k < w.length; k++) {
                        if (c.w[k] != c2.w[k]) {
                            circled_vertices.push(k);
                        }
                    }
                }
                c.draw(g, diam, [x_offset, y_offset+dy], {
                    "bold_idxs":[1],
                    "circled_vertices":circled_vertices,
                    "dotted_edges":dotted_edges
                });
                g.append("text")
                .attr("x", x_offset-2*diam-n*0.07)
                .attr("y", y_offset+dy)
                .attr("text-anchor", "middle")
                .text(arrstr(wi));
                this.last_codeword = c;
            }
            else {
                let c = new Codeword(wi);
                c.draw(g, diam, [x_offset, y_offset+dy], {
                    "min_idx":d,
                    "bold_idxs":[d]
                });
                let s = "";
                for (let k = d; k < wi.length; k++) {
                    s += wi[d];
                }
                g.append("text")
                .attr("x", x_offset)
                .attr("y", y_offset-diam/3.2+dy)
                .attr("text-anchor", "middle")
                .text(s);
                ni = this.make_stack_rec(wi, d-1, g_x_offset, y_offset, opts, g);
            }
            n_items += ni;
            y_offset += diam*1.5*ni
        }
        let y2 = y_offset - 1.4*diam/2;
        let x1 = x_offset - 1.5*diam/2;
        let x2 = x1 + 1.5*diam;
        if (d == 1) {
            x1 -= 0.1*n*diam;
        }

        /*let xbox = [x1, x1, x2, x2, x1]
        let ybox = [y1, y2, y2, y1, y1]
            ax.plot(xbox, ybox, c='k')
            if self.stack_index[d]%2 == 0:
                ax.fill(xbox, ybox, c=[0.9]*3)*/
        return n_items
    }

}
/*
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