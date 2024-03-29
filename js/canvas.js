class Canvas2D {
	constructor() {
		const container = document.getElementById("Canvas2DContainer");
		container.addEventListener("contextmenu", e => e.preventDefault());
		// Fix the width and height up front
		this.width = window.innerWidth * 0.9;
		this.height = window.innerHeight * 0.9;
		document.getElementById("info").width = this.width;
		this.container = container;
		this.canvas = d3.select("#Canvas2DContainer")
		.append("svg")
		.attr("width", this.width)
		.attr("height", this.height)
		.attr("style", "border-style: dotted;");
		this.canvas.on("mousedown", this.mouseDown.bind(this));
		this.container.obj = this;

		// Clear all graph elements if any exist
		this.canvas.selectAll("*").remove();

		let g = this.canvas.append("g").attr("class", "poly1");
		let c = new Codeword([0, 2, 1, 0]);
		c.draw(g, 100, [100, 100]);
	}
	/*updateGodzillaLine() {
		let Ps = this.getGodzillaLine();
		this.godzillaLineCollection.selectAll("line").each(function(){
			let sel = d3.select(this);
			sel.remove();
		});*/
	


	/**
	 * React to a mouse down event by adding a node
	 */
	mouseDown() {
		let point = d3.mouse(d3.event.currentTarget);
	}

	/**
	 * A function which toggles all of the visible elements to show
	 */
	show = function() {
		this.container.style("display", "block");
	}

	/**
	 * A function which toggles all of the visible elements to hide
	 */
	hide = function() {
		this.container.style("display", "none");
	}
}
