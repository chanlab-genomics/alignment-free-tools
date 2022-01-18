
//Constants for the SVG
var margin = {top: -5, right: -5, bottom: -5, left: -5},
    width = 1600 - margin.left - margin.right,
    height = 1200 - margin.top - margin.bottom;
    
//Set up the colour scale
var color = d3.scale.category20();

//Set up the force layout
var force = d3.layout.force()
    .charge(-120)
    .linkDistance(30)
    .size([width, height]);
    
var node_drag = d3.behavior.drag()
.on("dragstart", dragstart)
.on("drag", dragmove)
.on("dragend", dragend);
function dragstart(d, i) {
  force.stop() // stops the force auto positioning before you start dragging
}
function dragmove(d, i) {
  d.px += d3.event.dx;
  d.py += d3.event.dy;
  d.x += d3.event.dx;
  d.y += d3.event.dy;
}
function dragend(d, i) {
  d.fixed = true; // of course set the node to fixed so the force doesn't include the node in its auto positioning stuff
  force.resume();
}
function releasenode(d) {
  d.fixed = false; // of course set the node to fixed so the force doesn't include the node in its auto positioning stuff
  //force.resume();
}
var zoom = d3.behavior.zoom()
    .scaleExtent([1, 10])
    .on("zoom", zoomed);
    
//Append a SVG to the body of the html page. Assign this SVG as an object to svg
var svg = d3.select("body").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .attr("transform", "translate(" + margin.left + "," + margin.right + ")");
    //.call(zoom);

//Read the data from the mis element 
var mis = document.getElementById('mis').innerHTML;
graph = JSON.parse(mis);
graphRec=JSON.parse(JSON.stringify(graph)); //Add this line

//Creates the graph data structure out of the json data
force.nodes(graph.nodes)
    .links(graph.links)
    .start();

//Create all the line svgs but without locations yet
var link = svg.selectAll(".link")
    .data(graph.links)
    .enter().append("line")
    .attr("class", "link")
    .style("stroke-width", function (d) {
    //return Math.sqrt(d.value);
    return 2;
});

//Do the same with the circles for the nodes - no 
var node = svg.selectAll(".node")
    .data(graph.nodes)
    .enter().append("circle")
    .attr("class", "node")
    .attr("r", 5)
    .style("fill", function (d) {
    return color(d.group);
})
    //.call(force.drag)
    .call(node_drag)
    .on('click', releasenode)
    .on('dblclick', connectedNodes);
    

//Legend    
var legend = svg.selectAll(".legend")
    .data(color.domain())
    .enter().append("g")
    .attr("class", "legend")
    .attr("transform", function(d, i) { return "translate(0," + i * 20 + ")"; });

legend.append("rect")
    .attr("x", width - 30)
    .attr("y", 20)
    .attr("width", 18)
    .attr("height", 18)
    .style("fill", color);

legend.append("text")
    .attr("x", width - 35)
    .attr("y", 29)
    .attr("dy", ".35em")
    .style("text-anchor", "end")
    .text(function(d) { return d; });

    
node.append("title")
          .text(function(d) { return d.name+","+d.group; });
          
//Now we are giving the SVGs co-ordinates - the force layout is generating the co-ordinates which this code is using to update the attributes of the SVG elements
force.on("tick", function () {
    link.attr("x1", function (d) {
        return d.source.x;
    })
        .attr("y1", function (d) {
        return d.source.y;
    })
        .attr("x2", function (d) {
        return d.target.x;
    })
        .attr("y2", function (d) {
        return d.target.y;
    });

    node.attr("cx", function (d) {
        return d.x;
    })
        .attr("cy", function (d) {
        return d.y;
    });
    node.each(collide(0.5)); //Added 
});

/*
var fisheye = d3.fisheye.circular()
      .radius(120);


svg.on("mousemove", function() {
      //force.stop();
    
      fisheye.focus(d3.mouse(this));

      d3.selectAll("circle").each(function(d) { d.fisheye = fisheye(d); })
          .attr("cx", function(d) { return d.fisheye.x; })
          .attr("cy", function(d) { return d.fisheye.y; })
          .attr("r", function(d) { return d.fisheye.z * 8; });

      link.attr("x1", function(d) { return d.source.fisheye.x; })
          .attr("y1", function(d) { return d.source.fisheye.y; })
          .attr("x2", function(d) { return d.target.fisheye.x; })
          .attr("y2", function(d) { return d.target.fisheye.y; });

    
    });
*/

svg.on("contextmenu", function(){
	force.stop();
});  

var padding = 1, // separation between circles
    radius=8;

function collide(alpha) {
  var quadtree = d3.geom.quadtree(graph.nodes);
  return function(d) {
    var rb = 2*radius + padding,
        nx1 = d.x - rb,
        nx2 = d.x + rb,
        ny1 = d.y - rb,
        ny2 = d.y + rb;
    quadtree.visit(function(quad, x1, y1, x2, y2) {
      if (quad.point && (quad.point !== d)) {
        var x = d.x - quad.point.x,
            y = d.y - quad.point.y,
            l = Math.sqrt(x * x + y * y);
          if (l < rb) {
          l = (l - rb) / l * alpha;
          d.x -= x *= l;
          d.y -= y *= l;
          quad.point.x += x;
          quad.point.y += y;
        }
      }
      return x1 > nx2 || x2 < nx1 || y1 > ny2 || y2 < ny1;
    });
  };
}

//Toggle stores whether the highlighting is on
var toggle = 0;
//Create an array logging what is connected to what
var linkedByIndex = {};
for (i = 0; i < graph.nodes.length; i++) {
    linkedByIndex[i + "," + i] = 1;
};
graph.links.forEach(function (d) {
    linkedByIndex[d.source.index + "," + d.target.index] = 1;
});

//This function looks up whether a pair are neighbours
function neighboring(a, b) {
    return linkedByIndex[a.index + "," + b.index];
}

function connectedNodes() {
    if (toggle == 0) {
        //Reduce the opacity of all but the neighbouring nodes
        d = d3.select(this).node().__data__;
        node.style("opacity", function (o) {
            return neighboring(d, o) | neighboring(o, d) ? 1 : 0.1;
        });
        link.style("opacity", function (o) {
            return d.index==o.source.index | d.index==o.target.index ? 1 : 0.1;
        });
        //Reduce the op
        toggle = 1;
    } else {
        //Put them back to opacity=1
        node.style("opacity", 1);
        link.style("opacity", 1);
        toggle = 0;
    }
}

var optArray = [];
for (var i = 0; i < graph.nodes.length - 1; i++) {
    optArray.push(graph.nodes[i].name);
}
optArray = optArray.sort();
$(function () {
    $("#search").autocomplete({
        source: optArray
    });
});

function searchNode() {
    //find the node
    var selectedVal = document.getElementById('search').value;
    var node = svg.selectAll(".node");
    if (selectedVal == "none") {
        node.style("stroke", "white").style("stroke-width", "1");
    } else {
        var selected = node.filter(function (d, i) {
        	if (!d.name.includes(selectedVal)){
            	return d.name != selectedVal;
            }
        });
        selected.style("opacity", "0");
        var link = svg.selectAll(".link")
        link.style("opacity", "0");
        d3.selectAll(".node, .link").transition()
            .duration(5000)
            .style("opacity", 1);
    }
}


//adjust threshold
function threshold(thresh) {
    graph.links.splice(0, graph.links.length);

		for (var i = 0; i < graphRec.links.length; i++) {
			if (graphRec.links[i].value > thresh) {graph.links.push(graphRec.links[i]);}
		}
    restart();
}

function zoomed() {
  svg.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
}

//Restart the visualisation after any node and link changes
function restart() {
	
	link = link.data(graph.links);
	link.exit().remove();
	link.enter().insert("line", ".node").attr("class", "link").style("stroke-width", 2);
	node = node.data(graph.nodes);
	node.enter().insert("circle", ".cursor").attr("class", "node").attr("r", 5).call(force.drag);
	force.start();
}