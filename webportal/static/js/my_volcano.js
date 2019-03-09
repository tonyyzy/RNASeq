
var p_threshold = - Math.log10(0.05);
var log2_threshold = 1;

console.log('my_volcano js active')
var in_file = 'head_DGE_1.csv'

// var in_file = 'DGE_results_1_head.csv'



dataset = []
d3.csv(in_file, function(d) {
  return {
    gene_name : d[""],
    Base_Mean : +d.Base_Mean,
    log2FoldChange : +d.log2FoldChange,
    pvalue : +d.pvalue,
    padj : +d.padj
  };
}).then(function(data) {
  console.log(data[0]);
  dataset.push(data)
  arr_delay()
});


function dataFilter(){
  console.log('called')
  dataset = dataset[0]
  dataset = dataset.filter(function(d){return ! isNaN(d.padj)});
  myPlot()
}

function arr_delay(){
  console.log('delay')
  setTimeout(function() { dataFilter(); }, 1000);
}

function myPlot(){

  var margin = {top: 40, right: 20, bottom: 20, left: 40};
  //Width and height
  var w = 800 - margin.right - margin.left;
  var h = 800 - margin.top - margin.bottom;

  var xScale = d3.scaleLinear()
  .domain([d3.min(dataset, function(d) { return d.log2FoldChange - 0.5;}),
           d3.max(dataset, function(d) { return d.log2FoldChange + 0.5;})])
  .range([0, w]);

  var yScale = d3.scaleLinear()
  .domain([d3.min(dataset, function(d) { return - Math.log10(d.padj) - 0.5;}),
           d3.max(dataset, function(d) { return - Math.log10(d.padj) + 0.5;})])
  .range([h, 0]);

  var xAxis = d3.axisBottom()
  .scale(xScale);

  var yAxis = d3.axisLeft()
  .scale(yScale);

  var svg = d3.select("#painting")
  .append("svg")
  //.attr("width", w + margin.right + margin.left) // <-- Here
  //.attr("height", h + margin.top + margin.bottom)
  .attr("width", "100%") // <-- Here
  .attr("height", "100%")
  .append("g")
  .attr("transform", "translate(" + margin.left + "," + margin.bottom + ")");// <-- and here!
  svg.append("g")
  .attr("class", "axis") //Assign "axis" class
  .attr("transform", "translate(0," + h + ")") // <-- and here!
  .call(xAxis);
  svg.append("g")
  .attr("class", "axis") //Assign "axis" class
  .call(yAxis)

  var circles = svg.selectAll("circle")
  .data(dataset)
  .enter()
  .append("circle")
  .on("click", function() {
  //Do something mundane and annoying on click
  alert("Hey, don't click that!");
  })
  .on("mouseover", function(d) {
    d3.select("#label").remove();
    d3.select(this)
    .attr("fill", "orange");
    svg.append("text")
      .attr("id", "label")
      .attr("y", h + margin.bottom + margin.top/4)
      .text("Gene: " + d.gene_name + " L2fc: " + d.log2FoldChange + " p adj: " + d.padj);
  })
  .on("mouseout", function(d){
    d3.select(this).attr("fill", function(d){
      if(- Math.log10(d.padj) > p_threshold && (d.log2FoldChange > log2_threshold || d.log2FoldChange <  - log2_threshold)){
        return colours[d.dataset];
      } else {
        return "grey";
      };
    });
  });

  circles.attr("cx", function(d) {
  return xScale(d.log2FoldChange);
  })
  .attr("cy",function(d){
    return yScale(- Math.log10(d.padj));
  })
  .attr("r", 5)
  .attr("fill", function(d){
    if(- Math.log10(d.padj) > p_threshold && (d.log2FoldChange > log2_threshold || d.log2FoldChange <  - log2_threshold)){
      return colours[d.dataset];
    } else {
      return "grey";
    };
  });

  var v_lines_g = svg.append("g")
  .attr("class", "lines");

  var h_lines_g = svg.append("g")
  .attr("class", "lines");

  v_lines = [- log2_threshold, log2_threshold]
  var v_lines = v_lines_g.selectAll("line")
  .data(v_lines)
  .enter()
  .append("line")
  .attr("x1", function(d){
    return xScale(d);})
  .attr("y1", "0")
  .attr("x2", function(d){
    return xScale(d);})
  .attr("y2", "" + h)
  .attr("stroke", "red")
  .attr("stroke-width", "2");

  h_line = [p_threshold]
  var h_lines = h_lines_g.selectAll("line")
  .data(h_line)
  .enter()
  .append("line")
  .attr("x1", "0")
  .attr("y1", function(d){
    return yScale(d);})
  .attr("x2", "" + w)
  .attr("y2",function(d){ return yScale(d);})
  .attr("stroke", "red")
  .attr("stroke-width", "2");

}
