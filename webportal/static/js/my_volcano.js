
var p_threshold = - Math.log10(0.005);
var log2_threshold = 2;
console.log(p_threshold)


$("#id_p_threshold").change(function() {
  var p_val = $("#id_p_threshold").val()/100
  var p_thresh = - Math.log10(p_val);
  console.log(p_threshold)
});


// var initial = 10
// function one(initial){
//   intermediate = initial * 2
//   return intermediate
// }
//
// function two(param){
//   old = param /2
//   return old
// }

// note to morning self, calling the var param wrecks shit!
colours = [
  "#000000",
  "#009292",
  "#FF6DB6",
  "#FFB6DB",
  "#490092",
  "#006DDB",
  "#B66DFF",
  "#B6DBFF",
  "#6DB6FF",
  "#920000",
  "#924900",
  "#DB6D00",
  "#24FE23",
  "#074751",
  "#FFFF6D"];


dataset = []
function wf_select(param){
  if(param == 1 ){
    console.log('wf_1 chosen')
    var endpoint = 'wf_one/1'
  }
  if(param == 2 ){
    console.log('wf_2 chosen')
    var endpoint = 'wf_one/2'
  }
  function d3_run(){
    d3.json(endpoint).then(function(data) {
      console.log('d3_run called')
      console.log(p_threshold)
      dataset.push(data)
      dataFilter()
    });
  }
  console.log(p_threshold)
  d3_run()
}

function dataFilter(){
  dataset = dataset[0]
  dataset = dataset.filter(function(d){return ! isNaN(d.padj)});
  // myPlot()
  newPlot(dataset)
}

function newPlot(){
  var margin = {top: 40, right: 20, bottom: 20, left: 40};
      //Width and height
      var w = 800 - margin.right - margin.left;
      var h = 800 - margin.top - margin.bottom;

      var xScale = d3.scaleLinear()
      // .domain([d3.min(dataset, function(d) { return d.log2FoldChange - 0.5;}),
               // d3.max(dataset, function(d) { return d.log2FoldChange + 0.5;})])

      .domain([d3.min(dataset, function(d) { return - d.log2FoldChange - 5;}),
               d3.max(dataset, function(d) { return - d.log2FoldChange + 5;})])
      .range([0, w]);

      var yScale = d3.scaleLinear()
      .domain([d3.min(dataset, function(d) { return - Math.log10(d.padj) - 0.5;}),
               d3.max(dataset, function(d) { return - Math.log10(d.padj) + 0.5;})])
      .range([h, 0]);

      var xAxis = d3.axisBottom()
      .scale(xScale);
      // .ticks(5)

      var yAxis = d3.axisLeft()
      .scale(yScale);

      var svg = d3.select("#painting")
      .append("svg")
      .attr("width", w + margin.right + margin.left) // <-- Here
      .attr("height", h + margin.top + margin.bottom)
      // .attr("width", "100%") // <-- Here
      // .attr("height", "100%")
      .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.bottom + ")");// <-- and here!
      svg.append("g")
      .attr("class", "axis") //Assign "axis" class
      .attr("transform", "translate(0," + h + ")") // <-- and here!
      .call(xAxis);
      svg.append("g")
      .attr("class", "axis") //Assign "axis" class
      .call(yAxis)

      var circles = svg.selectAll('circle')
          .data(dataset)
          .enter().append('circle')
        .attr('cx',function (d) { return xScale(d.log2FoldChange) })
        .attr('cy',function (d) { return yScale(- Math.log10(d.padj)) })
        .attr('r','5')
        // .attr('cy',function (d) { return yScale(d.aror) })
        // .attr('stroke','black')
        // .attr('stroke-width',1)
        // .attr('fill',function (d,i) { return colorScale(i) })

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
