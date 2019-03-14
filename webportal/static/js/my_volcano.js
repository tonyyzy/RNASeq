var all_circles;

$( document ).ready(function() {
   p_val = 0.05
   p_threshold = - Math.log10(p_val);
   log2_threshold = 2;
  // console.log(p_threshold)
   p_out = document.getElementById("id_p_threshold_out");
  p_out.innerHTML = p_val;
   lfc = document.getElementById("id_log2_threshold_out");
  lfc.innerHTML = log2_threshold;
});

$("#id_p_threshold").change(function() {
  d3.select("#id_painting_volcano").selectAll("line").remove()
  p_val = $("#id_p_threshold").val()/100
  p_threshold = - Math.log10(p_val);
  var output = document.getElementById("id_p_threshold_out");
  output.innerHTML = p_val;
  // console.log(p_threshold)
  all_circles.attr("fill", function(d){
          if(- Math.log10(d.padj) > p_threshold && (d.log2FoldChange > log2_threshold || d.log2FoldChange <  - log2_threshold)){
            return colours[d.dataset_index];
          } else {
            return "grey";
          };
        })
  hLines()
  vLines()
  return p_threshold
});


$("#id_log2_threshold").change(function() {
  d3.select("#id_painting_volcano").selectAll("line").remove()
  log2_threshold = $("#id_log2_threshold").val()
  console.log(log2_threshold)
  var output = document.getElementById("id_log2_threshold_out");
  output.innerHTML = log2_threshold;
  all_circles.attr("fill", function(d){
          if(- Math.log10(d.padj) > p_threshold && (d.log2FoldChange > log2_threshold || d.log2FoldChange <  - log2_threshold)){
            return colours[d.dataset_index];
          } else {
            return "grey";
          };
        })
  vLines()
  hLines()
  return log2_threshold
});


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


function wf_select(param){
  dataset = []
  result = []
  var all_checkboxes = $('.checkbox_wf')
  for(var i = 0; i < all_checkboxes.length; i++){
    if (all_checkboxes[i].checked) {
      result.push(all_checkboxes[i].value)
     }
  }
  // console.log(result)
  var endpoint = 'wf_one/';
  for(var i = 0; i < result.length; i++){
    endpoint += result[i]
    endpoint += '_'
  }
  console.log(endpoint)
  console.log('hi there old sport')
  d3_run(endpoint)
}

// var volcanoplot = function() {
//   var dataset = [];
//     for (var i = 0; i < all_data.length; i++) {
//       dataset.push(d3.csv("DGE_results_" + all_data[i] + ".csv", function(dataset){
//         dataset.gene_name = dataset[""];
//         dataset.Base_Mean = parseFloat(dataset.baseMean);
//         dataset.log2FoldChange = parseFloat(dataset.log2FoldChange);
//         dataset.pvalue = parseFloat(dataset.pvalue);
//         dataset.padj = parseFloat(dataset.padj);
//         return dataset;
//       }));
// };


function d3_run(endpoint){
  console.log(endpoint)
  d3.json(endpoint).then(function(data) {
    console.log('d3_run called')
    dataset.push(data)
    dataFilter_volcano()
  });
}


function dataFilter_volcano(){
  dataset = dataset[0]
  dataset = dataset.filter(function(d){return ! isNaN(d.padj)});
  // myPlot()
  newPlot(dataset)
}

function newPlot(){
  console.log(dataset);

  var margin = {top: 40, right: 20, bottom: 20, left: 40};
      //Width and height
      w = 800 - margin.right - margin.left;
      h = 800 - margin.top - margin.bottom;

      xScale = d3.scaleLinear()
        .domain([d3.min(dataset, function(d) { return  d.log2FoldChange - 2;}),
                 d3.max(dataset, function(d) { return  d.log2FoldChange + 2;})])
      .range([0, w]);

      d3.select("id_p_threshold")
        .attr("min", "" + d3.min(dataset, function(d) { return d.log2FoldChange;}))
        .attr("max", "" + d3.max(dataset, function(d) { return d.log2FoldChange;}));

      yScale = d3.scaleLinear()
      .domain([d3.min(dataset, function(d) { return - Math.log10(d.padj) - 0.5;}),
               d3.max(dataset, function(d) { return - Math.log10(d.padj) + 0.5;})])
      .range([h, 0]);

      d3.select("id_p_threshold")
      .attr("min", "" + d3.min(dataset, function(d) { return - Math.log10(d.padj);}))
      .attr("max", "" + d3.max(dataset, function(d) { return - Math.log10(d.padj);}));

      xAxis = d3.axisBottom()
      .scale(xScale);
      // .ticks(5)

      yAxis = d3.axisLeft()
      .scale(yScale);
      // d3.select("#id_painting_volcano").selectAll("*").remove();
      svg = d3.select("#id_painting_volcano")
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
      circles()
      hLines()
      vLines()
}

function circles(){
  console.log(dataset)
  svg.selectAll("circle").remove()
        all_circles = svg.selectAll("circle")
        .data(dataset)
        .enter()
        .append("circle")
        .attr("cx", function(d) {
          return xScale(d.log2FoldChange);
        })
        .attr("cy",function(d){
          return yScale(- Math.log10(d.padj));
        })
        .attr("r", 5)
        .attr("fill", function(d){
          if(- Math.log10(d.padj) > p_threshold && (d.log2FoldChange > log2_threshold || d.log2FoldChange <  - log2_threshold)){
            return colours[d.dataset_index];
          } else {
            return "grey";
          };
        })
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
            .attr("y", h + 35 )
            .text("Gene: " + d[""] + " L2fc: " + d.log2FoldChange + " p adj: " + d.padj);
        })
        .on("mouseout", function(d){
          d3.select(this).attr("fill", function(d){
            if(- Math.log10(d.padj) > p_threshold && (d.log2FoldChange > log2_threshold || d.log2FoldChange <  - log2_threshold)){
              return colours[d.dataset_index];
            } else {
              return "grey";
            };
          });
        });
}




function vLines(){
  var v_lines_g = svg.append("g")
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
}

function hLines(){
  var h_lines_g = svg.append("g")
  .attr("class", "lines");

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
