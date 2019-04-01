$( document ).ready(function() {
   p_val = 0.05
   p_threshold = - Math.log10(p_val);
   log2_threshold = 2;
  // console.log(p_threshold)
   p_out = document.getElementsByClassName("p_threshold_out");
  p_out.innerHTML = p_val;
   lfc = document.getElementsByClassName("log2_threshold_out");
  lfc.innerHTML = log2_threshold;
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

  var dataFilter_barplot = function(dataset){
    dataset = dataset[0];
    dataset = dataset.filter(function(d){return d.p_adj != null; });
    dataset = dataset.filter(function(d){return - Math.log10(parseFloat(d.p_adj)) > p_threshold && (parseFloat(d.log2foldchange) > log2_threshold || parseFloat(d.log2foldchange) <  - log2_threshold);});
    var all_values  = [];
    for(var e = 0; e < dataset.length; e++){
      all_values.push(dataset[e].dataset_index)
    }
    var unique = d3.set(all_values).values();
    var values = [];
    for(var e in unique){
      temp_value = dataset.filter(function(d){return d.dataset_index == unique[e];});
      temp_value = temp_value.length
      values.push(temp_value);
    };
    return values
  };

  function wf_select_barplot(param){
    checkbox_selected = []
    var all_checkboxes = $('.checkbox_wf')
    for(var i = 0; i < all_checkboxes.length; i++){
      if (all_checkboxes[i].checked) {
        checkbox_selected.push(all_checkboxes[i].value)
       }
    }
    var endpoint = 'wf_data_mod/';
    for(var i = 0; i < checkbox_selected.length; i++){
      endpoint += checkbox_selected[i]
      endpoint += '_'
    }
    console.log(endpoint)
    d3_run_barplot(endpoint)
  }



function d3_run_barplot(endpoint){
  console.log(endpoint)
  d3.json(endpoint).then(function(data) {
    console.log('d3_run called')
    dataset.push(data)
    newPlot_Barplot(dataset)
  });
}

var newPlot_Barplot = function(){

    d3.select("#id_painting_barplot").selectAll("*").remove();

    dataset = dataFilter_barplot(dataset);

    // d3.select("#painting").selectAll("*").remove();

    var width = d3.select('#id_plotting_column_barplot').node().getBoundingClientRect().width;
    var height = d3.select('#id_plotting_column_barplot').node().getBoundingClientRect().height;
    console.log(width)
    console.log(height)

    var margin = {top: 40, right: 20, bottom: 20, left: 60};
    //Width and height
    w = width - margin.right - margin.left;
    h = height - margin.top - margin.bottom;


    var xScale = d3.scaleBand()
    .domain(d3.range(dataset.length))
    .rangeRound([0, w]) // <-- Also enables rounding
    .paddingInner(0.05);

    var yScale = d3.scaleLinear()
    .domain([0, d3.max(dataset, function(d) { return d;}) + 1])
    .range([h, 0]);

    var yAxis = d3.axisLeft()
    .scale(yScale);


    var svg = d3.select("#id_painting_barplot")
    //.attr("width", w + margin.right + margin.left) // <-- Here
    //.attr("height", h + margin.top + margin.bottom)
    .attr("width", "100%") // <-- Here
    .attr("height", "100%")
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.bottom + ")");// <-- and here!
    svg.append("g")
    .attr("class", "axis") //Assign "axis" class
    .call(yAxis);

    // text label for the y axis
    svg.append("text")
      .attr("transform", "rotate(-90)")
      .attr("y", 0 - margin.right - 40)
      .attr("x", 0 - (height / 2))
      .attr("dy", "1em")
      .attr("font-size", "1.5rem")
      .style("text-anchor", "middle")
      .text("Number of significant genes");

    var barPadding = 10;

    svg.selectAll("rect")
    .data(dataset)
    .enter()
    .append("rect")
    .attr("x", function(d, i) {
      return xScale(i); // <-- Set x values
    })
    .attr("y", function(d){return yScale(d);})
    .attr("width", xScale.bandwidth()/2)
    .attr("height", function(d){return h - yScale(d);})
    .attr("fill", function(d, i){return colours[i];})
    .on("mouseover", function(d){
      d3.select("#label").remove();
      svg.append("text")
        .attr("id", "label")
        .attr("x", w / 2 )
        .attr("y", h + 25)
        .text("Significant genes: " + d);
    })
    .on("mouseout", function(){
      d3.select("#label").remove();
    });
  };
