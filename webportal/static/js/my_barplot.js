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
    dataset = dataset.filter(function(d){return ! isNaN(d.padj);});
    dataset = dataset.filter(function(d){return - Math.log10(d.padj) > p_threshold && (d.log2FoldChange > log2_threshold || d.log2FoldChange <  - log2_threshold);});
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

    dataset = dataFilter_barplot(dataset)

    d3.select("#painting").selectAll("*").remove();

    var margin = {top: 40, right: 20, bottom: 20, left: 40};
    //Width and height
    var w = 800 - margin.right - margin.left;
    var h = 800 - margin.top - margin.bottom;

    var xScale = d3.scaleBand()
    .domain(d3.range(dataset.length))
    .rangeRound([0, w]) // <-- Also enables rounding
    .paddingInner(0.05);

    var yScale = d3.scaleLinear()
    .domain([0, d3.max(dataset, function(d) { return d;}) + 10])
    .range([h, 0]);

    var yAxis = d3.axisLeft()
    .scale(yScale);

    // d3.select("#painting_barplot").selectAll("*").remove();

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

    var barPadding = 10;

    svg.selectAll("rect")
    .data(dataset)
    .enter()
    .append("rect")
    .attr("x", function(d, i) {
      return xScale(i); // <-- Set x values
    })
    .attr("y", function(d){return yScale(d);})
    .attr("width", xScale.bandwidth())
    .attr("height", function(d){return h - yScale(d);})
    .attr("fill", function(d, i){return colours[i];})
    .on("mouseover", function(d){
    svg.append("text")
        .attr("id", "label")
        .attr("x", 400)
        .attr("y", h + 35)
        .text("Significant genes: " + d);
    })
    .on("mouseout", function(){
      d3.select("#label").remove();
    });
  });
}
