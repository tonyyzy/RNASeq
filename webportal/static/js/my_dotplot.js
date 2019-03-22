var dataset = []
var all_data = [1,2]

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

var get_data = function(){
  dataset = []
  for(var i = 0; i < all_data.length; i++){
     dataset.push(d3.csv("gene_count_matrix_" + all_data[i] + ".csv"));
  };
  Promise.all(dataset).then(function(dataset){
  for(var i = 0; i < dataset.length; i++){
    for(var e = 0; e < dataset[i].length; e++){
      dataset[i][e].dataset_index = i;
    };
  };
  dataset = [].concat.apply([], dataset);
  console.log(dataset);

  var value_selected = document.getElementById('selector').value;
  dataset = dataset.filter(function(d){return d.gene_id == value_selected});
  console.log(dataset);
  var fixed_dataset = [];

  for(var i = 0; i < dataset.length; i++){
    var keys = Object.keys(dataset[i]);
    for(var e = 1; e <  keys.length - 1; e++){
      fixed_dataset.push({gene_id: dataset[i].gene_id,
                          value: dataset[i][keys[e]],
                          condition: keys[e],
                          dataset_index: dataset[i].dataset_index});
    };
  };
  console.log(fixed_dataset);

      var margin = {top: 40, right: 20, bottom: 20, left: 40};
      //Width and height
      w = 800 - margin.right - margin.left;
      h = 800 - margin.top - margin.bottom;

      var yScale = d3.scaleLinear()
        .domain([d3.min(fixed_dataset, function(d) { return  d.value;}),
                 d3.max(fixed_dataset, function(d) { return  d.value;})])
      .range([0, h]);


      var xScale = d3.scaleBand()
        .domain(d3.range(fixed_dataset.length))
        .rangeRound([0, w]) // <-- Also enables rounding
        .paddingInner(0.05);

      xAxis = d3.axisBottom()
      .scale(xScale);
      // .ticks(5)

      yAxis = d3.axisLeft()
      .scale(yScale);

      d3.select("#dotplot").selectAll("*").remove();
      svg = d3.select("#dotplot")
      .append("svg")
      .attr("width", w + margin.right + margin.left)
      .attr("height", h + margin.top + margin.bottom)
      .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.bottom + ")");// <-- and here!
      svg.append("g")
      .attr("class", "axis") //Assign "axis" class
      .call(yAxis)
      svg.append("g")
      .attr("class", "axis") //Assign "axis" class
      .attr("transform", "translate(0," + h + ")") // <-- and here!
      .call(xAxis);

      var all_circles = svg.selectAll("circle")
        .data(fixed_dataset)
        .enter()
        .append("circle")
        .attr("cx",function(d, i) {
          return xScale(i); // <-- Set x values
        })
        .attr("cy",function(d){
          return yScale(parseFloat(d.value));
        })
        .attr("r", 5)
        .attr("fill", function(d){
            return colours[d.dataset_index];
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
            .text("Gene: " + d.gene_id + " condition: " + d.condition);
        })
        .on("mouseout", function(d){
          d3.select(this).attr("fill",function(d){
            return colours[d.dataset_index];
        });
        });
});
}

var dataFilter_dotplot = function(){

}

var plotting = function(){
  console.log("test");
}
