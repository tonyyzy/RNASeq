var dataset;
var metadata;
var data_sorted;
var layout;
function sort_data_dot_plot() {
  var ids = dataset.gene_id
  var newData = dataset.filter(function(obj){return obj.gene_id === document.getElementById("gene").value})[0];
  var tmp = [];
  var i;
  for(i = 1; i <= metadata.name.length; i++){
    var tmp2 = metadata.name[i-1];
    tmp.push(newData[tmp2]);
  };
  var trace1 = {
    type: 'scatter',
    y: tmp,
    x: metadata.condition,
    mode: 'markers',
    name: 'condition',
    marker: {
      color: 'rgba(0, 0, 0, 0.95)',
      line: {
        color: 'rgba(0, 0, 0, 1.0)',
        width: 1,
      },
      symbol: 'circle',
      size: 16
    }
  };

  var data = [trace1];

  var layout = {
    title: {
      text: newData.gene_id,
      font: {
        family: document.getElementById("font").value,
        size: parseFloat(document.getElementById("text_size").value),
        color: 'rgba(0, 0, 0, 1.0)'
      }
    },
    yaxis: {
      showgrid: false,
      showline: true,
      linecolor: 'rgb(102, 102, 102)',
      titlefont: {
        font: {
          color: 'rgb(204, 204, 204)'
        }
      },
      tickfont: {
        font: {
          color: 'rgb(102, 102, 102)'
        }
      },
      autotick: false,
      dtick: 10,
      ticks: 'outside',
      tickcolor: 'rgb(102, 102, 102)',
      title: {
        text: 'gene count',
        font: {
          family: document.getElementById("font").value,
          size: parseFloat(document.getElementById("text_size").value),
          color: 'rgba(0, 0, 0, 1.0)'
        }
      },
    },
    xaxis: {
       type: 'category',
      showline: true,
      title: {
        text: 'condition',
        font: {
          family: document.getElementById("font").value,
          size: parseFloat(document.getElementById("text_size").value),
          color: 'rgba(0, 0, 0, 1.0)'
        }
      }
    },
    margin: {
      l: 140,
      r: 40,
      b: 50,
      t: 80
    },
    legend: {
      font: {
        size: 10,
      },
      yanchor: 'middle',
      xanchor: 'right'
    },
    width: 600,
    height: 600,
    paper_bgcolor: 'rgb(255, 255, 255)',
    plot_bgcolor: 'rgb(255, 255, 255)',
    hovermode: 'closest'
  };
  return [data, layout]
}
function sort_data_volcano() {
  data = dataset.filter(function(d){
    if(isNaN(d.padj)){
        return false;
    }
    return true;
  });
  function log10(val) {
    return Math.log(val) / Math.LN10;
  }
  var name = [],
  log2FoldChange = [],
  log_padj = [],
  sig_log2FoldChange = [],
  sig_log_padj = [];
  data.map(function(d) {
    name.push(d.name);
    if(d.padj <= parseFloat(document.getElementById("pvalue").value) && Math.abs(d.log2FoldChange) > parseFloat(document.getElementById("foldchange").value) ){
      sig_log2FoldChange.push(d.log2FoldChange);
      sig_log_padj.push(log10(d.padj)*-1);
    } else {
      log2FoldChange.push(d.log2FoldChange);
      log_padj.push(log10(d.padj)*-1);
    }
  })
  function getMaxOfArray(numArray) {
    return Math.max.apply(null, numArray);
  }
  var max = getMaxOfArray(sig_log_padj)
  var trace1 = {
    x: sig_log2FoldChange,
    y: sig_log_padj,
    mode: 'markers',
    name: 'sig',
    marker: {
      color: 'rgba(0, 0, 0, 0.95)',
      line: {
        color: document.getElementById("sig_colour").value,
        width: 1,
      },
      symbol: 'circle',
      size: 1
    }
  };
  if(document.getElementById("nonsig").checked == true){
    var trace2 = {
      x: log2FoldChange,
      y: log_padj,
      mode: 'markers',
      name: 'unsig',
      marker: {
        color: 'rgba(0, 0, 0, 0.95)',
        line: {
          color: 'rgba(0, 0, 0, 1.0)',
          width: 1,
        },
        symbol: 'circle',
        size: 1
      }
    };

    var data = [trace1, trace2];
  } else {
    var data = [trace1];
  }

  var layout = {
    title: {
      text: document.getElementById("title").value,
      font: {
        family: document.getElementById("font").value,
        size: parseFloat(document.getElementById("text_size").value),
        color: 'rgba(0, 0, 0, 1.0)'
      }
    },
    yaxis: {
      showline: true,
      showgrid: false,
      title: {
        text: 'negative log10 p_adj',
        font: {
          family: document.getElementById("font").value,
          size: parseFloat(document.getElementById("text_size").value),
          color: 'rgba(0, 0, 0, 1.0)'
        }
      },
    },
    xaxis: {
      showline: false,
      showgrid: false,
      title: {
        text: 'log2FoldChange',
        font: {
          family: document.getElementById("font").value,
          size: parseFloat(document.getElementById("text_size").value),
          color: 'rgba(0, 0, 0, 1.0)'
        }
      }
    },
    margin: {
      l: 140,
      r: 40,
      b: 50,
      t: 80
    },
    legend: {
      font: {
        size: 10,
      },
      yanchor: 'middle',
      xanchor: 'right'
    },
    width: 600,
    height: 600,
    shapes: [
      {
        type: 'line',
        x0: parseFloat(document.getElementById("foldchange").value),
        y0: 0,
        x1: parseFloat(document.getElementById("foldchange").value),
        y1: max,
        line: {
          color: document.getElementById("line_colour").value,
          width: 2
        }
      },
      {
        type: 'line',
        x0: -1*parseFloat(document.getElementById("foldchange").value),
        y0: 0,
        x1: -1*parseFloat(document.getElementById("foldchange").value),
        y1: max,
        line: {
          color: document.getElementById("line_colour").value,
          width: 2
        }
      }
    ]
  };
  return [data, layout]
}
function read_data(file) {
  var e = document.getElementById(file);
  var strUser = e.options[e.selectedIndex].value;
  d3.csv(strUser).then(function(data){
      dataset = data;
  });
}
d3.csv("metadata.csv").then(function(data){
    metadata = data;
    var name = [];
    var condition = [];
    var libtype = [];
    metadata.map(function(d) {
      name.push(d.name);
      condition.push(d.condition);
      libtype.push(d.libtype);
    });
    metadata = {name:name, condition:condition, libtype:libtype};
});
function create_plot(func){
  out = func();
  Plotly.newPlot('main', out[0], out[1], {displayModeBar: false});
}
function download_plot(func){
  var img_jpg= d3.select('#jpg-export');
  out = func();
  Plotly.plot('divDownload', out[0], out[1]).then(function(gd){
    Plotly.downloadImage(gd,{format:'jpeg',height:300,width:300, filename:'newPlot'})
  });
}
