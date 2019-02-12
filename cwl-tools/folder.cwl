#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: ExpressionTool
requirements:
  InlineJavascriptRequirement: {}

inputs:
  item: Any
  name: string

outputs:
  out: Directory

expression: "${
    if (inputs.item.class == 'File'){
        var arr = [inputs.item];
        }
    else {
        var arr = inputs.item;
    }
    return {
        'out': {
            'class': 'Directory',
            'basename': inputs.name,
            'listing': arr
        }
    }
}"