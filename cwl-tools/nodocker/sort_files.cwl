cwlVersion: v1.0
class: ExpressionTool
requirements:
  InlineJavascriptRequirement: {}
inputs:
  files: File[]
  input_name: string
outputs:
  foldering_out: Directory[]
expression: |
  ${
      var dirs = [];
      var samples = [];
      for (var i = 0; i < inputs.files.length; i++) {
          var file = inputs.files[i];
          samples.push(file);
        }
        for (var e in samples) {
          dirs.push({"class": "Directory",
                  "basename": inputs.input_name,
                  "listing": [samples[e]]});
                }
      return {"foldering_out": dirs};
      }
