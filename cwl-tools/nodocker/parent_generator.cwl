cwlVersion: v1.0
class: ExpressionTool
requirements:
  InlineJavascriptRequirement: {}
inputs:
  sub_directory: Directory[]
  name: string
outputs:
  parenting_out: Directory[]
expression: |
  ${
      var dirs = [];
      dirs.push({"class": "Directory",
                 "basename": "sample_" + inputs.name,
                 "listing": inputs.sub_directory});
      return {"parenting_out": dirs};
      }
