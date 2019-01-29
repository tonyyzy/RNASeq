cwlVersion: v1.0
class: ExpressionTool
requirements:
  InlineJavascriptRequirement: {}
inputs:
  sub_directory: Directory[]
  out_location: string
outputs:
  parenting_out: Directory
expression: |
  ${
      var dirs = [];
      dirs.push({"class": "Directory",
                 "basename": inputs.out_location,
                 "listing": inputs.sub_directory});
      return {"parenting_out": dirs};
      }
