import yaml

cwl_workflow = {"cwl": {"cwlVersion": "v1.0",
                    "class": "Workflow",
                    "requirements":{
                        "ScatterFeatureRequirement": {},
                        "MultipleInputFeatureRequirement": {},
                        "StepInputExpressionRequirement": {},
                        "InlineJavascriptRequirement": {}},
                    "inputs": {"threads": "string"},
                    "outputs": {},
                    "steps": {}
                    },
                "yml": {
                    "threads": 2
                }}


def star_readmap(input_files, yaml):
    for i in range(len(input_files)):
        yaml["cwl"]["inputs"]["genomeDir"] = "Directory"
        yaml["cwl"]["inputs"]["subject_name{0}".format(i + 1)] = "File[]"
        yaml["cwl"]["inputs"]["fastq{0}".format(i + 1)] = "string"
        yaml["cwl"]["outputs"]["star_readmap_1_out"] = {
        "type": "Directory",
        "outputSource": "star_readmap_{0}/star_read_out".format(i + 1)}
        yaml["cwl"]["steps"]["star_readmap_{0}".format(i + 1)] = {
        "run": "../cwl-tools/docker/STAR_readmap.cwl",
        "in": {
            "threads": "threads",
            "genomeDir": "genomeDir",
            "readFilesIn": "fastq{0}".format(i + 1),
            "outFileNamePrefix": "subject_name{0}".format(i + 1)},
        "out": ["sam_output", "star_read_out"]}
        yaml["yml"]["genomeDir"] = {
            "class": "Directory",
            "path": "../tests/GenomeIndex"}
        yaml["yml"]["subject_name{0}".format(i + 1)] = "test{0}".format(i + 1)
        yaml["yml"]["fastq{0}".format(i + 1)] = [{"class": "File", "path": "path"} for e in range(len(input_files))]

    yaml["cwl"]["steps"]["star_folder"] = {
    "run": {
        "class": "ExpressionTool",
        "requirements": {
            "InlineJavascriptRequirement": {}},
        "inputs": dict([("dir{0}".format(i + 1), "Directory") for i in range(len(input_files))]),
        "outputs": {
            "out": "Directory"},
        "expression": "|${{return{{\"out\":{{\"class\":\"Directory\",\"basename\":\"star\",\"listing\":[{}]}};}}}}".format(",".join(["inputs.dir{0}".format(i + 1) for i in range(len(input_files))]))},
    "in": dict([("dir{0}".format(i + 1), "star_readmap_{0}/star_read_out".format(i + 1)) for i in range(len(input_files))]),
    "out": ["out"]
    }

    return yaml

def samtools(input_files, yaml):
    for i in range(len(input_files)):
        yaml["cwl"]["steps"]["samtools_{0}".format(i + 1)] = {
        "run": "../cwl-tools/docker/samtools.cwl",
        "in": {
            "samfile": "star_readmap_{0}/sam_output".format(i + 1),
            "threads": "threads",
            "outfilename": {
                "source": ["subject_name{0}".format(i + 1)],
                "valueFrom": "$(self + \".bam\")"}},
        "out": ["samtools_out"]}

    yaml["cwl"]["steps"]["samtools_folder"] = {
    "run": {
        "class": "ExpressionTool",
        "requirements": {
            "InlineJavascriptRequirement": {}},
        "inputs": dict([("file{0}".format(i + 1), "File") for i in range(len(input_files))]),
        "outputs": {
            "out": "Directory"},
        "expression": "|${{return{{\"out\":{{\"class\":\"Directory\",\"basename\":\"samtools\",\"listing\":[{}]}};}}}}".format(",".join(["inputs.file{0}".format(i + 1) for i in range(len(input_files))]))},
    "in": dict([("file{0}".format(i + 1), "samtools_{0}/samtools_out".format(i + 1)) for i in range(len(input_files))]),
    "out": ["out"]
    }

    return yaml


def prepDE(input_files, yaml):
    yaml["cwl"]["inputs"] = {"prepDE_script": "File"}
    yaml["cwl"]["outputs"] = {"prepDE_out":{
            "type": "Directory",
            "outputSource": "prepDE_folder/out"}}
    yaml["cwl"]["steps"]["prepDE"] = {
        "run": "../cwl-tools/docker/prepDE.cwl",
        "in": {"program": "prepDE_script", "gtfs": ["stringtie_{0}/stringtie_out".format(i + 1) for i in range(len(input_files))]},
    "out": ["gene_output", "transcript_output"]}

    yaml["cwl"]["steps"]["prepDE_folder"] = {
    "run": {
        "class": "ExpressionTool",
        "requirements": {
            "InlineJavascriptRequirement": {}},
        "inputs": dict([("file{0}".format(i + 1), "File") for i in range(len(input_files))]),
        "outputs": {
            "out": "Directory"},
        "expression": "|${{return{{\"out\":{{\"class\":\"Directory\",\"basename\":\"prepDE\",\"listing\":[{}]}};}}}}".format(",".join(["inputs.file{0}".format(i + 1) for i in range(len(input_files))]))},
    "in": {
        "file1": "prepDE/gene_output",
        "file2": "prepDE/transcript_output"},
    "out": ["out"]
    }

    yaml["yml"]["prepDE_script"] = {
        "class": "File",
        "path": "../scripts/prepDE.py"}
    return yaml

def stringtie():
    pass

def deseq2():
    pass



input_file = {"test1":{"type": "paired_end", "path":{1: "test.fastq", 2: "test.fastq"}},
              "test2":{"type": "paired_end", "path":{1: "test.fastq", 2: "test.fastq"}}}
ale = star_readmap(input_file, cwl_workflow)
ale = samtools(input_file, ale)
ale = prepDE(input_file, ale)
with open("test.cwl", "w+") as outfile:
    yaml.dump(ale['cwl'], outfile, default_flow_style=False)
with open("test.yml", "w+") as outfile:
    yaml.dump(ale['yml'], outfile, default_flow_style=False)
