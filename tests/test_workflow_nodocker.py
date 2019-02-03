import subprocess
import filecmp

def test_workflow_nodocker():

    subprocess.run(["cwl-runner",
                    "--outdir=./test_workflow_nodocker",
                    "./workflows/nodocker/star_samtools_stringtie_prepDE_DESeq2.cwl",
                    "./tests/star_samtools_stringtie_prepDE_DESeq2.yml"])
        
    assert filecmp.cmp("./test_workflow_nodocker/DGE_results.csv", "./tests/DGE_results.csv")


if __name__ == "__main__":
    test_workflow_nodocker()
