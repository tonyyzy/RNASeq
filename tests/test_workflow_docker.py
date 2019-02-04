import subprocess
import filecmp

def test_workflow_docker():

    subprocess.run(["cwl-runner",
                    "--outdir=./test_workflow_docker",
                    "./workflows/docker/star_samtools_stringtie_prepDE_DESeq2.cwl",
                    "./tests/star_samtools_stringtie_prepDE_DESeq2.yml"])
        
    assert filecmp.cmp("./test_workflow_docker/DGE_results.csv", "./tests/DGE_results.csv")


if __name__ == "__main__":
    test_workflow_docker()
