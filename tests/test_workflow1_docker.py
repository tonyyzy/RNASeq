import subprocess
import filecmp

def test_workflow1_docker():
    """test_workflow1_docker"""
    subprocess.run(["cwl-runner",
                    "--outdir=./test_workflow1_docker",
                    "./workflows/docker/star_samtools_stringtie_prepDE_DESeq2.cwl",
                    "./tests/star_samtools_stringtie_prepDE_DESeq2.yml"])
        
    assert filecmp.cmp("./test_workflow1_docker/DESeq2/groupuntreated-grouptreated_DGE_results.csv", "./tests/DGE_results.csv")


if __name__ == "__main__":
    test_workflow1_docker()
