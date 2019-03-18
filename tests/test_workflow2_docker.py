import subprocess
import filecmp

def test_workflow2_docker():

    subprocess.run(["cwl-runner",
                    "--outdir=./test_workflow2_docker",
                    "./workflows/docker/hisat2_htseq_dexseq.cwl",
                    "./tests/hisat2_htseq_dexseq.yml"])
        
    assert filecmp.cmp("./test_workflow2_docker/DEXSeq/untreated-treated_DEE_results.csv", "./tests/DEE_results.hisat2_dexseq.csv")


if __name__ == "__main__":
    test_workflow2_docker()
