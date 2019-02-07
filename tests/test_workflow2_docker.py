import subprocess
import filecmp

def test_workflow2_docker():

    subprocess.run(["cwl-runner",
                    "--outdir=./test_workflow2_docker",
                    "./workflows/docker/hisat2_htseq_dexseq_four_inputs.cwl",
                    "./tests/hisat2_htseq_dexseq_four_inputs.yml"])
        
    assert filecmp.cmp("./test_workflow2_docker/DEXSeq/DEE_results.csv", "./tests/DEE_results.csv")


if __name__ == "__main__":
    test_workflow2_docker()
