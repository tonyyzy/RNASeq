import subprocess
import filecmp

def test_DESeq2_docker():

    subprocess.run(["cwl-runner",
                    "--outdir=./test_DESeq2_out",
                    "./cwl-tools/docker/DESeq2.cwl",
                    "./tests/DESeq2.yml"])

    assert filecmp.cmp("./tests/DGE_res.star_prepde.csv", "./test_DESeq2_out/untreated-treated_DGE_res.csv")

if __name__ == "__main__":
    test_DESeq2_docker()
