import subprocess
import filecmp

def test_hisat2_index_docker():
    """test_hisat2_index_docker"""
    subprocess.run(["cwl-runner",
                    "--outdir=./test_hisat2_index",
                    "./workflows/docker/hisat2_index.cwl",
                    "./tests/hisat2_index.yml"])
    assert filecmp.cmp("./tests/HISAT2Index/test.1.ht2", 
                        "./test_hisat2_index/HISAT2Index/test.1.ht2")

if __name__ == "__main__":
    test_hisat2_index_docker()
