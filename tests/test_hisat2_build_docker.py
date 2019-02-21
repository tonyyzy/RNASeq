import subprocess
import filecmp

def test_hisat2_build_docker():
    """test_hisat2_build_docker"""
    subprocess.run(["cwl-runner",
                    "--outdir=./test_hisat2_build",
                    "./cwl-tools/docker/hisat2_build.cwl",
                    "./tests/hisat2_build.yml"])
    assert filecmp.cmp("./tests/ht2_index/test_index.1.ht2", 
                        "./test_hisat2_build/test_index.1.ht2")

if __name__ == "__main__":
    test_hisat2_build_docker()
