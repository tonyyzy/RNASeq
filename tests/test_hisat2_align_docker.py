import subprocess
import filecmp

def test_hisat2_align_docker():
    """test_hisat2_align_docker"""
    subprocess.run(["cwl-runner",
                    "--outdir=./test_hisat2_align",
                    "./cwl-tools/docker/hisat2_align.cwl",
                    "./tests/hisat2_align.yml"])
    subprocess.run(["tail -n +5 ./test_hisat2_align/test1/test1.sam > ./test_hisat2_align/test1.tail.sam"], shell=True)
    assert filecmp.cmp("./tests/test1.hisat2.tail.sam","./test_hisat2_align/test1.tail.sam")

if __name__ == "__main__":
    test_hisat2_align_docker()
