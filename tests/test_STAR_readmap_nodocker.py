import subprocess
import filecmp

def test_STAR_readmap_nodocker():
    """test_STAR_readmap_nodocker"""
    subprocess.run(["cwl-runner",
                    "./cwl-tools/nodocker/STAR_readmap.cwl",
                    "./tests/STAR_readmap_nodocker.yml"])
    subprocess.run(["tail -n +5 ./test1_STARAligner/test1Aligned.out.sam > ./test1_STARAligner/test1.tail.sam"], shell=True)
    assert filecmp.cmp("./test1_STARAligner/test1.tail.sam", "./tests/test1.tail.sam")


if __name__ == "__main__":
    test_STAR_readmap_nodocker()
