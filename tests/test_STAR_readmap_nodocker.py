import subprocess
import filecmp

def test_STAR_readmap_nodocker():
    """test_STAR_readmap_nodocker"""
    subprocess.run(["cwl-runner",
                    "./cwl-tools/nodocker/STAR_readmap.cwl",
                    "./tests/STAR_readmap.yml"])
    subprocess.run(["tail -n +5 ./test1/test1Aligned.out.sam > ./test1/test1.tail.sam"], shell=True)
    subprocess.run(["cat -A ./test1/test1.tail.sam"], shell=True)
    assert filecmp.cmp("./test1/test1.tail.sam", "./tests/test1.tail.sam")


if __name__ == "__main__":
    test_STAR_readmap_nodocker()
