import subprocess
import filecmp

def test_STAR_readmap_docker():
    """test_STAR_readmap_nodocker"""
    subprocess.run(["cwl-runner",
                    "--outdir=./test1d",
                    "./cwl-tools/docker/STAR_readmap.cwl",
                    "./tests/STAR_readmap.yml"])
    subprocess.run(["tail -n +5 ./test1d/test1Aligned.out.sam > ./test1d/test1.tail.sam"], shell=True)
    assert filecmp.cmp("./test1d/test1.tail.sam", "./tests/test1.tail.sam")


if __name__ == "__main__":
    test_STAR_readmap_docker()
