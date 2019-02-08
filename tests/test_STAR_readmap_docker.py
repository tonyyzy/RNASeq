import subprocess
import filecmp

def test_STAR_readmap_docker():
    """test_STAR_readmap_nodocker"""
    subprocess.run(["cwl-runner",
                    "./cwl-tools/docker/STAR_readmap.cwl",
                    "./tests/STAR_readmap_docker.yml"])
    subprocess.run(["tail -n +5 ./test3/test3Aligned.out.sam > ./test3/test3.tail.sam"], shell=True)
    assert filecmp.cmp("./test3/test3.tail.sam", "./tests/test3.tail.sam")


if __name__ == "__main__":
    test_STAR_readmap_docker()
