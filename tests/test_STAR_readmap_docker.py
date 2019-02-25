import subprocess
import filecmp

def test_STAR_readmap_docker():
    """test_STAR_readmap_nodocker"""
    subprocess.run(["cwl-runner",
                    "--outdir=./test_star_readmap_docker",
                    "./cwl-tools/docker/STAR_readmap.cwl",
                    "./tests/STAR_readmap_docker.yml"])
    subprocess.run(["tail -n +5 ./test_star_readmap_docker/test3/test3Aligned.out.sam "
                    "> ./test_star_readmap_docker/test3/test3.tail.sam"], shell=True)
    assert filecmp.cmp("./test_star_readmap_docker/test3/test3.tail.sam",
                        "./tests/test3.tail.sam")


if __name__ == "__main__":
    test_STAR_readmap_docker()
