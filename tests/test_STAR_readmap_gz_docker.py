import subprocess
import filecmp

def test_STAR_readmap_gz_docker():
    """test_STAR_readmap_gz_docker"""
    subprocess.run(["cwl-runner",
                    "--outdir=./test_star_readmap_gz",
                    "./cwl-tools/docker/STAR_readmap.cwl",
                    "./tests/STAR_readmap_gz.yml"])
    subprocess.run(["tail -n +5 ./test_star_readmap_gz/test1gz/test1gzAligned.out.sam "
                    "> ./test_star_readmap_gz/test1gz/test1gz.tail.sam"], shell=True)
    assert filecmp.cmp("./test_star_readmap_gz/test1gz/test1gz.tail.sam",
                        "./tests/test1.star.tail.sam")


if __name__ == "__main__":
    test_STAR_readmap_gz_docker()
