import subprocess
import filecmp

def test_STAR_readmap_gz_docker():
    """test_STAR_readmap_gz_docker"""
    subprocess.run(["cwl-runner",
                    "./cwl-tools/docker/STAR_readmap.cwl",
                    "./tests/STAR_readmap_gz.yml"])
    subprocess.run(["tail -n +5 ./test1gz_STARAligner/test1gzAligned.out.sam > ./test1gz_STARAligner/test1gz.tail.sam"], shell=True)
    assert filecmp.cmp("./test1gz_STARAligner/test1gz.tail.sam", "./tests/test1.tail.sam")


if __name__ == "__main__":
    test_STAR_readmap_gz_docker()
