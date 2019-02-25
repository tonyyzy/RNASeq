import subprocess
import filecmp

def test_salmon_index_docker():
	subprocess.run(["cwl-runner",
					"--outdir=./test_salmon_index",
					"./cwl-tools/docker/salmon_index.cwl",
					"./tests/salmon_index.yml"])
	assert filecmp.cmp("./tests/Salmonindex/sa.bin",
						"./test_salmon_index/Salmonindex/sa.bin")

if __name__ == "__main__":
	test_salmon_index_docker()
