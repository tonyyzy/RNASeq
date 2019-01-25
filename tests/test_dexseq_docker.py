import subprocess
import filecmp

def test_dexseq_docker():
	subprocess.run(["cwl-runner",
			"--outdir=./test_dexseq_out",
			"./cwl-tools/docker/dexseq.cwl",
			"./tests/dexseq.yml"])
	assert filecmp.cmp("./tests/test_DEE_results.csv", "./test_dexseq_out/DEE_results.csv")

if __name__ == "__main__":
	test_dexseq_docker()
