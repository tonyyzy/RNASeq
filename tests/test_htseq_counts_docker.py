import subprocess
import filecmp

def test_htseq_counts_docker():
	subprocess.run(["cwl-runner",
			"--outdir=./test_htseq_counts_out",
			"./cwl-tools/docker/htseq.cwl",
			"./tests/htseq_counts.yml"])
	assert filecmp.cmp("./tests/test1_htseq_count.csv", "./test_htseq_counts_out/test1_htseq_count.csv")

if __name__ == "__main__":
	test_htseq_counts_docker()
