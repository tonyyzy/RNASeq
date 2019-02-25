import subprocess
import filecmp

def test_salmon_quant_paired():
	subprocess.run(["cwl-runner",
					"--outdir=./test_salmon_quant",
					"./cwl-tools/docker/salmon_quant.cwl",
					"./tests/salmon_quant.yml"])
	subprocess.run(["""tail -n +2 ./test_salmon_quant/test2/quant.sf |
					 awk 'BEGIN{OFS=FS="\t"}{$3=sprintf("%3.0f",$3);
						 $4=sprintf("%3.0f",$4)}1' > ./test_salmon_quant/test2.sf"""], shell=True)
	assert filecmp.cmp("./tests/salmon_quant/test2.sf",
						"./test_salmon_quant/test2.sf")

def test_salmon_quant_single():
	subprocess.run(["cwl-runner",
					"--outdir=./test_salmon_quant",
					"./cwl-tools/docker/salmon_quant.cwl",
					"./tests/salmon_quant.single.yml"])
	subprocess.run(["""tail -n +2 ./test_salmon_quant/test3/quant.sf |
					 awk 'BEGIN{OFS=FS="\t"}{$3=sprintf("%3.0f",$3);
						 $4=sprintf("%3.0f",$4)}1' > ./test_salmon_quant/test3.sf"""], shell=True)
	assert filecmp.cmp("./tests/salmon_quant/test3.sf",
						"./test_salmon_quant/test3.sf")


if __name__ == "__main__":
	test_salmon_quant_paired()
	test_salmon_quant_single()