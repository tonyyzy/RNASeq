import subprocess
import filecmp

def test_fgsea_docker():
	subprocess.run(["cwl-runner",
			"--outdir=./test_fgsea_out",
			"./cwl-tools/docker/fgsea.cwl",
			"./tests/fgsea.yml"])
	assert filecmp.cmp("./tests/gsea_res.star_prepde_deseq2.csv", "./test_fgsea_out/gsea_res.csv")

if __name__ == "__main__":
	test_fgsea_docker()
