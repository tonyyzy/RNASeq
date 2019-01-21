# Script to generate fastq files
import random


def fake(exon, num):
	global sam
	for i in range(num):
		file1.write(f"@sample{sam}/1\n")
		file2.write(f"@sample{sam}/2\n")
		left_start = random.randint(0, len(exon) - 302)
		file1.write(exon[left_start - 1:left_start + 100] + "\n+\n")
		file1.write("".join(chr(i) for i in random.choices(range(94, 110), k=101)) + "\n")
		file2.write(exon[left_start + 200:left_start + 301] + "\n+\n")
		file2.write("".join(chr(i) for i in random.choices(range(94, 110), k=101)) + "\n")
		sam += 1


if __name__ == "__main__":
	seq = ""
	sam = 1
	with open("./test.fa") as fasta:
		fasta.readline()
		for i in fasta.readlines():
			seq += i.strip()

	exon1 = seq[499:999]
	exon2 = seq[2999:3299] + seq[3499:3999] + seq[4399:4999] + seq[5049:5299]

	file1 = open("./test1.1.fastq", "w")
	file2 = open("./test1.2.fastq", "w")
	fake(exon1, 20)
	fake(exon2, 30)

	file1 = open("./test2.1.fastq", "w")
	file2 = open("./test2.2.fastq", "w")
	fake(exon1, 40)
	fake(exon2, 60)

