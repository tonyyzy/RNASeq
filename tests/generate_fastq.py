# Script to generate fastq files
import random


def fake_pairend(exon, num):
	global sam
	for i in range(num):
		file1.write(f"@sample{sam}/1\n")
		file2.write(f"@sample{sam}/2\n")
		left_start = random.randint(0, len(exon) - 222)
		file1.write(exon[left_start:left_start + 101] + f"\n+{left_start}\n")
		file1.write("".join(chr(i) for i in random.choices(range(94, 110), k=101)) + "\n")
		file2.write("".join([complement[base] for base in exon[left_start + 221:left_start + 120:-1]]) + f"\n+{left_start}\n")
		file2.write("".join(chr(i) for i in random.choices(range(94, 110), k=101)) + "\n")
		sam += 1


def fake_single(exon, num):
	global sam
	for i in range(num):
		file1.write(f"@sample{sam}\n")
		left_start = random.randint(0, len(exon) - 101)
		file1.write(exon[left_start:left_start + 101] + f"\n+{left_start}\n")
		file1.write("".join(chr(i) for i in random.choices(range(94, 110), k=101)) + "\n")
		sam += 1


if __name__ == "__main__":
	complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
	seq = ""
	sam = 1
	with open("./test.fa") as fasta:
		fasta.readline()
		for i in fasta.readlines():
			seq += i.strip()

	exon1 = seq[499:999]
	# exon2 = seq[2999:3299] + seq[3499:3999] + seq[4399:4999] + seq[5049:5299]
	exon2_1 = seq[2999:3299] + seq[4399:4999] + seq[5049:5299]
	exon2_2 = seq[2999:3299] + seq[3499:3999] + seq[5049:5299]

	file1 = open("./test1.1.fastq", "w")
	file2 = open("./test1.2.fastq", "w")
	fake_pairend(exon1, 10)
	fake_pairend(exon2_1, 90)
	fake_pairend(exon2_2, 30)

	file1 = open("./test2.1.fastq", "w")
	file2 = open("./test2.2.fastq", "w")
	fake_pairend(exon1, 70)
	fake_pairend(exon2_2, 91)
	fake_pairend(exon2_1, 30)

	file1 = open("./test3.fastq", "w")
	fake_single(exon1, 15)
	fake_single(exon2_1, 87)
	fake_single(exon2_2, 30)

	file1 = open("./test4.fastq", "w")
	fake_single(exon1, 77)
	fake_single(exon2_2, 93)
	fake_single(exon2_1, 30)

