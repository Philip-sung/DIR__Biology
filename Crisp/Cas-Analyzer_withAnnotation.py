#모듈 가져오기
from needle import *

# 변수 설정
ref_seq = 'ATGCATTTTTTTGAGACAAGGTCTTGCTCTATTGTCCAGGCTGGAGTGCAGTGGCACAATCACAGTTCACTCCAGCCTCAACATCCTGCACTAAAGTGATTTTCCCACCTCACCTCTCAAGTAGCTGGGACTACAGGTACATGCTACCATGCCTGGCTAATTTTTTTTTTTTTGCAGGCATGGGGTCTCACTGTATTGCCCAGGTTGGTGTGGAAGTTTAATGACTAAGAGGTGTTTGTTATAAAGTTTAATGTATGAAACTTTCTATTAAATTCCTGATTTTATTTCTGTAGGACTGAACGTCTTGCTCGAGATGTGATGAAGGAGATGGGAGGCCATCACATTGTAGCCCTCTGTGTGCTCAAGGGGGGCTATAAATTCTTTGCTGACCTGCTGGATTACATCAAAGCACTGAATAGAAATAGTGATAGATCCATTCCTATGACTGTAGATTTTATCAGACTGAAGAGCTATTGTGTGAGTATATTTAATATATGATTCTTTTTAGTGGCAACAGTAGGTTTTCTTATATTTTCTTTGAATCTCTGCAAACCATACTTGCTTTCATTTCACTTGGTTACAGTGAGATTTTTCTAACATATTCACTAGTACTTTACATCAAAGCCAATACTGTTTTTTTAAAACTAGTCACCTTGGAGGATATATACTTATTTTACAGGTGTGTGTGGTTTTTTAAATAAACTCCTTTTAGGAATTGCTGTTG'.upper()
target_seq = 'GCCCCCCTTGAGCACACAGAGGG'.upper()
R = 70
min_n = 1
r = 5

#fastq파일 오픈하여 10만줄까지 읽기
f = open('85.fastqjoin').readlines()[:100000]

#
def rev_comp(s): return s.translate(s.maketrans('ATGC','TACG'))[::-1]

if ref_seq.find(target_seq) != -1:
	cv_p = ref_seq.find(target_seq) + 16
elif ref_seq.find(rev_comp(target_seq)) != -1:
	cv_p = ref_seq.find(rev_comp(target_seq)) + 6
else:
	print('ERROR!! can not find target sequence in ref')

st = cv_p - R
ed = cv_p + R

if st < 0:
	st = 0

if ed > len(ref_seq):
	ed = len(ref_seq)

range_seq = ref_seq[st:ed]
wt_marker = ref_seq[cv_p - r: cv_p + r]

indi_f = ref_seq[st: st + 15]
indi_r = ref_seq[ed - 15: ed]

indi_f_list = []
indi_r_list = []

for i in range(15):
	for nt in 'ATGC':
		indi_f_list.append(indi_f[:i] + nt + indi_f[i + 1:])
		indi_r_list.append(indi_r[:i] + nt + indi_r[i + 1:])

seq_dict = {}
cnt_dict = {'all': 0, 'indi': 0, 'min': 0, 'wt': 0, 'del': 0, 'ins': 0, 'other': 0}

fw_ins = open('ins.txt', 'w')
fw_del = open('del.txt', 'w')

for line_n, line in enumerate(f):
	
	if line_n % 4 != 1:
		continue
	
	cnt_dict['all'] += 1
	
	line = line.strip()

	seq_st = False
	seq_ed = False

	for i in indi_f_list:
		if line.find(i) != -1:
			seq_st = line.find(i)
			break
	
	for i in indi_r_list:
		if line.find(i) != -1:
			seq_ed = line.find(i) + 15
			break
	
	if seq_st == False or seq_ed == False:
		continue
	
	cnt_dict['indi'] += 1

	seq = line[seq_st: seq_ed]

	if seq not in seq_dict:
		seq_dict[seq] = 1
	else:
		seq_dict[seq] += 1

for seq, cnt in seq_dict.items():

	if cnt <= min_n:
		continue

	cnt_dict['min'] += cnt

	needle_res = needle(range_seq, seq, 10, 10, 0.5, 0.5)

	if seq.find(wt_marker) != -1:
		cnt_dict['wt'] += cnt
	else:
		if len(seq) == len(range_seq):
			cnt_dict['wt'] += cnt
		elif len(seq) < len(range_seq):
			cnt_dict['del'] += cnt
			fw_del.write(needle_res[0] + '\t' + str(cnt) + '\n' + '\n'.join(needle_res[1:3]) + '\n')
		elif len(seq) > len(range_seq):
			cnt_dict['ins'] += cnt
			fw_ins.write(needle_res[0] + '\t' + str(cnt) + '\n' + '\n'.join(needle_res[1:3]) + '\n')
		else:
			cnt_dict['other'] += cnt

fw_del.close()
fw_ins.close()

print(cnt_dict)

