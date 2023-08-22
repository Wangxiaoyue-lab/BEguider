import itertools

def generate_combination():
    '''
    产生sgRNA的pos4-pos8位，所有碱基序列的组合
    num2seq = {0:'AAAAA', ...}
    seq2num = {'AAAAA':0, ...}
    '''
    size = 5
    list1 = ['A', 'T', 'C', 'G']
    seqs = [''.join(i) for i in itertools.product(list1, repeat=size)]
    nums = list(range(4**size))
    num2seq = {}
    seq2num = {}
    for m,n in zip(seqs,nums):
        num2seq[n] = m
        seq2num[m] = n
    return num2seq, seq2num

num2seq, seq2num = generate_combination()
if 'GCTCT' in list(seq2num.keys()):
	print(seq2num['GCTCT'])
