import os
import pandas as pd
import pickle
import numpy as np
import sys
from sklearn.model_selection import train_test_split

def SplitData(onehot, eff_labels):
    random_state=682
    test_size = 0.1
    x_train_onehot, x_test_onehot, eff_y_train, eff_y_test = train_test_split(onehot, eff_labels, test_size=test_size, random_state=random_state)
    data = {'eff-input':
            {'train':{'onehot':x_train_onehot},
             'test':{'onehot':x_test_onehot}
             },
            'eff-label':
            {'train':eff_y_train,
             'test':eff_y_test
             }
             }
           
    return data

def encode(sgrna):
    length = 24
    code_dict = {'A': [[1], [0], [0], [0]], 'T': [[0], [1], [0], [0]], 'G': [[0], [0], [1], [0]], 'C': [[0], [0], [0], [1]]}     
    onehot = []
    for i in range(len(sgrna)):    
        mers = []
        for j in range(length):
            mers.append(code_dict[sgrna[i][j]])
        onehot.append(mers)
    np_onehot = np.array(onehot).reshape(-1,length,4,1)    
    return np_onehot

def saveonehot(inputs):
    length = 24
    sgrnas = []
    with open(inputs, 'r') as fb:
        next(fb)
        lines = fb.readlines()
    
    for line in lines:
        a = line.strip().split(',')
        if len(a[0])==length:
            sgrnas.append(a[0])

    #print('sgrnas: ', len(sgrnas))
    np_onehot = encode(sgrnas)
    #print('np_onehot: ', np_onehot.shape)

    return np_onehot

def saveeff(inputs,index=1):
    length = 24
    f=open(inputs,"r")
    line = f.readline()
    #print('line: ', line)
    i = 0
    eff = []
    while line:
        i += 1
        line = f.readline().replace('\n','')
        a = line.split(',')
        if len(a[0])==length:
            if len(a)<=index or a[index]=='':
                eff.append(-1)
            else:
                eff.append(float(a[index]))
    np_eff = np.array(eff,dtype = float)
    #print('np_eff: ', np_eff.shape)
    f.close()

    return np_eff
    
def sequencing(sgrnas):
    length = 24
    seq_dict = {'T': 1, 'A': 2, 'C': 3, 'G': 4,'START': 0}     
    seqs = []
    for i in range(len(sgrnas)):    
        mers = []
        #mers.append(seq_dict['START'])
        for j in range(length):
            mers.append(seq_dict[sgrnas[i][j]])
        seqs.append(mers)
    #print('mers: ',len(mers))
    np_seqs = np.array(seqs).reshape(-1,length)
    return np_seqs

def saveseq(inputs):
    length = 24
    sgrnas = []
    with open(inputs, 'r') as fb:
        next(fb)
        lines = fb.readlines()
    
    for line in lines:
        a = line.strip().split(',')
        if len(a[0])==length:
            sgrnas.append(a[0])
    #print(len(sgrnas))

    np_seqs = sequencing(sgrnas)

    #f.close()
    #fo = open('./sequence.pkl','wb')
    #pickle.dump(np_seqs,fo)
    #fo.close()
    return np_seqs

if __name__ == '__main__':
    assert len(sys.argv)==3

    inputs = sys.argv[1]
    outputs = sys.argv[2]
    seq = saveseq(inputs)
    eff_labels = saveeff(inputs)
    onehot = saveonehot(inputs)
    print('onehot.shape: ',onehot.shape, ', eff.shape: ', eff_labels.shape)
    #savebiofeature()
    #data = {'onehot':onehot, 'seq':seq, 'eff':eff}
    data = SplitData(onehot, eff_labels)
    with open(outputs, 'wb') as fb:
        pickle.dump(data, fb)