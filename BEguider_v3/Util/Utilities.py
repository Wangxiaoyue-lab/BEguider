#!/usr/bin/env python3
# -*- coding utf-8 -*-

import itertools
import re
import os
import sys
import platform
import numpy as np
import pandas as pd
from subprocess import Popen
from tensorflow.keras.models import load_model
from fetch_single_rsid import get_snp
from fetch_genomic_seq import *

def GetParams(dataset):
    if dataset == 'ABE7.10-NGG':
        effpath = './saved_model/efficiency/ABE710NGG_BestCNN.h5'
        proppath = './saved_model/proportion/ABE710NGG_BestCNN.h5'
    elif dataset == 'BE4-NGG':
        effpath = './saved_model/efficiency/BE4NGG_BestCNN.h5'
        proppath = './saved_model/proportion/BE4NGG_BestCNN.h5'
    elif dataset == 'ABEmax-SpRY':
        effpath = './saved_model/efficiency/ABEmax_BestENS.h5'
        proppath = './saved_model/proportion/ABEmax_BestENS.h5'
    elif dataset == 'ABE8e-P(AP)3-SpRY':
        effpath = './saved_model/efficiency/ABE8e-7_BestENS.h5'
        proppath = './saved_model/proportion/ABE8e-7_BestENS.h5'
    elif dataset == 'ABE8e-NL-SpRY':
        effpath = './saved_model/efficiency/ABE8eNL_BestENS.h5'
        proppath = './saved_model/proportion/ABE8eNL_BestENS.h5'
    elif dataset == 'BE4max-SpRY':
        effpath = './saved_model/efficiency/CBE4max_BestENS.h5'
        proppath = './saved_model/proportion/CBE4max_BestENS.h5'
    elif dataset == 'FNLS-YE1-BE4max-SpRY':
        effpath = './saved_model/efficiency/CBE-FNLS_BestENS.h5'
        proppath = './saved_model/proportion/CBE-FNLS_BestENS.h5'
    elif dataset == 'YE1-BE4max-SpRY':
        effpath = './saved_model/efficiency/CBE-YE1_BestENS.h5'
        proppath = './saved_model/proportion/CBE-YE1_BestENS.h5'
    else:
        print('unsupport the base editor')
        sys.exit()
    return effpath, proppath

class Result(object):
    Best = -1
    

def off_txt_header(offpath):
    #写入hg38所在的文件夹路径，不能写到文件名
    if '\\' in offpath:
        hg38 = '.\\Util\\'
    else:
        hg38 = './Util/'
    
    pattern = 'NNNNNNNNNNNNNNNNNNNNNGG'

    with open(offpath, 'w') as files:
        files.write(hg38+'\n')
        files.write(pattern+'\n')

def process_dir(BaseEditor, seq_path, outpath):
    if outpath is not None:
        cwd = outpath
    else:
        cwd = os.getcwd()

    if (platform.system()=='Windows'):
        name = seq_path.split('\\')[-1].split('.')[0]
        effout = cwd + '\\outputs\\' + name +  '_Designed_sgRNA_for'+'_'+BaseEditor+'.csv'
        propout = cwd + '\\outputs\\' + name + '_Predicted_Editing_Proportion_for'+'_'+BaseEditor+'.csv'
        outpath = cwd + '\\outputs\\'
        tempdir = cwd+ '\\outputs\\temp'
        offpath = tempdir + '\\seqs_to_offtarget.txt'
    else:
        name = seq_path.split('/')[-1].split('.')[0]
        effout = cwd + '/outputs/' + name +  '_designed_sgRNA_for'+'_'+BaseEditor+'.csv'
        propout = cwd + '/outputs/' + name + '_Predicted_Editing_Proportion_for'+'_'+BaseEditor+'.csv'
        outpath = cwd + '/outputs/'
        tempdir = cwd+ '/outputs/temp'
        offpath = tempdir + '/seqs_to_offtarget.txt'

    outfile = (effout, propout)
    if "outputs" not in os.listdir(cwd):
        os.makedirs(outpath)
    
    if "temp" not in os.listdir(outpath):
        os.makedirs(tempdir)

    off_txt_header(offpath)

    return outfile, tempdir, offpath, name


def get_offtarget(fname, tempdir):
    '''
    cmd -- C:CPU
    '''
    dirs =  tempdir[:-4]
    inpath = tempdir + '/'+fname+'_seqs_to_offtarget.txt'
    outpath = dirs +'/'+fname+'_offtarget_result.txt'

    if (platform.system()=='Windows'):
        inpath = inpath.replece('/', '\\')
        outpath = outpath.replece('/', '\\')
        cmd = '.\\Util\\win-cas-offinder.exe ' + inpath + ' C' + ' ' + outpath
        #print('Windows:', cmd)
    elif (platform.system()=='Linux'):
        cmd = './Util/linux-cas-offinder ' + inpath + ' C' + ' ' + outpath
        cmd2 = 'chmod u+x ../Util/linux-cas-offinder'
        chmo = Popen(args=cmd2, shell=True)
        chmo.wait()
        #print('Linux:', cmd)
    elif (platform.system()=='Darwin'):
        cmd = './Util/mac-cas-offinder ' + inpath + ' C' + ' ' + outpath
        cmd2 = 'chmod u+x ../Util/mac-cas-offinder'
        chmo = Popen(args=cmd2, shell=True)
        chmo.wait()
        #print('MacOS:', cmd)
    else:
        print('The system does not support searching for off-target sites.')
        print('Please use Windows/Linux/MacOS platform.')
        sys.exit()

    off = Popen(args=cmd, shell=True)
    off.wait() #等待进程完成
    #print('off.returncode:' ,off.returncode)
    if off.returncode == 0:
        print('Off-target finished.')
    return outpath

def general(path):
    with open(path, 'r') as files:
        next(files)
        lines = files.readlines()

    return lines

def get_genomic_seq(pam, path, tempdir):
    lines = general(path)
    outpathABE = tempdir +'/genomic_seqs_for_ABE.txt'
    outpathCBE = tempdir +'/genomic_seqs_for_CBE.txt'

    if (platform.system()=='Windows'):
        outpath = outpath.replace('/', '\\')

    seqs = {}
    for line in lines:
        info = line.strip().split(',')
        chrom = info[0]
        site = int(info[1])
        try:
            seq = fetch_dna_coordinates('hg38', chrom, site-20, site+20, './')
            #print('seq:\n', seq)
        except:
            print('The infomation of SNP '+chrom+':'+info[1]+' was not found.')
            continue
        #etype = 'r'
        etype = info[2]
        seq = snp_seq(pam, etype, seq)
        #print('seq:',seq)
        if seq[0] == 0:
            continue
        else:
            key = chrom+':'+str(site)
            seqs[key] = seq
    #print('seqs[key]:', key,seq)
    tag = split_seq_by_BE(seqs, outpathABE, outpathCBE)

    return (outpathABE, outpathCBE), tag

def split_seq_by_BE(seqs, outpathABE, outpathCBE):
    fa = open(outpathABE, 'w')
    fc = open(outpathCBE, 'w')
    col = 'Info,Seqs\n'
    seqABE = []
    seqCBE = []
    for key, value in seqs.items():
        x = [','.join(i)+'\n' for i in itertools.product([key],value[0])]
        if value[1] == 'ABE':
            seqABE += x
        elif value[1] == 'CBE':
            seqCBE += x
    
    if len(seqABE) > 0 :
        fa.write(col)
        fa.writelines(seqABE)
        tagABE = 1
    else:
        tagABE = 0
    if len(seqCBE) > 0:
        fc.write(col)
        fc.writelines(seqCBE)
        tagCBE = 1
    else:
        tagCBE = 0

    fa.close()
    fc.close()
    #print('tagABE, tagCBE:',tagABE, tagCBE)
    return tagABE, tagCBE

def identify_strand(pam, seq):
    site = seq[20]
    rev = {'G':'C','T':'A'}
    chooseBE = {'C':'CBE','G':'CBE','A':'ABE','T':'ABE'}
    if site in ['C', 'A']:
        #print('seq[13:]:',seq[13:])
        be = chooseBE[site]
        if pam.upper() == 'NGG':
            new_seq = get_pam_seq(seq[13:], substrate = site)
        elif pam.upper() == 'FREE':
            new_seq = get_pamfree_seq(seq[13:], substrate = site)
    
    elif site in ['T', 'G']:
        #print('seq[:28]:',seq[:28])
        be = chooseBE[site]
        rev_seq = reverse_seq(seq[:28])
        
        #print('new_seq:', new_seq)
        if pam.upper() == 'NGG':
            new_seq = get_pam_seq(rev_seq, substrate = rev[site])
        elif pam.upper() == 'FREE':
            new_seq = get_pamfree_seq(rev_seq, substrate = rev[site])

        if new_seq !=0:
            temp = [ reverse_seq(i) for i in new_seq]
            new_seq = temp
            #print('temp:', temp)
    else:
        new_seq = 0
        be = 0
        
    return [new_seq,be]

def snp_seq(pam, etype, seq):
    
    allele = {'A':'G', 'G':'A', 'C':'T', 'T':'C'}
    if etype == 'r':
        new_seq = identify_strand(pam, seq)
        
    elif etype == 'a':
        new_seq = seq[:20] + allele[seq[20]] + seq[21:]
        new_seq = identify_strand(pam, new_seq)
    return new_seq

def get_snp_seq(pam, path, tempdir):
    lines = general(path)
    outpathABE = tempdir +'/snp_seqs_for_ABE.txt'
    outpathCBE = tempdir +'/snp_seqs_for_CBE.txt'

    if (platform.system()=='Windows'):
        outpath = outpath.replace('/', '\\')

    seqs = {}
    for line in lines:
        line = line.strip().split(',')
        snp = line[0]
        etype = line[1]
        try:
            seqinfo = get_snp(snp, tempdir)
        except:
            print('The infomation of SNP '+snp+' was not found.')
            continue
        chrom = seqinfo[0]
        site = seqinfo[1]
        seq = snp_seq(pam, etype, seqinfo[2])
        #print('seqinfo:',seqinfo)
        if seq[0] == 0:
            continue
        else:
            seqs[snp+'_chr'+chrom+':'+str(site)] = seq
        #print('seqinfo:',seqinfo)
    #print(seqs)
    tag = split_seq_by_BE(seqs, outpathABE, outpathCBE)

    return (outpathABE, outpathCBE), tag

def encode(sgrna):
    length  = 24
    code_dict = {'A': [[1], [0], [0], [0]], 'T': [[0], [1], [0], [0]], 'G': [[0], [0], [1], [0]], 'C': [[0], [0], [0], [1]]}
    onehot = []
    for i in range(len(sgrna)):
        mers = []
        for j in range(length):
            mers.append(code_dict[sgrna[i][j]])
        onehot.append(mers)
    np_onehot = np.array(onehot).reshape(-1,length,4,1)
    return np_onehot

def reverse_seq(sequence):
    temp = []
    rev = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    for i in reversed(sequence):
        temp.append(rev[i])

    rev_seq = ''.join(i for i in temp)
    return rev_seq

def activate_win(seq, substrate):
    
    win_start = 4
    win_end = 8

    #print('substrate:', substrate)
    if substrate in seq[win_start-1:win_end]:
        tag = 1
    else:
        tag = 0
    return tag

def get_pam_seq(sequence, substrate):
    '''
    0 is a label, which means there is no NGG PAM in the seq
    '''
    #print('sequence:',sequence)
    length = 24
    indexs = [m.start() for m in re.finditer('(?=GG)', sequence) if m.start() >= 21]
    seqs = []
    #print('index: ',indexs)
    if len(indexs) == 0:
        return 0
    elif len(indexs) > 0:
        for idx in indexs:
            seq = sequence[idx-21:idx+3]
            tag = activate_win(seq, substrate)
            
            if len(seq) != length or tag != 1:
                continue
            elif len(seq) == length and tag == 1:
                seqs.append(seq)

        if len(seqs) == 0:
            return 0
        else:
            return seqs

def get_pamfree_seq(sequence, substrate):
    length = 24
    indexs = [m.start() for m in re.finditer(substrate, sequence) if m.start() >= 3 ]
    #print('sequence:',sequence)
    #print('substrate:',substrate)
    
    seqs = []
    if len(indexs) == 0:
        return 0
    elif len(indexs) > 0:
        for idx in indexs:
            for pos in range(3,8):
                seq = sequence[idx-pos:idx+24-pos]
                if len(seq) != length:
                    continue
                elif len(seq) == length:
                    seqs.append(seq)
        
        if len(seqs) == 0:
            return 0
        else:
            return seqs

def getfile_inference(variant, inpath, substrate, pam):
    geneid = {}
    lines = general(inpath)
    if len(lines) == 0:
        print('No suitable sgRNA was found.')
    else:
        index = 0
        for line in lines:
            index += 1
            #print('each line: ', line.strip())
            line = line.strip().split(',')
            #print(line)
            name = line[0]
            line = line[1].upper()
            #print(name, line)
            rev_line= reverse_seq(line)

            if pam.upper() == 'NGG':
                seqs = get_pam_seq(line, substrate)
                rev_seqs = get_pam_seq(rev_line, substrate)
            elif pam.upper() == 'NNN':
                #print('PAM:', pam)
                seqs = get_pamfree_seq(line, substrate)
                rev_seqs = get_pamfree_seq(rev_line, substrate)
                #print('seqs:',seqs)
                #print('rev_seqs: ',rev_seqs)
            elif pam.upper() == 'FREE':
                if 'GG' in variant:
                    seqs = get_pam_seq(line, substrate)
                    rev_seqs = get_pam_seq(rev_line, substrate)
                else:
                    seqs = get_pamfree_seq(line, substrate)
                    rev_seqs = get_pamfree_seq(rev_line, substrate)
            #print('variant, substrate, pam:', variant, substrate, pam)
            #print('seqs:', seqs)
            #print('rev_seqs:', rev_seqs)

            if seqs == 0 and rev_seqs == 0:
                continue
            elif seqs == 0:
                geneid[name+'||'+str(index)] = [rev_seqs, - len(rev_seqs)]

            elif rev_seqs == 0:
                geneid[name+'||'+str(index)] = [seqs, len(seqs)]
            else:    
                geneid[name+'||'+str(index)] = [seqs + rev_seqs, len(seqs)]
        #print('geneid:\n', geneid)
        if len(geneid) == 0:
            print('There is no NGG PAM in these genomic sequences.')
                
        return geneid

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

def addinfo(dataset):
    '''
    找出变异，加注释，没写完
    '''	
	#mutations
    mutations = []
	
    for idx, row in dataset.iterrows():
        query = np.array(list(row['Predicted-Editing-Outcomes']))
        reference = np.array(list(row['Designed-sgRNA']))
	
        temp = []
        for ref, alt in zip(reference, query):
            temp.append((ref, alt))
		
        var_list = align_mut(temp)
		
		#print('len(var_list):',len(var_list))
        if len(var_list) > 0:
            new_var = '|'.join( x for x in var_list)
        else:
            new_var = 'no-edit'
		
        mutations.append(new_var)

    dataset['mutations'] = mutations
    return dataset

def align_mut(temp):
    var_list = []
    for pos, pair in enumerate(temp):
        ref = pair[0]
        alt = pair[1]
        if pair[0] != alt:
            var = str(pos+1)+':'+ ref + '>' + alt
            var_list.append(var)
        else:
            continue
    return var_list

def get_prop_result(dataset, y_pred):
    #print('dataset:\n',dataset)
    threshold = 0.005
    num2seq, seq2num = generate_combination()
    row, col = np.where(y_pred>threshold)
    locate = pd.DataFrame(columns=['row', 'col'])
    locate['row'] = row
    locate['col'] = col
    #print('y_pred:', len(row), len(col))
    #print('y_pred:', row, col)
    output = pd.DataFrame()
    for i,dfdata in locate.groupby(by='row'):
        temp = pd.DataFrame(columns=['Base-Editor', 'SNP-Site', 'Strand','Designed-sgRNA','pos1-3','pos4-8','pos9-20', 'PAM', 'Predicted-Editing-Outcomes', 'Pred-Proportion'])
        refseq = dataset.at[int(i),'Designed-sgRNA']
        #print('Designed-sgRNA:', refseq)
        pred_prop = list(y_pred[i][dfdata['col']])
        pred_seq = [num2seq[z] for z in dfdata['col'] ]
        #print(pred_seq, len(pred_freq))
        
        ref4to8 = seq2num[refseq[3:8]] #ref seq index    

        temp['pos4-8'] = pred_seq
        temp['Pred-Proportion'] = pred_prop
        temp['Designed-sgRNA'] = refseq
        temp['pos1-3'] = refseq[:3]
        temp['pos9-20'] = refseq[8:]
        temp['Base-Editor'] = dataset.at[int(i),'Base-Editor']
        temp['SNP-Site'] = dataset.at[int(i),'SNP-Site']
        temp['Strand'] = dataset.at[int(i),'Strand']
        temp['PAM'] = dataset.at[int(i),'PAM']

        temp['Predicted-Editing-Outcomes'] = temp['pos1-3'] + temp['pos4-8'] + temp['pos9-20']
        

        if refseq in temp['Predicted-Editing-Outcomes']:
            ref_prop = temp[temp['Predicted-Editing-Outcomes'] == refseq]['Pred-Proportion'].to_list()[0]
            pred_freq = np.array(pred_freq) / (1 - ref_prop)
        
        temp = temp[temp['Predicted-Editing-Outcomes'] != refseq]
        
        output = pd.concat([output, temp[['Base-Editor', 'SNP-Site', 'Strand', 'Designed-sgRNA', 'PAM', 'Predicted-Editing-Outcomes', 'Pred-Proportion']]])
        
    output = output[output['Pred-Proportion'] > threshold]
    
    #output = addinfo(output) #该函数未写完
    #print('output:\n',output)
    #print('locate:',locate.shape)
    return output

def do_predict(baseeditor, genes_info, model_path):
    effpath = model_path[0]
    propath = model_path[1]
    effmodel = load_model(effpath)
    propmodel = load_model(propath)
    effdata = pd.DataFrame(columns=['Base-Editor','SNP-Site','Strand', 'Designed-sgRNA', 'PAM', 'Pred-Efficiency'], dtype=str)
    propdata = pd.DataFrame()
    effdata['Pred-Efficiency'] = effdata['Pred-Efficiency'].astype(float)
    for key, value in genes_info.items():
        plus = value[1]
        tempe = pd.DataFrame(index=range(len(value[0])), columns=['tempseq','Base-Editor', 'SNP-Site', 'Strand', 'Designed-sgRNA', 'PAM', 'Pred-Efficiency'], dtype=str)
        tempe['Pred-Efficiency'] = tempe['Pred-Efficiency'].astype(float)
        tempe['SNP-Site'] = key.split('||')[0]
        tempe['Base-Editor'] = baseeditor
        seq = value[0]
        #print('seq:', seq)
        tempe['tempseq'] = seq
        tempe['Designed-sgRNA'] = tempe['tempseq'].str[:20]
        tempe['PAM'] = tempe['tempseq'].str[20:]
        #print('value[0] :', len(value[0]))
        temp_onehot = encode(value[0])
        tempeff = np.squeeze(effmodel.predict(temp_onehot))
        tempe['Pred-Efficiency'] = tempeff

        #只在负链找到sgRNA
        if plus < 0:
            tempe['Strand'] = '-'

        #只在正链找到sgRNA
        elif plus > 0 and plus == len(value[0]):
            tempe['Strand'] = '+'

        #正、负链均找到sgRNA
        elif plus > 0 and plus < len(value[0]):
            negative = len(value[0])-plus
            strand = ['+']*plus + ['-']*negative
            tempe['Strand'] = strand
        
        del tempe['tempseq']
        effdata = pd.concat([effdata, tempe], axis=0)

        tempprop_pred = propmodel.predict(temp_onehot)
        #print('tempprop_pred', tempprop_pred.shape)
        tempp = get_prop_result(tempe, tempprop_pred)
        propdata = pd.concat([propdata, tempp], axis=0)
        #print('tempp:\n',tempp)
    #print('effdata:',effdata.shape)
    effdata['Pred-Efficiency'].loc[effdata['Pred-Efficiency']<0] = 0.0
    effdata = effdata.sort_values(by='Pred-Efficiency', ascending=False)
    effdata['Pred-Efficiency'] = effdata['Pred-Efficiency'].apply(lambda x:round(x,4))
    propdata['Pred-Proportion'] = propdata['Pred-Proportion'].apply(lambda x:round(x,4))
    #print('effdata:',effdata.shape)
    return effdata, propdata

if __name__ == '__main__':
    path = sys.argv[1]
    #get_snp_seq(path)
    get_offtarget(path)