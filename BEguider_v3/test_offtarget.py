import os
import sys
import pandas as pd
from subprocess import Popen
import argparse
import platform
import itertools

def get_parser():
    desc = """
    Program: BEguider
    Version: 0.1
    
           """

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=desc)
    parser.add_argument('-i','--input',help="a csv")
    parser.add_argument('-b','--BaseEditor', type=str, required = True, help="ABE / CBE / ALL. Design sgRNA for ABE(ABE7.10-NGG) or CBE(BE4-NGG) or ALL(ABE7.10-NGG, BE4-NGG)" )
    parser.add_argument('-f','--offtarget', type=bool, default=False, help="True / False. True: predict offtarget sites in hg38 genome. (Default=False)" )
    parser.add_argument('-m','--mismatch',type=int, default=4, help='Allowed maximum mismatch site between sgRNA and genome while finding off-target sites.(Default=4)')
    parser.add_argument('-o','--output', type=str, default=None, help="output directory (Default = current directory)" )
        
    return parser.parse_args()

def process_dir(BaseEditor, seq_path, outpath):
    if outpath is not None:
        cwd = outpath
    else:
        cwd = os.getcwd()
    
    if (platform.system()=='Windows'):
        name = seq_path.split('\\')[-1].split('.')[0]
        outfile = cwd + '\\outputs\\' + name +  '_designed_sgRNA_for'+'_'+BaseEditor+'.csv'
        outpath = cwd + '\\outputs\\'
        tempdir = cwd+ '\\outputs\\temp'
        offpath = tempdir + '\\seqs_to_offtarget.txt'
    else:
        name = seq_path.split('/')[-1].split('.')[0]
        outfile = cwd + '/outputs/' + name +  '_designed_sgRNA_for'+'_'+BaseEditor+'.csv'
        outpath = cwd + '/outputs/'
        tempdir = cwd+ '/outputs/temp'
        offpath = tempdir + '/seqs_to_offtarget.txt'

    if "outputs" not in os.listdir(cwd):
        os.makedirs(outpath)
    
    if "temp" not in os.listdir(outpath):
        os.makedirs(tempdir)

    off_txt_header(offpath)

    return outfile, tempdir, offpath

def off_txt_header(offpath):
    
    if '\\' in offpath:
        hg38 = '..\\Util\\'
    else:
        hg38 = '../Util/'
    
    pattern = 'NNNNNNNNNNNNNNNNNNNNNGG'

    with open(offpath, 'w') as files:
        files.write(hg38+'\n')
        files.write(pattern+'\n')

def offtarget(tempdir):
    dirs =  tempdir[:-4]
    inpath = tempdir + '/seqs_to_offtarget.txt'
    outpath = dirs +'/offtarget_result.txt'

    if (platform.system()=='Windows'):
        inpath = tempdir + '\\seqs_to_offtarget.txt'
        outpath = dirs +'\\offtarget_result.txt'
        cmd = '..\\Util\\win-cas-offinder.exe ' + inpath + ' C' + ' ' + outpath
        print('Windows:', cmd)
    elif (platform.system()=='Linux'):
        #cmd = '../Util/linux-cas-offinder ' + inpath + ' C' + ' ' + outpath
        cmd = 'cas-offinder ' + inpath + ' C' + ' ' + outpath
        print('Linux:', cmd)
    elif (platform.system()=='Darwin'):
        cmd = '../Util/mac-cas-offinder ' + inpath + ' C' + ' ' + outpath
        #print('MacOS:', cmd)
    else:
        print('The system edition does not support for finding off-target sites.')
        print('Please use Windows/Linux/MacOS platform.')
        sys.exit()

    off = Popen(args=cmd, shell=True)
    off.wait() #等待进程完成
    #print('off.returncode:' ,off.returncode)
    if off.returncode == 0:
        print('Off-target finished.')

def main(BaseEditor, parse=None):

    outdata = pd.read_csv(parse.input, header=0)
    inp = '/home/gjj/software/jupyter/BEguider_model/New/test_chrom.txt'
    outfile, tempdir, offpath = process_dir(BaseEditor, inp, parse.output)

    if parse.offtarget:
        outdata['temp1'] = outdata['Designed-sgRNA'].str[:-3]
        outdata['temp2'] = 'NNN'
        outdata['mis-match'] = parse.mismatch
        outdata['cas-offinder_seq'] = outdata['temp1'] + outdata['temp2']
        to_off = outdata[['SNP-Site', 'cas-offinder_seq', 'mis-match']]
        #to_off['mis-match'] = 4
        to_off[['cas-offinder_seq', 'mis-match']].to_csv(offpath, mode='a', sep=' ', index=False, header=False)
        
        offtarget(tempdir)


if __name__ == '__main__':
    
    parse = get_parser()
    
    BaseEditor = parse.BaseEditor
    BaseEditor = BaseEditor.upper()
    assert BaseEditor in ['ABE', 'CBE', 'ALL']

    main(BaseEditor=BaseEditor, parse=parse)