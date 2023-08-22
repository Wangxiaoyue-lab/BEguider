import pandas as pd
import sys
import os

def main(effpath, proppath, outpath):
    effdata = pd.read_csv(effpath)
    propdata = pd.read_csv(proppath)

    effdata = effdata[effdata['Strand'] == '+']
    propdata = propdata[propdata['Strand'] == '+']

    mergedata = pd.merge(effdata, propdata, how='outer', on=['Base-Editor', 'SNP-Site', 'Strand', 'Designed-sgRNA', 'PAM'])
    mergedata['pred_abs_freq'] = mergedata['Pred-Proportion'] * mergedata['Pred-Efficiency']
    #print('mergedata:\n', mergedata)
    mergedata.to_csv(outpath, index=False)

if __name__ == '__main__':
    effpath = sys.argv[1]
    proppath = sys.argv[2]
    outpath = sys.argv[3]
    main(effpath, proppath, outpath)