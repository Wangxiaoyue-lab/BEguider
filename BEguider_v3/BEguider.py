#!/usr/bin/env python3
# -*- coding utf-8 -*-

import pandas as pd
import argparse
import shutil
import platform
import sys
import os

add = "./Util"
if platform.system() == "Windows":
    add = ".\\Util"
sys.path.append(add)
# from Utilities import *
from Utilities import *


def get_parser():
    desc = """
    Program: BEguider
    Version: 0.2
    Author : Jingjing Gao
    Email  : <gaojingjing@ibms.pumc.edu.cn>

    Optional Base Editors: ['ABEmax-SpRY', 'ABE8e-SL-SpRY', 'ABE8e-NL-SpRY', 'BE4max-SpRY', 
                            'FNLS-YE1-SpRY', 'YE1-SpRY', 'ABE7.10-NGG', 'BE4-NGG']
    
    Base Editors with Strict Limits of NGG PAM: ['ABE7.10-NGG', 'BE4-NGG']
    PAM-less Base Editors: ['ABEmax-SpRY', 'ABE8e-SL-SpRY', 'ABE8e-NL-SpRY', 
                            'BE4max-SpRY', 'FNLS-YE1-SpRY', 'YE1-SpRY'] 

    The pairings of Base Editors and SNPs:
    CBE -- C>T ; G>A       ABE -- A>G ; T>C

    Examples:
        1) use genes as input file:
                Genes,Seqs
                TP53,CTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCT
                BRCA1,TGGCTGAAGAATTTGCTAAGCAATCAGGAAAGCTGGTGG
                ARID1A,CTCCACCGAAGGGAGGACCCACTGCCCCCAGCCGGGGTCTCG

        2) use chromosomes and coordinates as input file:
                Chrom,Coordinate,Type
                chr1,145634,r
                chrX,87632,a

        3) use rsIDs as input file:
                SNP,Type
                rs5297,r
                rs12603332,a
                rs75456785,r

           """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=desc
    )
    igroup = parser.add_mutually_exclusive_group(required=True)
    igroup.add_argument(
        "-g",
        "--genes",
        help="A txt file including genes and sequences separated by comma(,).",
    )
    igroup.add_argument(
        "-c",
        "--chromosome",
        help="A txt file including chromosomes, coordinates and genetic type separated by comma(,). Genetic type: 'r' means editing wild genes, 'a' means editing mutant genes.",
    )
    igroup.add_argument(
        "-s",
        "--rsID",
        help="A txt file including rsID and genetic type separated by comma(,).",
    )
    parser.add_argument(
        "-b",
        "--BaseEditor",
        type=str,
        required=True,
        help="Base Editors: ALL / A Specific BE name. Design sgRNA for a specific BE in optional BEs or for all BEs",
    )
    parser.add_argument(
        "-f",
        "--offtarget",
        type=str,
        default=False,
        help="True / False. True: predict off-target sites in hg38 genome. (Default = False)",
    )
    parser.add_argument(
        "-m",
        "--mismatch",
        type=int,
        default=3,
        help="Allowed maximum mismatch site between sgRNA and genome while searching for off-target sites.(Default = 3)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="./",
        help="Output directory. (Default = current directory)",
    )

    return parser.parse_args()


def get_predict_data(variant, pam, inp):
    if "ABE" in variant:
        substrate = "A"
    else:
        substrate = "C"

    genes_info = getfile_inference(variant, inp, substrate, pam)
    # print('substrate, genes_info:\n', substrate, genes_info)
    model_path = GetParams(variant)
    # print(model_path)
    if platform.system() == "Windows":
        model_path = model_path.replace("/", "\\")

    if len(genes_info) > 0:
        effout, propout = do_predict(variant, genes_info, model_path)
    else:
        effout = pd.DataFrame()
        propout = pd.DataFrame()
    return effout, propout


def get_PAM(BaseEditor):
    nnnnbe = [
        "ABEmax-SpRY",
        "ABE8e-SL-SpRY",
        "ABE8e-NL-SpRY",
        "BE4max-SpRY",
        "FNLS-YE1-SpRY",
        "YE1-SpRY",
    ]
    nggbe = ["ABE7.10-NGG", "BE4-NGG"]
    if BaseEditor == "ALL":
        pam = "FREE"
    elif BaseEditor in nggbe:
        pam = "NGG"
    elif BaseEditor in nnnnbe:
        pam = "NNN"
    return pam


def choose_variants(belist, pam, inp):
    alleff = pd.DataFrame()
    allprop = pd.DataFrame()
    for variant in belist:
        # print('variant:',variant)
        tempeff, tempprop = get_predict_data(variant, pam, inp)
        alleff = pd.concat([alleff, tempeff])
        allprop = pd.concat([allprop, tempprop])
    return alleff, allprop


def main(BaseEditor, parse=None):
    genes = parse.genes
    chrom = parse.chromosome
    rs = parse.rsID
    pam = get_PAM(BaseEditor)
    # print('PAM:', pam)
    print("genes:", genes, " chrom:", chrom, " rs:", rs)
    # print('pam:',pam)
    allbe = [
        "ABEmax-SpRY",
        "ABE8e-SL-SpRY",
        "ABE8e-NL-SpRY",
        "BE4max-SpRY",
        "FNLS-YE1-SpRY",
        "YE1-SpRY",
        "ABE7.10-NGG",
        "BE4-NGG",
    ]
    abes = ["ABEmax-SpRY", "ABE8e-SL-SpRY", "ABE8e-NL-SpRY", "ABE7.10-NGG"]
    cbes = ["BE4max-SpRY", "FNLS-YE1-SpRY", "YE1-SpRY", "BE4-NGG"]

    # inp[0], tags[0] -- ABE
    # inp[1], tags[1] -- CBE
    if genes:
        inp = [genes, genes]
        outfile, tempdir, offpath, fname = process_dir(BaseEditor, genes, parse.output)
        tags = [2, 2]

    elif chrom:
        outfile, tempdir, offpath, fname = process_dir(BaseEditor, chrom, parse.output)
        inp, tags = get_genomic_seq(pam, chrom, tempdir)
        print(inp, tags)

    elif rs:
        outfile, tempdir, offpath, fname = process_dir(BaseEditor, rs, parse.output)
        inp, tags = get_snp_seq(pam, rs, tempdir)

    # print('BaseEditor:',BaseEditor,' tag:',tags)
    if BaseEditor == "ALL":
        if tags[0] == 2:
            # print('all tagABE, tagCBE')
            alleff, allprop = choose_variants(allbe, pam, inp[0])

        elif tags[0] == 1 and tags[1] == 1:
            # print('1 tagABE, tagCBE')
            abe_eff, abe_prop = choose_variants(abes, pam, inp[0])
            cbe_eff, cbe_prop = choose_variants(cbes, pam, inp[1])

            alleff = pd.concat([abe_eff, cbe_eff])
            allprop = pd.concat([abe_prop, cbe_prop])
            if alleff.shape[0] == 0:
                print("No suitable BaseEditor.")
                sys.exit()
        elif tags[0] == 0 and tags[1] == 1:
            alleff, allprop = choose_variants(cbes, pam, inp[1])
        elif tags[0] == 1 and tags[1] == 0:
            # print('tags[0]==1 and tags[1]==0')
            alleff, allprop = choose_variants(abes, pam, inp[0])
    elif BaseEditor != "ALL":
        if tags[0] == 2:
            alleff, allprop = get_predict_data(BaseEditor, pam, inp[0])
        elif tags[0] != 2:
            if (BaseEditor in abes and tags[0] == 0) or (
                BaseEditor in cbes and tags[1] == 0
            ):
                print("No suitable Base Editors.")
                sys.exit()
            elif BaseEditor in abes :#and tags[0] == 1:
                alleff, allprop = get_predict_data(BaseEditor, pam, inp[0])
            elif BaseEditor in cbes:# and tags[0] == 1:
                alleff, allprop = get_predict_data(BaseEditor, pam, inp[1])
    elif BaseEditor not in allbe:
        print("The base editor was wrong.")
        print("Please use command: python BEguider -h")
        sys.exit()
    alleff = alleff.drop_duplicates()
    alleff.to_csv(outfile[0], header=True, index=False)
    allprop = allprop.drop_duplicates(
        subset=["Base-Editor", "SNP-Site", "Predicted-Editing-Outcomes"]
    )
    allprop.to_csv(outfile[1], header=True, index=False)

    tag = False
    if str(parse.offtarget).lower() == "false":
        tag = False
    elif str(parse.offtarget).lower() == "true":
        tag = True

    if tag:
        alleff["temp1"] = alleff["Designed-sgRNA"]
        alleff["temp2"] = "NNN"
        alleff["mismatch"] = parse.mismatch
        alleff["cas-offinder_seq"] = alleff["temp1"] + alleff["temp2"]
        to_off = alleff[["SNP-Site", "cas-offinder_seq", "mismatch"]]
        to_off[["cas-offinder_seq", "mismatch"]].to_csv(
            offpath, mode="a", sep=" ", index=False, header=False
        )

        offoutpath = get_offtarget(fname, tempdir)
        # offinfo = pd.read_csv(offoutpath, header=None, sep='\t')
        # offinfo = offinfo.rename(columns={0:'ref-seq',1:'chrom',2:'coordinate',3:'offtarget-seq',4:'strand',5:'mismatch-num'})
    else:
        print("")

    shutil.rmtree(tempdir)


if __name__ == "__main__":
    parse = get_parser()

    BaseEditor = parse.BaseEditor
    if BaseEditor.upper() == "ALL":
        BaseEditor = "ALL"
    else:
        BaseEditor = BaseEditor

    assert BaseEditor in [
        "ALL",
        "ABEmax-SpRY",
        "ABE8e-SL-SpRY",
        "ABE8e-NL-SpRY",
        "BE4max-SpRY",
        "FNLS-YE1-SpRY",
        "YE1-SpRY",
        "ABE7.10-NGG",
        "BE4-NGG",
    ]

    main(BaseEditor=BaseEditor, parse=parse)
