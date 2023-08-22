from curses import beep
import streamlit as st
import pandas as pd

st.set_page_config(page_title="BEguider", layout="wide")
st.markdown("# BEguider")

st.sidebar.header("BEguider")

import sys

sys.path.insert(0, "../BEguider_v3")
import BEguider_web as bw

st.markdown("## Input\n**Examples:**\n")

tab1, tab2, tab3 = st.tabs(["genes", "chromosomes", "rsIDs"])


with tab1:
    st.markdown(
        """
    ```
    Genes,Seqs
    TP53,CTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCT
    ```
    """
    )
    filename = "Gene_example.txt"
    st.download_button(
        label="Download sample data",
        data=f"../BEguider_v3/examples/{filename}",
        file_name=filename,
        mime="text/csv",
    )

with tab2:
    st.markdown(
        """
    ```
    Chrom,Coordinate,Type
    chr1,145634,r
    ```
    """
    )
    filename = "Chrom_example.txt"
    st.download_button(
        label="Download sample data",
        data=f"../BEguider_v3/examples/{filename}",
        file_name=filename,
        mime="text/csv",
    )

with tab3:
    st.markdown(
        """
    ```
    SNP,Type
    rs5297,r
    ```
    """
    )
    filename = "SNP_example.txt"
    st.download_button(
        label="Download sample data",
        data=f"../BEguider_v3/examples/{filename}",
        file_name=filename,
        mime="text/csv",
    )


BaseEditor = [
    "ALL",
    "ABEmax-SpRY",
    "ABE8e-P(AP)3-SpRY",
    "ABE8e-NL-SpRY",
    "BE4max-SpRY",
    "FNLS-YE1-BE4max-SpRY",
    "YE1-BE4max-SpRY",
    "ABE7.10-NGG",
    "BE4-NGG",
]


class Parse:
    def __init__(
        self,
        genes=None,
        chromosome=None,
        rsID=None,
        offtarget=False,
        mismatch=3,
        output="./",
    ):
        self.genes = genes
        self.chromosome = chromosome
        self.rsID = rsID
        self.offtarget = offtarget
        self.mismatch = mismatch
        self.output = output


parse = Parse(genes=gn, chromosome=cs, rsID=ri, offtarget=ot, mismatch=mm)


alleff, allprop = bw.main_pre(be, parse=parse)
