import streamlit as st
import pandas as pd

# 设置页面标题
st.set_page_config(page_title="BEguider", layout="wide")
st.sidebar.header("Intro")

# 添加页面标题
st.write("# Welcome to BEguider!")
st.markdown(
    """
    **BEguider is a deep-learning model, can design sgRNA for multiple gene loci in batches to predict the editing efficiency 
    and editing products of each sgRNA.**

    ## BEguider
    - Base editors with Strict Limits of NGG PAM: `ABE7.10-NGG`, `BE4-NGG`
    - PAM-less base editors: `ABEmax-SpRY`, `ABE8e-P(AP)3-SpRY`, `ABE8e-NL-SpRY`,`BE4max-SpRY`, `FNLS-YE1-BE4max-SpRY`, `YE1-BE4max-SpRY`
    - The pairs of Base Editors and SNVs: 
        - For **ABE**: _pos_  :orange[A] > :green[G] ,  _neg_  :blue[T] > :red[C] 
        - For **CBE**: _pos_  :red[C] > :blue[T] ,  _neg_  :green[G] > :orange[A]       
        
    - Jump into our [source code](https://docs.streamlit.io) if you want to run model locally or get help!

    ## Clinvar data
    - Data from clinvar has been predicted by PAM-less base editors and can be used for reference. (Clinvar edition:2021/03/02)
    - Choose more sites from [clinvar database](https://www.ncbi.nlm.nih.gov/clinvar/)
    
    ## About
    - Thanks for your using.
    - Please cite our paper if you find our model useful.
        - To be published.
"""
)
