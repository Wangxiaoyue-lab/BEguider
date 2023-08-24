import streamlit as st
import pandas as pd
import os
import shutil
from pathlib import Path

script_dir = os.path.dirname(__file__)


# Set page title
st.set_page_config(page_title="BEguider", layout="wide")


def local_css(file_name):
    with open(file_name) as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)


local_css(script_dir + "/style/style.css")

st.sidebar.header("Intro")

# Add page title
st.write("# Welcome to BEguider!")
st.markdown(
    """
    **BEguider is a deep-learning model, could design sgRNA for multiple gene loci in batches to predict the efficiency 
    and products of each sgRNA.**

    ## BEguider
    - Base editors with Strict Limits of NGG PAM: `ABE7.10-NGG`, `BE4-NGG`
    - PAM-less base editors: `ABEmax-SpRY`, `ABE8e-P(AP)3-SpRY`, `ABE8e-NL-SpRY`,`BE4max-SpRY`, `FNLS-YE1-BE4max-SpRY`, `YE1-BE4max-SpRY`
    - The pairs of Base Editors and SNVs: 
        - For **ABE**: _pos_  :orange[A] > :green[G] ,  _neg_  :blue[T] > :red[C] 
        - For **CBE**: _pos_  :red[C] > :blue[T] ,  _neg_  :green[G] > :orange[A]       
        
    - Jump to our [code repository](https://github.com/Wangxiaoyue-lab/BEguider) if you want to run model locally or get help!

    ## Clinvar data
    - Data from clinvar were predicted by PAM-less base editor for reference. (Clinvar edition:2021/03/02)
    - Choose more sites from [clinvar database](https://www.ncbi.nlm.nih.gov/clinvar/)
    
    ## About
    - Thanks for your using.
    - Please cite our paper if you find our model useful.
        - To be published.
"""
)


def get_size(path):
    # Calculates the total size of all files in the specified path
    total_size = 0
    for entry in Path(path).rglob("*"):
        if entry.is_file():
            total_size += entry.stat().st_size
    return total_size


def monitor_and_cleanup_folder(path, max_folders, num_folders_to_delete, max_size):
    # Gets all subfolders in the specified path
    folders = [entry for entry in Path(path).iterdir() if entry.is_dir()]

    # If the number of subfolders exceeds the threshold
    if len(folders) > max_folders:
        # Sort by creation time

        folders.sort(key=lambda folder: folder.stat().st_ctime)
        for folder in folders[:num_folders_to_delete]:
            shutil.rmtree(folder)

    # Calculates the total size of all subfolders in the specified path
    total_size = get_size(path)
    # If the total size exceeds the threshold

    if total_size > max_size * 1024 * 1024 * 1024:
        # sort by size
        folders.sort(key=lambda folder: get_size(folder), reverse=True)

        # Delete the maximum 10 subfolders
        for folder in folders[:10]:
            shutil.rmtree(folder)


max_folders = 2000
num_folders_to_delete = 500
max_size = 5  # G
# Monitor specified folder


path = script_dir + "/temp_dirs"
monitor_and_cleanup_folder(path, max_folders, num_folders_to_delete, max_size)
