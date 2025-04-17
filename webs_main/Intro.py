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
    **BEguider is a deep-learning model that can design sgRNA for multiple gene 
    loci in batches to predict the efficiency and products of different base editors.**

    ## BEguider
    - Base editors with strict limits of NGG PAM: `ABE7.10-NGG`, `BE4-NGG`
    - PAM-less base editors: `ABEmax-SpRY`, `ABE8e-SL-SpRY`, `ABE8e-NL-SpRY`,`BE4max-SpRY`, `FNLS-YE1-SpRY`, `YE1-SpRY`
    - The pairs of base editors and SNVs: 
        - For **ABE**: _pos_  :orange[A] > :green[G] ,  _neg_  :blue[T] > :red[C] 
        - For **CBE**: _pos_  :red[C] > :blue[T] ,  _neg_  :green[G] > :orange[A]       
    - In addition to directly entering sequence information, [hg38 chromosome loci](http://genome.ucsc.edu/) or [rsID](https://www.ncbi.nlm.nih.gov/snp/) could also be supported.
    

    ## ClinVar data
    - Data from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) had already been predicted by PAM-less base editor for reference (ClinVar edition: 2021/03/02).
    
    
    ## About
    - Thank you for using our model.
    - Jump to our [code repository](https://github.com/Wangxiaoyue-lab/BEguider) if you want to run the model locally or get help!
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


####################################################################
import sqlite3


# Initialize the SQLite database connection
conn = sqlite3.connect("view_count.db")
c = conn.cursor()

# Create a table to store the view count
c.execute(
    """
    CREATE TABLE IF NOT EXISTS ViewCount (
        id INTEGER PRIMARY KEY,
        count INTEGER
    )
"""
)

# Initialize the count if it's the first run
c.execute("SELECT count FROM ViewCount WHERE id=1")
row = c.fetchone()
if row is None:
    c.execute("INSERT INTO ViewCount (id, count) VALUES (1, 0)")
    conn.commit()


# Function to update and retrieve the view count
def update_view_count():
    c.execute("UPDATE ViewCount SET count = count + 1 WHERE id=1")
    conn.commit()
    c.execute("SELECT count FROM ViewCount WHERE id=1")
    return c.fetchone()[0]


# Increment the view count and get the updated value
view_count = update_view_count()


# CSS animation for "+1" effect
st.markdown(
    """
<style>
@keyframes slideUp {
  from {opacity: 1; transform: translateY(0);}
  to {opacity: 0; transform: translateY(-20px);}
}
.count {
  font-size: 16px;
  color: #000000;
  position: relative;
  display: inline-block;
}
.plus-one {
  font-size: 14px;
  color: #000000;
  position: relative;
  display: inline-block;
  margin-left: 5px;
  animation: slideUp 1s ease-out forwards;
  pointer-events: none;
}
</style>
""",
    unsafe_allow_html=True,
)

# Display the view count with "+1" animation
st.markdown(
    f"""
<div class="count">
  Page View Count: {view_count}
  <span class="plus-one">+1</span>
</div>
""",
    unsafe_allow_html=True,
)
