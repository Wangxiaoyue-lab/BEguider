from curses import beep
import streamlit as st
import pandas as pd
import uuid
import os

session_id = str(uuid.uuid4())
session_dir = os.path.join("temp_dirs", session_id)
if not os.path.exists(session_dir):
    os.makedirs(session_dir)

st.set_page_config(page_title="BEguider", layout="wide")
st.markdown("# BEguider")

st.sidebar.header("BEguider")

import sys

sys.path.insert(0, "../BEguider_v3")
import BEguider_web as bw

st.markdown("\n**Examples:**\n")

tab1, tab2, tab3 = st.tabs(["genes", "chromosomes", "rsIDs"])


@st.cache_data
def load_data(filepath, **kwargs):
    df = pd.read_csv(filepath, **kwargs)
    return df


@st.cache_data
def convert_df(df):
    return df.to_csv(index=False).encode("utf-8")


def ld_sp():
    st.download_button(
        label="Download sample .csv",
        data=convert_df(load_data(f"../BEguider_v3/examples/{filename}")),
        file_name=filename,
        mime="text/csv",
    )


with tab1:
    st.markdown(
        """
    ```csv
    Genes,Seqs # Including genes and sequences separated by comma(,)
    TP53,CTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCT
    ```
    """
    )
    filename = "Gene_example.csv"
    ld_sp()

with tab2:
    st.markdown(
        """
    ```csv
    Chrom,Coordinate,Type # Including chromosomes, coordinates and genetic type separated by comma(,)
    chr1,145634,r # Genetic type: 'r' means editing wild genes, 'a' means editing mutant genes.
    ```
    """
    )
    filename = "Chrom_example.csv"
    ld_sp()

with tab3:
    st.markdown(
        """
    ```csv
    SNP,Type # Including rsID and genetic type separated by comma(,)
    rs5297,r  # Genetic type: 'r' means editing wild genes, 'a' means editing mutant genes.
    ```
    """
    )
    filename = "SNP_example.csv"
    ld_sp()


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


option = st.radio(
    "Please choose the input types. The uploaded file must have a header like the examples provided.",
    ("genes", "chromosomes", "rsIDs"),
)
col01, col02 = st.columns(2)
# gn, cs, ri = None, None, None


if option == "genes":
    st.session_state.gn_value = "Genes,Seqs\n"
    gn = col01.text_area(
        "Including genes and sequences separated by comma(,)",
        value=st.session_state.get("gn_value"),
    ).strip()
    st.session_state.gn_value = gn

elif option == "chromosomes":
    st.session_state.cs_value = "Chrom,Coordinate,Type\n"
    cs = col01.text_area(
        "Including chromosomes, coordinates and genetic type separated by comma(,)",
        value=st.session_state.get("cs_value"),
    ).strip()
    st.session_state.cs_value = cs

elif option == "rsIDs":
    st.session_state.ri_value = "SNP,Type\n"
    ri = col01.text_area(
        "Including rsID and genetic type separated by comma(,)",
        value=st.session_state.get(
            "ri_value",
        ),
    ).strip()
    st.session_state.ri_value = ri

uploaded_file = col02.file_uploader(
    "**Or** choose a .csv file. Text box has a higher priority.",
    type=["csv"],
)

with st.form(key="my_form"):
    BaseEditor = [
        "ABEmax-SpRY",
        "ABE8e-SL-SpRY",
        "ABE8e-NL-SpRY",
        "BE4max-SpRY",
        "FNLS-YE1-SpRY",
        "YE1-SpRY",
        "ABE7.10-NGG",
        "BE4-NGG",
        "ALL",
    ]

    (
        col1,
        col2,
        col3,
    ) = st.columns(3)

    be = col1.selectbox("Please choose an editor", BaseEditor)
    # ot = col2.selectbox(
    #     "Predict off-target sites in hg38 genome",
    #     [
    #         False,
    #         True,
    #     ],
    # )
    # mm = col3.selectbox(
    #     "Allowed maximum mismatch site",
    #     [
    #         0,
    #         1,
    #         2,
    #         3,
    #         4,
    #         5,
    #     ],
    # )
    ot = False
    mm = 0
    submit_button = st.form_submit_button(label="Submit")

if submit_button:
    file_path = None
    gn_path, cs_path, ri_path = None, None, None

    if uploaded_file is not None:
        file_path = os.path.join(session_dir, uploaded_file.name)
        with open(file_path, "wb") as f:
            f.write(uploaded_file.getvalue())
        st.info("Your file has been uploaded.", icon="ℹ️")

    if option == "genes":
        if gn != "Genes,Seqs":
            # 写入genes.csv
            file_path = os.path.join(session_dir, "genes.csv")
            with open(file_path, "w") as f:
                f.write(gn)
        gn_path = file_path

    elif option == "chromosomes":
        if cs != "Chrom,Coordinate,Type":
            # 写入chromosomes.csv
            file_path = os.path.join(session_dir, "chromosomes.csv")
            with open(file_path, "w") as f:
                f.write(cs)
        cs_path = file_path

    elif option == "rsIDs":
        if ri != "SNP,Type":
            # 写入rsIDs.csv
            file_path = os.path.join(session_dir, "rsIDs.csv")
            with open(file_path, "w") as f:
                f.write(ri)
        ri_path = file_path
    if gn_path is None and cs_path is None and ri_path is None:
        st.warning("Your input is empty, please check your input!")
        st.stop()

    parse = Parse(
        genes=gn_path,
        chromosome=cs_path,
        rsID=ri_path,
        offtarget=ot,
        mismatch=mm,
        output=session_dir,
    )

    try:
        with st.spinner("Running..."):
            alleff, allprop = bw.main_pre(be, parse=parse)
            all_files= pd.merge(alleff,allprop, how="inner",on=["Base-Editor","SNP-Site","Strand","Designed-sgRNA","PAM"])
            all_files=all_files[all_files["Pred-Efficiency"]!=0]
            all_files['Pred-Frequency']=all_files["Pred-Efficiency"] * all_files["Pred-Proportion"]
    except SystemExit as e:
        st.error("Sorry, no suitable Base Editors.")
    except Exception as e:
        if str(e) in (
            "local variable 'new_seq' referenced before assignment",
            "4",
            "R",
        ):
            st.error("Please check your input format!", icon="🚨")
        else:
            st.error(f"ERROR: {e}", icon="🚨")

    else:
        st.success("SUCCESS!")
        
        col11, col12, _, _, _ = st.columns(5)

        col11.download_button(
            label="Download results",
            data=convert_df(all_files),
            file_name="Frequency.csv",
            mime="text/csv",
        )

        # col12.download_button(
        #     label="Download **proption** files",
        #     data=convert_df(allprop),
        #     file_name="proption.csv",
        #     mime="text/csv",
        # )

st.markdown(
    """
- Some concepts:
    - Base editing **efficiency** refers to the proportion of all edited outcomes in the editing window to all outcomes.
    - Base editing outcome **proportion** refers to the predicted proportion of the edited substrate 
    (i.e., the proportion of the edited substrate base to all edited outcomes in the editing window).
    - Absolute **frequency** shows the absolute predicted frequency of the edited substrate 
    (i.e., the proportion of the edited substrate base to all outcomes).
"""
)
