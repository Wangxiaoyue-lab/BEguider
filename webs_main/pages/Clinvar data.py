import streamlit as st
import pandas as pd

st.set_page_config(page_title="ClinVar data", layout="wide")
st.markdown("## ClinVar data")
st.sidebar.header("ClinVar data")

# Add prompt
st.markdown(
    """
    ##### _Note: All sites are set at the sixth position of the sgRNA._
    - **'REF'** means that you could try to **repair the pathogenic site.**
    - **'ALT'** means that you could try to **generate the disease model.**
    """
)


# Create two columns
col1, col2 = st.columns(2)

# Place select boxes in columns
type1 = col1.selectbox("Base editors", ["ABE", "CBE"])
type2 = col2.selectbox("Allele alteration", ["REF", "ALT"])


# Add explanatory text

# Determine the file name to read based on the selected type
filename = f"{type1[0].lower()}_{type2.lower()}_anno.csv"
filepath = f"../datas/{filename}"


# csv
@st.cache_data
def load_data(filepath):
    df = pd.read_csv(filepath)
    return df


@st.cache_data
def convert_df(df):
    return df.to_csv().encode("utf-8")


df = load_data(filepath)
col_names = df.columns[-9:]
num_col = [x for i, x in enumerate(col_names, 1) if i not in [3, 6, 9]]
edit_col = col_names[[3, 5, 8]]

#df = df[list(df.columns[-9:]) + list(df.columns[:-9])]  # exchange the order

st.download_button(
    label="Download all the data for local query",
    data=convert_df(df),
    file_name=filename,
    mime="text/csv",
)

df = df.drop(
    [
        "index",
        "lib",
        "base",
        "ref_position",
        "significance",
    ],
    axis=1,
)

st.markdown(
    """
    ## Online query
    """
)


(
    col_0,
    col_1,
    col_2,
    col_3,
) = st.columns(4)


# Create multiselect widgets for the chrom and site columns
site_values = col_0.multiselect("Site", df["site"].unique())
symbol = col_1.multiselect("Symbol", df["symbol"].unique())
options = set()
for value in df["caused_disease"].unique():
    options.update(value.split("|"))
disease = col_2.multiselect("Disease", options)
editable = col_3.multiselect("Editable", ["Y", "N"])

def highlight_Y(s):
    return ['background-color: green' if v == 'Y' else '' for v in s]
def format_number(val):
    return f"{val:.4f}"


cols = st.columns(8)
# Filter data frame based on selected values
if cols[0].button("Update"):
    # Create a mask to filter rows in the data frame
    mask = pd.Series(True, index=df.index)

    # Apply filters for columns with selected values
    if symbol:
        mask &= df["symbol"].isin(symbol)
    if site_values:
        mask &= df["site"].isin(site_values)
    if disease:
        mask &= df["caused_disease"].apply(lambda x: any(d in x for d in disease))
    if editable:
        mask &= (
            df[col_names[2]].isin(editable)
            | df[col_names[5]].isin(editable)
            | df[col_names[8]].isin(editable)
        )

    # Filter data frame based on mask
    filtered_df = df[mask]

    # Display filtered data frame or message if no data found
    if filtered_df.empty:
        st.write("Sorry, no data were found.")
    else:
        st.dataframe(filtered_df.style.apply(highlight_Y, subset=edit_col).format({col: format_number for col in num_col}))
        cols[1].download_button(
            label="Download",
            data=convert_df(filtered_df),
            file_name="filtered_data.csv",
            mime="text/csv",
        )


st.markdown(
    """
    #####  \* Click the **Update** button to check the results based on the base editors and your proposals. 
    Sites that were defined as editable are highlighted.
    - Some concepts:
        - Base editing **efficiency** refers to the proportion of all edited outcomes in the editing window to all outcomes.
        - Base editing outcome **proportion** (pred_prop) refers to the predicted proportion of the edited substrate 
        (i.e., the proportion of the edited sixth base to all edited outcomes in the editing window).
        - Absolute **frequency** (pred_freq) shows the absolute predicted frequency of the edited substrate 
        (i.e., the proportion of the edited sixth base to all outcomes).
        - The threshold for editable sites is defined as **outcome proportion >= 90% & efficiency >= 5%**.
    """
)
