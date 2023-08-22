import streamlit as st
import pandas as pd

st.set_page_config(page_title="Clinvar data", layout="wide")
st.markdown("## Clinvar data")
st.sidebar.header("Clinvar data")

# 添加提示信息
st.markdown(
    """
    ##### _Note: All sites were set at the sixth position of the sgrna._
    - **'REF'** means that you could try to **repair the pathogenic site.**
    - **'ALT'** means that you could try to **generate the disease model.**
    """
)


# Create two columns
col1, col2 = st.columns(2)

# Place select boxes in columns
type1 = col1.selectbox("Base editors", ["ABE", "CBE"])
type2 = col2.selectbox("Proposes", ["REF", "ALT"])


# Add explanatory text

# 根据选择的类型确定要读取的文件名
filename = f"{type1[0].lower()}_{type2.lower()}_anno.csv"
filepath = f"../datas/{filename}"

# Create download button
st.download_button(
    label="Download all the data for local query",
    data=filepath,
    file_name=filename,
    mime="text/csv",
)


st.markdown(
    """
    ## Online query
    """
)


# 读取csv文件
@st.cache_data
def load_data(filepath):
    df = pd.read_csv(filepath)
    return df.drop(
        [
            "index",
            "lib",
            "base",
            "ref_position",
            "significance",
            "20nt_sgRNA",
            "4nt_PAM",
        ],
        axis=1,
    )


df = load_data(filepath)
col_names = df.columns[-9:]
result = [x for i, x in enumerate(col_names, 1) if i not in [3, 6, 9]]
df = df[list(df.columns[-9:]) + list(df.columns[:-9])]  # 调换顺序


def highlight_greater_than_threshold(val, threshold):
    color = "yellow" if val > threshold else ""
    return f"background-color: {color}"


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


cols = st.columns(10)
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
        st.dataframe(filtered_df.style.highlight_max(subset=result, axis=0))
        cols[1].download_button(
            label="Download",
            data=filtered_df.to_csv(index=False),
            file_name="filtered_data.csv",
            mime="text/csv",
        )

st.markdown(
    """
    #####  \* Click the **Update** button to check the results based on the base editors and your proposes. Sites that performed best on each prediction were highlighted.
    #####  \* Download the data for more **sequence information**, including 20nt-sgRNAs and 4nt-PAMs.
    - Some concepts:
        - Base editing efficiency (i.e. proportion of all edited outcomes in editing window to all outcomes).
        - Base editing outcome proportion(**pred_prop**) shows the predicted proportion of edited substrate(i.e. proportion of edited sixth base to all edited outcomes in editing window).
        - Absolute frequency(**pred_freq**) shows the absolute predicted frequency of edited substrate(i.e. proportion of edited sixth base to all outcomes).
        - The threshold of Editable sites were defined as **outcome proportion >= 90% & efficiency>=5%**.
        
    """
)
