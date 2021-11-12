from pathlib import Path
from typing import List

import streamlit as st
from streamlit.delta_generator import DeltaGenerator

from gbkviz.genbank import read_upload_gbk_file
from gbkviz.genome_diagram import cds_feature_list2fig

st.header("GBKviz: Genbank Data Visualization Tool")

###########################################################
# Sidebar Parameters Widgets
###########################################################
with st.sidebar:

    # Genbank files upload widgets
    upload_file_list = st.file_uploader(
        label="Upload your genbank files (*.gb|*.gbk)",
        type=["gb", "gbk"],
        accept_multiple_files=True,
    )

    # Features select option
    target_features = st.multiselect(
        label="Features select option",
        options=["CDS", "gene", "tRNA", "misc_feature"],
        default=["CDS"],
    )

    # Colorpicker widgets
    color_cols: List[DeltaGenerator]
    color_cols = st.columns(4)
    cds_color = color_cols[0].color_picker(label="CDS", value="#FFA500")
    gene_color = color_cols[1].color_picker(label="gene", value="#0FE8E4")
    trna_color = color_cols[2].color_picker(label="tRNA", value="#E80F0F")
    misc_color = color_cols[3].color_picker(label="misc", value="#E80FC6")
    feature2color = {
        "CDS": cds_color,
        "gene": gene_color,
        "tRNA": trna_color,
        "misc_feature": misc_color,
    }

    #
    label_type = st.selectbox(
        label="Features Label",
        options=["gene", "protein_id", "locus_tag", "product"],
        index=0,
    )
    symbol = st.selectbox(
        label="Features Symbol",
        options=["BIGARROW", "BOX", "ARROW", "OCTO", "JAGGY"],
        index=0,
    )

    check_cols = List[DeltaGenerator]
    check_cols = st.columns(3)
    show_label = check_cols[0].checkbox("Label", True)
    show_ticks = check_cols[1].checkbox("Ticks", True)
    show_scale = check_cols[2].checkbox("Scale", True)

    slider_cols = List[DeltaGenerator]
    slider_cols = st.columns(2)
    fig_width = slider_cols[0].slider(
        label="Fig Width(cm)",
        min_value=10,
        max_value=100,
        value=30,
        step=5,
    )
    fig_track_height = slider_cols[1].slider(
        label="Fig Track Height(cm)",
        min_value=5,
        max_value=100,
        value=5,
        step=5,
    )
    label_angle = slider_cols[0].slider(
        label="Label Angle",
        min_value=0,
        max_value=90,
        value=30,
        step=15,
    )


###########################################################
# Main Screen Widgets
###########################################################
if upload_file_list:

    gbk_info_list = []
    gbk_info_placeholder = st.empty()
    image_placeholder = st.empty()

    with st.form(key="form"):

        st.form_submit_button(label="Update Figure")
        input_cols: List[DeltaGenerator]
        input_cols = st.columns(2)

        for upload_gbk_file in upload_file_list:
            gbk = read_upload_gbk_file(upload_gbk_file)
            cds_feature_list = gbk.extract_features(target_features)
            max_length = gbk.max_length
            cds_count = len(cds_feature_list)

            gbk_name = Path(gbk.name).stem

            # Min-Max range input widget
            range_label = f"{gbk_name} Min-Max Range (Max={max_length:,} bp)"
            min_value = input_cols[0].number_input(
                label=range_label,
                min_value=0,
                max_value=max_length,
                value=0,
                step=1000,
            )
            max_value = input_cols[1].number_input(
                label="",
                min_value=0,
                max_value=max_length,
                value=10000,
                step=1000,
            )
            gbk_info_list.append(f"{gbk_name} ({min_value:,} - {max_value:,} bp)")

    all_gbk_info = ""
    for cnt, gbk_info in enumerate(gbk_info_list, 1):
        all_gbk_info += f"Track{cnt:02d}: {gbk_info}  \n"
    gbk_info_placeholder.markdown(all_gbk_info)

    fig_bytes = cds_feature_list2fig(
        cds_feature_list=cds_feature_list,
        start_pos=int(min_value),
        end_pos=int(max_value),
        feature2color=feature2color,
        fig_width=fig_width,
        fig_track_height=fig_track_height,
        show_label=show_label,
        show_ticks=show_ticks,
        show_scale=show_scale,
        label_type=label_type,
        symbol=symbol,
        label_angle=label_angle,
    )

    image_placeholder.image(fig_bytes, use_column_width="never")
