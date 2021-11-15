import os
from pathlib import Path
from typing import List

import streamlit as st
from streamlit.delta_generator import DeltaGenerator

from gbkviz.align_coord import AlignCoord
from gbkviz.draw_param import DrawParam
from gbkviz.genbank import Genbank
from gbkviz.genbank_diagram import gbk2fig
from gbkviz.genome_align import GenomeAlign

# Page basic configuration
st.set_page_config(
    page_title="GBKviz: Genbank Data Visualization Tool",
    layout="centered",
    initial_sidebar_state="expanded",
    menu_items={
        "Get Help": "https://github.com/moshi4/GBKviz",
        "Report a Bug": "https://github.com/moshi4/GBKviz",
    },
)

st.header("GBKviz: Genbank Data Visualization Tool")

###########################################################
# Sidebar Parameters Widgets
###########################################################

# GitHub repository hyperlink
repo_hyperlink = "[GBKviz (GitHub)](https://github.com/moshi4/GBKviz/)"
st.sidebar.markdown(repo_hyperlink)

# Genbank files upload widgets
upload_files = st.sidebar.file_uploader(
    label="Upload your genbank files (*.gb|*.gbk)",
    type=["gb", "gbk"],
    accept_multiple_files=True,
)

if upload_files:
    # Download format selectbox widgets
    fig_format = st.sidebar.selectbox(
        label="Download Figure Format (JPG/PNG/SVG/PDF)",
        options=["JPG", "PNG", "SVG", "PDF"],
        index=0,
    )
    fig_format = str(fig_format).lower()

    # Visibility control checkbox widgets
    check_cols = List[DeltaGenerator]
    check_cols = st.sidebar.columns(3)
    show_label = check_cols[0].checkbox("Label", False)
    show_scale = check_cols[1].checkbox("Scale", True)
    show_ticks = check_cols[2].checkbox("ScaleTicks", True)

    # Figure appearence control widgets
    fig_appearence_cols: List[DeltaGenerator]
    fig_appearence_cols = st.sidebar.columns(2)
    label_type = fig_appearence_cols[0].selectbox(
        label="Feature Label Type",
        options=["gene", "protein_id", "locus_tag", "product"],
        index=0,
    )
    feature_symbol = fig_appearence_cols[1].selectbox(
        label="Feature Symbol",
        options=["BIGARROW", "BOX", "ARROW", "OCTO"],
        index=0,
    )

    label_angle_str = fig_appearence_cols[0].selectbox(
        label="Label Angle",
        options=["0", "15", "30", "45", "60", "75", "90"],
        index=2,
    )
    label_angle = int(label_angle_str)

    sint2int = {
        "1Kbp": 1000,
        "5Kbp": 5000,
        "10Kbp": 10000,
        "50Kbp": 50000,
        "100Kbp": 100000,
        "1Mbp": 1000000,
    }
    scaleticks_sint = fig_appearence_cols[1].selectbox(
        label="ScaleTicks Interval",
        options=list(sint2int.keys()),
        index=2,
    )
    scaleticks_interval = sint2int[scaleticks_sint]

    label_fsize = fig_appearence_cols[0].number_input(
        label="Label Font Size",
        min_value=0,
        max_value=100,
        value=10,
        step=1,
    )
    scaleticks_fsize = fig_appearence_cols[1].number_input(
        label="ScaleTicks Font Size",
        min_value=0,
        max_value=100,
        value=8,
        step=1,
    )

    # Figure "Width", "TrackHeight", "TrackSizeRatio" control slider widgets
    slider_cols = List[DeltaGenerator]
    slider_cols = st.sidebar.columns(2)
    fig_width = slider_cols[0].slider(
        label="Fig Width(cm)",
        min_value=10,
        max_value=100,
        value=25,
        step=5,
    )
    fig_track_height = slider_cols[1].slider(
        label="Fig Track Height(cm)",
        min_value=1,
        max_value=10,
        value=3,
        step=1,
    )
    fig_track_size = slider_cols[0].slider(
        label="Fig Track Size",
        min_value=0.1,
        max_value=1.0,
        value=0.5,
        step=0.1,
    )

    # Target features selection widget
    target_features = st.sidebar.multiselect(
        label="Target Feature Selection",
        options=["CDS", "gene", "tRNA", "misc_feature"],
        default=["CDS"],
    )

    # Features colorpicker widgets
    color_cols: List[DeltaGenerator]
    color_cols = st.sidebar.columns(4)
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

    genome_comparison = None
    if len(upload_files) >= 2 and GenomeAlign.check_requirements():
        # Genome comparison type selectbox widget
        genome_comparison = st.sidebar.selectbox(
            label="Comparison Type",
            options=[
                None,
                "Nucleotide One-to-One",
                "Nucleotide Many-to-Many",
                "Protein One-to-One",
                "Protein Many-to-Many",
            ],
            index=0,
        )

    ###########################################################
    # Main Screen Widgets
    ###########################################################
    gbk_list: List[Genbank] = []
    min_value_list: List[int] = []
    max_value_list: List[int] = []
    gbk_info_list: List[str] = []

    gbk_info_placeholder = st.empty()
    download_btn_placeholder = st.empty()
    fig_placeholder = st.empty()

    with st.form(key="form"):

        st.form_submit_button(label="Update Figure")
        input_cols: List[DeltaGenerator]
        input_cols = st.columns(2)

        for upload_gbk_file in upload_files:
            gbk = Genbank.read_upload_gbk_file(upload_gbk_file)
            gbk_list.append(gbk)
            max_length = gbk.max_length

            # Min-Max range input widget
            range_label = f"{gbk.name} Min-Max Range (Max={max_length:,} bp)"
            default_max_value = max_length if max_length <= 50000 else 50000
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
                value=default_max_value,
                step=1000,
            )
            length = int(max_value - min_value)
            gbk_info_list.append(
                f"{gbk.name} ({min_value:,} - {max_value:,} bp), Length={length:,} bp"
            )
            min_value_list.append(int(min_value))
            max_value_list.append(int(max_value))

    # Genome comparison
    align_coords_list: List[AlignCoord] = []
    if genome_comparison is not None:
        genome_fasta_files: List[Path] = []
        gbkviz_tmpdir = Path("tmp_gbkviz")
        os.makedirs(gbkviz_tmpdir, exist_ok=True)
        # with tempfile.TemporaryDirectory() as gbkviz_tmpdir:
        for gbk in gbk_list:
            # Make genome fasta file
            genome_fasta_file = Path(gbkviz_tmpdir) / (gbk.name + ".fa")
            if not genome_fasta_file.exists():
                gbk.write_genome_fasta(genome_fasta_file)
            genome_fasta_files.append(genome_fasta_file)
        seqtype, maptype = genome_comparison.split(" ")
        genome_align = GenomeAlign(genome_fasta_files, seqtype, maptype)
        align_coords_list = genome_align.run()
        print(gbkviz_tmpdir)
        print(align_coords_list)

    # Uploaded genbank file information
    all_gbk_info = ""
    for cnt, gbk_info in enumerate(gbk_info_list, 1):
        all_gbk_info += f"Track{cnt:02d}: {gbk_info}  \n"
    gbk_info_placeholder.markdown(all_gbk_info)

    # Display genbank visualization figure
    jpg_bytes, format_bytes = gbk2fig(
        gbk_list=gbk_list,
        start_pos_list=min_value_list,
        end_pos_list=max_value_list,
        feature2color=feature2color,
        fig_width=fig_width,
        fig_track_height=fig_track_height,
        show_label=show_label,
        show_ticks=show_ticks,
        show_scale=show_scale,
        label_type=label_type,
        feature_symbol=feature_symbol,
        label_angle=label_angle,
        target_features=target_features,
        scaleticks_interval=scaleticks_interval,
        fig_track_size=fig_track_size,
        label_fsize=int(label_fsize),
        scaleticks_fsize=int(scaleticks_fsize),
        fig_format=fig_format,
    )
    fig_placeholder.image(jpg_bytes, use_column_width="never")

    download_btn_placeholder.download_button(
        label=f"Download Figure ({fig_format.upper()} Format)",
        data=format_bytes,
        file_name=f"gbkviz_figure.{fig_format}",
    )
