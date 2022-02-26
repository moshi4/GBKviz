from pathlib import Path
from typing import List

import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from streamlit.uploaded_file_manager import UploadedFile

from gbkviz.align_coord import AlignCoord
from gbkviz.draw_genbank_fig import DrawGenbankFig
from gbkviz.genbank import Genbank
from gbkviz.genome_align import GenomeAlign
from gbkviz.util import load_files, make_session_dir, remove_olddir

# Page basic configuration
st.set_page_config(
    page_title="GBKviz: Genbank Data Visualization WebApp",
    layout="centered",
    initial_sidebar_state="expanded",
    menu_items={
        "Get Help": "https://github.com/moshi4/GBKviz",
        "Report a Bug": "https://github.com/moshi4/GBKviz",
    },
)

st.header("GBKviz: Genbank Data Visualization WebApp")

###########################################################
# Sidebar Parameters Widgets
###########################################################

# GitHub repository hyperlink
repo_hyperlink = "[GBKviz (GitHub)](https://github.com/moshi4/GBKviz/)"
st.sidebar.markdown(repo_hyperlink)

# Load example files or Upload files
if st.sidebar.checkbox(label="Load example genbank files", value=False):
    genbank_dir = Path(__file__).parent / "genbank"
    gbk_files = list(genbank_dir.glob("*.gbk"))
    upload_files = load_files(gbk_files)
else:
    with st.sidebar.expander(label="Toggle Genbank Upload Box", expanded=True):
        # Genbank files upload widgets
        upload_files: List[UploadedFile]
        upload_files = st.file_uploader(
            label="Upload your genbank files (*.gb|*.gbk)",
            type=["gb", "gbk"],
            accept_multiple_files=True,
        )

if upload_files:
    # Visibility control checkbox widgets
    check_cols = List[DeltaGenerator]
    check_cols = st.sidebar.columns(3)
    show_label = check_cols[0].checkbox("Label", False)
    show_scale = check_cols[1].checkbox("Scale", True)
    show_ticks = check_cols[2].checkbox("ScaleTicks", False)

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

    # Target feature types widget
    target_feature_types = st.sidebar.multiselect(
        label="Target Feature Types",
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
    cross_link_color = ""
    inverted_cross_link_color = ""
    if len(upload_files) >= 2 and GenomeAlign.check_requirements():
        # Genome comparison type selectbox widget
        genome_comparison = st.sidebar.selectbox(
            label="Genome Comparison Type",
            options=[
                None,
                "Nucleotide One-to-One",
                "Nucleotide Many-to-Many",
                "Protein One-to-One",
                "Protein Many-to-Many",
            ],
            index=0,
        )

        # Genome comparison colorpicker widgets
        cross_link_color_cols: List[DeltaGenerator]
        cross_link_color_cols = st.sidebar.columns(2)
        cross_link_color = cross_link_color_cols[0].color_picker(
            label="Cross Link (Normal)", value="#0000FF"
        )
        inverted_cross_link_color = cross_link_color_cols[1].color_picker(
            label="Cross Link (Inverted)", value="#FF0000"
        )

    ###########################################################
    # Main Screen Widgets
    ###########################################################
    gbk_list: List[Genbank] = []
    min_ranges: List[int] = []
    max_ranges: List[int] = []
    gbk_info_list: List[str] = []

    gbk_info_placeholder = st.empty()
    dl_btn_cols = st.columns([3, 3, 5])
    dl_png_btn_placeholder = dl_btn_cols[0].empty()
    dl_svg_btn_placeholder = dl_btn_cols[1].empty()
    dl_align_coords_btn_placeholder = dl_btn_cols[2].empty()
    fig_placeholder = st.empty()

    with st.form(key="form"):

        st.form_submit_button(label="Update Figure")

        st.markdown("**Display Genome Min-Max Range & Reverse Option**")

        range_cols: List[DeltaGenerator]
        range_cols = st.columns([3, 3, 1])

        for upload_gbk_file in upload_files:
            gbk = Genbank.read_upload_gbk_file(upload_gbk_file)
            gbk_list.append(gbk)

            # Min-Max range input widget
            range_label = f"{gbk.name} (Max={gbk.max_length:,} bp)"
            min_range = range_cols[0].number_input(
                label=range_label,
                min_value=1,
                max_value=gbk.max_length,
                value=1,
                step=1000,
                key=gbk.name,
            )
            default_max_range = gbk.max_length if gbk.max_length <= 100000 else 100000
            max_range = range_cols[1].number_input(
                label="",
                min_value=1,
                max_value=gbk.max_length,
                value=default_max_range,
                step=1000,
                key=gbk.name,
            )
            reverse = range_cols[2].selectbox(
                label="Reverse",
                options=["Yes", "No"],
                index=1,
                key=gbk.name,
            )
            gbk.reverse = True if reverse == "Yes" else False

            gbk_info = (
                f"{gbk.name} ({min_range:,} - {max_range:,} bp), "
                + f"Length={int(max_range - min_range + 1):,} bp, "
                + f"CDS={gbk.count_feature(int(min_range), int(max_range)):,}"
            )
            gbk_info_list.append(gbk_info)
            min_ranges.append(int(min_range))
            max_ranges.append(int(max_range))

    # Genome alignment
    align_coords: List[AlignCoord] = []
    gbkviz_tmpdir = Path.home() / ".gbkviz"
    gbkviz_tmpdir.mkdir(exist_ok=True)
    if genome_comparison is not None:
        genome_fasta_files: List[Path] = []
        gbkviz_session_tmpdir = make_session_dir(gbkviz_tmpdir)
        for gbk in gbk_list:
            # Make genome fasta file
            suffix = "_reverse.fa" if gbk.reverse else ".fa"
            genome_fasta_file = Path(gbkviz_session_tmpdir) / (gbk.name + suffix)
            if not genome_fasta_file.exists():
                gbk.write_genome_fasta(genome_fasta_file)
            genome_fasta_files.append(genome_fasta_file)
        # Run MUMmer genome alignment
        seqtype, maptype = genome_comparison.split(" ")
        genome_align = GenomeAlign(
            genome_fasta_files, gbkviz_session_tmpdir, seqtype, maptype
        )
        align_coords = genome_align.run()

    # Remove old genome comparison result directory
    for session_dir in gbkviz_tmpdir.iterdir():
        remove_olddir(session_dir)

    # Uploaded genbank file information
    all_gbk_info = ""
    for cnt, gbk_info in enumerate(gbk_info_list, 1):
        all_gbk_info += f"Track{cnt:02d}: {gbk_info}  \n"
    gbk_info_placeholder.markdown(all_gbk_info)

    # Create visualization and comparison figure
    dgf = DrawGenbankFig(
        gbk_list=gbk_list,
        min_ranges=min_ranges,
        max_ranges=max_ranges,
        align_coords=align_coords,
        show_label=show_label,
        show_scale=show_scale,
        show_ticks=show_ticks,
        label_type=label_type,
        feature_symbol=feature_symbol,
        label_angle=label_angle,
        scaleticks_interval=scaleticks_interval,
        label_fsize=int(label_fsize),
        scaleticks_fsize=int(scaleticks_fsize),
        fig_width=fig_width,
        fig_track_height=fig_track_height,
        fig_track_size=fig_track_size,
        cross_link_color=cross_link_color,
        inverted_cross_link_color=inverted_cross_link_color,
        target_feature_types=target_feature_types,
        feature2color=feature2color,
    )

    # Show figure
    png_bytes = dgf.get_figure("png")
    fig_placeholder.image(png_bytes, use_column_width="never")

    # Download figure button widget
    dl_png_btn_placeholder.download_button(
        label="Download PNG Figure",
        data=png_bytes,
        file_name="gbkviz_figure.png",
    )
    dl_svg_btn_placeholder.download_button(
        label="Download SVG Figure",
        data=dgf.get_figure("svg"),
        file_name="gbkviz_figure.svg",
    )

    # Download align coords button widget
    if align_coords:
        header = (
            "REF_START\tREF_END\tQUERY_START\tQUERY_END\tREF_LENGTH\t"
            + "QUERY_LENGTH\tIDENTITY\tREF_NAME\tQUERY_NAME\n"
        )
        dl_align_coords_btn_placeholder.download_button(
            label="Download Comparison Result",
            data=header + "\n".join([ac.as_tsv_format for ac in align_coords]),
            file_name="gbkviz_comparison.tsv",
        )
else:
    # No Uploaded files, display toppage contents
    tutorial_md_file = Path(__file__).parent / "toppage.md"
    with open(tutorial_md_file) as f:
        tutorial_content = f.read()
    st.markdown(tutorial_content)
