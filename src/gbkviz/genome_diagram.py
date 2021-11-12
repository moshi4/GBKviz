from typing import List

from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import FeatureSet
from Bio.SeqFeature import SeqFeature
from reportlab.lib import colors
from reportlab.lib.units import cm


def cds_feature_list2fig(
    cds_feature_list: List[SeqFeature],
    start_pos: int,
    end_pos: int,
    color: str,
    fig_width: int,
    fig_track_height: int,
    show_label: bool,
    show_ticks: bool,
    show_scale: bool,
    label_type: str,
    symbol: str,
    label_angle: int,
) -> bytes:

    gd = GenomeDiagram.Diagram("Genbank Genome Diagram")
    gd_feature_set: FeatureSet = gd.new_track(
        1,
        # name=gbk_file.stem,
        name="test",
        greytrack=False,
        greytrack_labels=0,
        greytrack_fontcolor=colors.black,
        # start=10000,
        # end=30000,
        # Scale & Scale Ticks property
        scale=show_scale,
        scale_fontsize=8,
        scale_fontangle=0,
        scale_format="SInt",
        scale_color=colors.black,
        scale_ticks=show_ticks,
        scale_smallticks=0.4,
        scale_smalltick_interval=5000,
        scale_smalltick_labels=True,
        # scale_largeticks=0.6,
        # scale_largetick_interval=10000,
        # scale_largetick_labels=False,
        axis_label=True,
    ).new_set()

    for cds_feature in cds_feature_list:
        qualifiers = cds_feature.qualifiers
        if label_type in ("gene", "protein_id", "locus_tag", "product"):
            label_name = qualifiers.get(label_type, [""])[0]
        else:
            raise ValueError(f"Invalid label type `{label_type}` detected!!")

        strand_label_angle = (
            label_angle if cds_feature.strand == 1 else 180 - label_angle
        )

        gd_feature_set.add_feature(
            cds_feature,
            color=color,
            name=label_name,
            label=show_label,
            label_size=10,
            label_angle=strand_label_angle,
            label_position="middle",  # "start", "middle", "end"
            sigil=symbol,  # "BOX", "ARROW", "OCTO", "JAGGY", "BIGARROW"
            arrowhead_length=0.5,  # Default: 0.5
            arrowshaft_height=0.3,
        )

    one_track_size = (fig_width * cm, fig_track_height * cm)

    gd.draw(
        format="linear",
        orientation="landscape",
        pagesize=one_track_size,  # X, Y
        fragments=1,
        start=start_pos,
        end=end_pos,
        tracklines=False,
        # track_size=0.05,
        track_size=0.3,
    )

    return gd.write_to_string(output="jpg")

    # plot_file = Path("test.jpg")
    # file_ext = plot_file.suffix.replace(".", "")
    # gd.write(plot_file, file_ext)
