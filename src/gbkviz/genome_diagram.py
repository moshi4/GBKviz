from typing import Dict, List

from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import FeatureSet
from reportlab.lib import colors
from reportlab.lib.units import cm

from gbkviz.genbank import Genbank


def cds_feature_list2fig(
    gbk_list: List[Genbank],
    start_pos_list: List[int],
    end_pos_list: List[int],
    feature2color: Dict[str, str],
    fig_width: int,
    fig_track_height: int,
    show_label: bool,
    show_ticks: bool,
    show_scale: bool,
    label_type: str,
    symbol: str,
    label_angle: int,
    target_features: List[str],
) -> bytes:

    gd = GenomeDiagram.Diagram("Genbank Genome Diagram")

    print("####################")
    for track_cnt, (gbk, start_pos, end_pos) in enumerate(
        zip(gbk_list, start_pos_list, end_pos_list), 1
    ):
        print(gbk.name, start_pos, end_pos)
        track_cnt = len(gbk_list) - 1

        gd_feature_set: FeatureSet = gd.new_track(
            track_level=track_cnt,
            name=gbk.name,
            greytrack=False,
            greytrack_labels=0,
            greytrack_fontcolor=colors.black,
            start=start_pos,
            end=end_pos,
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

        for cds_feature in gbk.extract_features(target_features):
            qualifiers = cds_feature.qualifiers
            if label_type in ("gene", "protein_id", "locus_tag", "product"):
                label_name = qualifiers.get(label_type, [""])[0]
                color = feature2color[cds_feature.type]
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

    one_track_size = (fig_width * cm, fig_track_height * len(gbk_list) * cm)

    gd.draw(
        format="linear",
        orientation="landscape",
        pagesize=one_track_size,  # X, Y
        fragments=1,
        start=min(start_pos_list),
        end=max(end_pos_list),
        tracklines=False,
        # track_size=0.05,
        track_size=0.3,
    )

    return gd.write_to_string(output="jpg")
