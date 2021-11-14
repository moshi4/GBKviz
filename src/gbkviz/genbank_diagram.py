from io import StringIO
from typing import Dict, List, Tuple, Union

from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import FeatureSet
from Bio.SeqFeature import FeatureLocation, SeqFeature
from reportlab.lib import colors
from reportlab.lib.units import cm

from gbkviz.genbank import Genbank


def gbk2fig(
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
    scaleticks_interval: int,
    fig_track_size: float,
    label_fsize: int,
    scaleticks_fsize: int,
    fig_format: str,
) -> Tuple[bytes, Union[str, bytes]]:

    gd = GenomeDiagram.Diagram("Genbank Genome Diagram")

    for gbk, start_pos, end_pos in zip(gbk_list, start_pos_list, end_pos_list):

        gd_feature_set: FeatureSet = gd.new_track(
            track_level=0,
            name=gbk.name,
            greytrack=False,  # Disable greytrack
            greytrack_labels=0,
            greytrack_fontcolor=colors.black,
            start=0,
            end=end_pos - start_pos,
            scale=show_scale,
            scale_fontsize=scaleticks_fsize,
            scale_fontangle=0,
            scale_format="SInt",
            scale_color=colors.black,
            scale_ticks=show_ticks,
            scale_smallticks=0.4,
            scale_smalltick_interval=scaleticks_interval,
            scale_smalltick_labels=True,
            scale_largeticks=0.4,
            scale_largetick_interval=9999999999,  # Set large value to disable largetick
            scale_largetick_labels=False,
            axis_label=True,
        ).new_set()

        for feature in gbk.extract_features(target_features):
            loc = feature.location
            start = loc.start - start_pos
            end = loc.end - start_pos
            strand = feature.strand
            feature = SeqFeature(
                location=FeatureLocation(start, end, strand),
                type=feature.type,
                qualifiers=feature.qualifiers,
            )

            qualifiers = feature.qualifiers
            if label_type in ("gene", "protein_id", "locus_tag", "product"):
                label_name = qualifiers.get(label_type, [""])[0]
                color = feature2color[feature.type]
            else:
                raise ValueError(f"Invalid label type `{label_type}` detected!!")

            strand_label_angle = (
                label_angle if feature.strand == 1 else 180 - label_angle
            )

            gd_feature_set.add_feature(
                feature=feature,
                color=color,
                name=label_name,
                label=show_label,
                label_size=label_fsize,
                label_angle=strand_label_angle,
                label_position="middle",  # "start", "middle", "end"
                sigil=symbol,  # "BOX", "ARROW", "OCTO", "JAGGY", "BIGARROW"
                arrowhead_length=0.5,  # Default: 0.5
                arrowshaft_height=0.3,
            )

    size_list = [e - s for s, e in zip(start_pos_list, end_pos_list)]

    gd.draw(
        format="linear",
        orientation="landscape",
        pagesize=(fig_width * cm, fig_track_height * len(gbk_list) * cm),  # X, Y
        fragments=1,
        start=0,
        end=max(size_list),
        tracklines=False,
        track_size=fig_track_size,
    )

    jpg_bytes = gd.write_to_string(output="jpg")
    if fig_format == "svg":
        # format_bytes = gd.write_to_string(output="SVG")
        handle = StringIO()
        gd.write(handle, fig_format)
        return jpg_bytes, handle.getvalue()
    else:
        format_bytes = gd.write_to_string(output=fig_format)
        return jpg_bytes, format_bytes
