from io import StringIO
from typing import Dict, List, Tuple, Union

from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import FeatureSet, Track
from Bio.SeqFeature import FeatureLocation, SeqFeature
from reportlab.lib import colors
from reportlab.lib.units import cm

from gbkviz.align_coord import AlignCoord
from gbkviz.genbank import Genbank


def draw_gbk_fig(
    gbk_list: List[Genbank],
    min_ranges: List[int],
    max_ranges: List[int],
    align_coords: List[AlignCoord],
    fig_format: str,
    show_label: bool,
    show_scale: bool,
    show_ticks: bool,
    label_type: str,
    feature_symbol: str,
    label_angle: int,
    scaleticks_interval: int,
    label_fsize: int,
    scaleticks_fsize: int,
    fig_width: int,
    fig_track_height: int,
    fig_track_size: float,
    target_features: List[str],
    feature2color: Dict[str, str],
    cross_link_color: str,
    inverted_cross_link_color: str,
) -> Tuple[bytes, Union[str, bytes]]:
    """Draw genbank features figure

    Args:
        gbk_list (List[Genbank]): Genbank object list
        min_ranges (List[int]): Min range list
        max_ranges (List[int]): Max range list
        align_coords (List[AlignCoord]): Alignment coordinate list
        fig_format (str): Figure format (JPG/PNG/SVG/PDF)
        show_label (bool): Show label or not
        show_scale (bool): Show scale or not
        show_ticks (bool): Show scaleticks or not
        label_type (str): Label type ('gene', 'protein_id', 'locus_tag', 'product')
        feature_symbol (str): Feature symbol ('BIGARROW', 'BOX', 'ARROW', 'OCTO')
        label_angle (int): Label angle
        scaleticks_interval (int): Scaleticks interval
        label_fsize (int): Label font size
        scaleticks_fsize (int): Scaleticks font size
        fig_width (int): Figure width (cm)
        fig_track_height (int): Figure one track height (cm)
        fig_track_size (float): Figure track size ratio (0.0 - 1.0)
        target_features (List[str]): Target features \
                                    ('CDS', 'gene', 'tRNA', 'misc_feature')
        feature2color (Dict[str, str]): Feature color dict (e.g. 'CDS' -> '#FFA500')
        cross_link_color (str): Cross link hexcolor
        inverted_cross_link_color (str): Inverted cross link hexcolor

    Returns:
        Tuple[bytes, Union[str, bytes]]: JPG image, Specified format image

    Note:
        SVG format image type = string
        Other format image type = bytes
    """
    # Create GenomeDiagram object
    gd = GenomeDiagram.Diagram("Genbank Genome Diagram")

    for gbk, min_range, max_range in zip(gbk_list, min_ranges, max_ranges):

        # Add track of one genbank
        gd_feature_set: FeatureSet = gd.new_track(
            track_level=0,
            name=gbk.name,
            greytrack=False,  # Disable greytrack
            greytrack_labels=0,
            greytrack_fontcolor=colors.black,
            start=0,
            end=max_range - min_range,
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
            # Make location fixed feature
            start = feature.location.start - min_range
            end = feature.location.end - min_range
            strand = feature.strand
            feature = SeqFeature(
                location=FeatureLocation(start, end, strand),
                type=feature.type,
                qualifiers=feature.qualifiers,
            )

            # Get draw feature 'label_name', 'feature_color', 'label_angle'
            if label_type in ("gene", "protein_id", "locus_tag", "product"):
                label_name = feature.qualifiers.get(label_type, [""])[0]
                color = feature2color[feature.type]
            else:
                raise ValueError(f"Invalid label type `{label_type}` detected!!")
            strand_label_angle = (
                label_angle if feature.strand == 1 else 180 - label_angle
            )

            # Add feature to genbank track
            gd_feature_set.add_feature(
                feature=feature,
                color=color,
                name=label_name,
                label=show_label,
                label_size=label_fsize,
                label_angle=strand_label_angle,
                label_position="middle",  # "start", "middle", "end"
                sigil=feature_symbol,  # "BOX", "ARROW", "OCTO", "JAGGY", "BIGARROW"
                arrowhead_length=0.5,  # Default: 0.5
                arrowshaft_height=0.3,
            )

    # Get cross links
    name2track: Dict[str, Track] = {}
    for track in gd.get_tracks():
        name2track[track.name] = track
    name2start: Dict[str, int] = {}
    for gbk, min_range in zip(gbk_list, min_ranges):
        name2start[gbk.name] = min_range
    cross_links = []
    for align_coord in align_coords:
        cross_link = align_coord.get_cross_link(
            name2track=name2track,
            name2start=name2start,
            normal_color=cross_link_color,
            inverted_color=inverted_cross_link_color,
        )
        cross_links.append(cross_link)

    # Set figure draw settings
    gd.draw(
        format="linear",
        orientation="landscape",
        pagesize=(fig_width * cm, fig_track_height * len(gbk_list) * cm),  # X, Y
        fragments=1,
        start=0,
        end=max([max_r - min_r for min_r, max_r in zip(min_ranges, max_ranges)]),
        tracklines=False,
        track_size=fig_track_size,
        cross_track_links=cross_links,
    )

    # Make 1. JPG image, 2. Specified format image
    jpg_bytes = gd.write_to_string(output="jpg")
    if fig_format == "svg":
        handle = StringIO()
        gd.write(handle, fig_format)
        return jpg_bytes, handle.getvalue()
    else:
        format_bytes = gd.write_to_string(output=fig_format)
        return jpg_bytes, format_bytes
