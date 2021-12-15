from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Dict, List, Tuple, Union

from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import FeatureSet, Track
from Bio.SeqFeature import FeatureLocation, SeqFeature
from reportlab.lib import colors
from reportlab.lib.units import cm

from gbkviz.align_coord import AlignCoord
from gbkviz.genbank import Genbank


@dataclass
class DrawGenbankFig:
    """Draw Genbank Figure Class"""

    gbk_list: List[Genbank]
    min_ranges: List[int]
    max_ranges: List[int]
    align_coords: List[AlignCoord]
    show_label: bool
    show_scale: bool
    show_ticks: bool
    label_type: str
    feature_symbol: str
    label_angle: int
    scaleticks_interval: int
    label_fsize: int
    scaleticks_fsize: int
    fig_width: int
    fig_track_height: int
    fig_track_size: float
    target_feature_types: List[str]
    feature2color: Dict[str, str]
    cross_link_color: str
    inverted_cross_link_color: str

    def __post_init__(self):
        self.gd = self._setup_genome_diagram()

    @property
    def max_range_length(self) -> int:
        """Max range length"""
        return max(
            [max_r - min_r for min_r, max_r in zip(self.min_ranges, self.max_ranges)]
        )

    @property
    def draw_pagesize(self) -> Tuple[float, float]:
        """Draw width * height pagesize (cm)"""
        width = self.fig_width * cm
        height = self.fig_track_height * len(self.gbk_list) * cm
        return (width, height)

    def get_figure(self, format: str) -> Union[str, bytes]:
        """Get genome diagram figure

        Args:
            format (str): Figure format ('jpg'|'png'|'pdf'|'svg')

        Returns:
            Union[str, bytes]: Figure string or bytes
        """
        format = format.lower()
        if format == "svg":
            handle = StringIO()
            self.gd.write(handle, format)
            return handle.getvalue()
        else:
            return self.gd.write_to_string(format)

    def write_figure(self, outfile: Union[str, Path]) -> None:
        """Write genome diagram figure

        Args:
            outfile (Union[str, Path]): Output file path
        """
        format = Path(outfile).suffix.replace(".", "")
        self.gd.write(outfile, format)

    def _setup_genome_diagram(self) -> GenomeDiagram.Diagram:
        # Create GenomeDiagram.Diagram object
        gd = GenomeDiagram.Diagram("Genbank Genome Diagram")

        for gbk, min_range, max_range in zip(
            self.gbk_list, self.min_ranges, self.max_ranges
        ):
            # Add track of one genbank
            gd_feature_set: FeatureSet = gd.new_track(
                track_level=0,
                name=gbk.name,
                greytrack=False,  # Disable greytrack
                greytrack_labels=0,
                greytrack_fontcolor=colors.black,
                start=0,
                end=max_range - min_range,
                scale=self.show_scale,
                scale_fontsize=self.scaleticks_fsize,
                scale_fontangle=0,
                scale_format="SInt",
                scale_color=colors.black,
                scale_ticks=self.show_ticks,
                scale_smallticks=0.4,
                scale_smalltick_interval=self.scaleticks_interval,
                scale_smalltick_labels=True,
                scale_largeticks=0.4,
                scale_largetick_interval=9999999999,  # Set large value to disable
                scale_largetick_labels=False,
                axis_labels=True,
            ).new_set()

            for feature in gbk.extract_features(self.target_feature_types):
                # Make location fixed feature
                feature = SeqFeature(
                    location=FeatureLocation(
                        feature.location.start - min_range,
                        feature.location.end - min_range,
                        feature.strand,
                    ),
                    type=feature.type,
                    qualifiers=feature.qualifiers,
                )

                # Get draw feature 'label_name', 'feature_color', 'label_angle'
                if self.label_type in ("gene", "protein_id", "locus_tag", "product"):
                    label_name = feature.qualifiers.get(self.label_type, [""])[0]
                    color = self.feature2color[feature.type]
                else:
                    raise ValueError(
                        f"Invalid label type `{self.label_type}` detected!!"
                    )
                strand_label_angle = (
                    180 - self.label_angle if feature.strand == -1 else self.label_angle
                )

                # Add feature to genbank track
                gd_feature_set.add_feature(
                    feature=feature,
                    color=color,
                    name=label_name,
                    label=self.show_label,
                    label_size=self.label_fsize,
                    label_angle=strand_label_angle,
                    label_position="middle",  # "start", "middle", "end"
                    sigil=self.feature_symbol,  # "BOX", "ARROW", "OCTO", "BIGARROW"
                    arrowhead_length=0.5,  # Default: 0.5
                    arrowshaft_height=0.3,
                )

        # Get cross links
        name2track: Dict[str, Track] = {track.name: track for track in gd.get_tracks()}
        name2start: Dict[str, int] = {
            gbk.name: mr for gbk, mr in zip(self.gbk_list, self.min_ranges)
        }
        cross_links = []
        for align_coord in self.align_coords:
            cross_link = align_coord.get_cross_link(
                name2track=name2track,
                name2start=name2start,
                normal_color=self.cross_link_color,
                inverted_color=self.inverted_cross_link_color,
            )
            cross_links.append(cross_link)

        # Set figure draw settings
        gd.draw(
            format="linear",
            orientation="landscape",
            pagesize=self.draw_pagesize,
            fragments=1,
            start=0,
            end=self.max_range_length,
            tracklines=False,
            track_size=self.fig_track_size,
            cross_track_links=cross_links,
        )
        return gd
