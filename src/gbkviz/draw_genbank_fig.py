from io import StringIO
from pathlib import Path
from typing import Dict, List, Tuple, Union

from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import FeatureSet
from Bio.SeqFeature import FeatureLocation, SeqFeature
from reportlab.lib import colors
from reportlab.lib.units import cm

from gbkviz.align_coord import AlignCoord
from gbkviz.genbank import Genbank


class DrawGenbankFig:
    """Draw Genbank Figure Class"""

    def __init__(
        self,
        gbk_list: List[Genbank],
        align_coords: List[AlignCoord] = [],
        show_label: bool = False,
        show_scale: bool = True,
        show_ticks: bool = False,
        label_type: str = "gene",
        feature_symbol: str = "BIGARROW",
        label_angle: int = 30,
        scaleticks_interval: int = 10000,
        label_fsize: int = 10,
        scaleticks_fsize: int = 8,
        fig_width: int = 25,
        fig_track_height: int = 3,
        fig_track_size: float = 0.5,
        cross_link_color: str = "#0000FF",
        inverted_cross_link_color: str = "#FF0000",
        target_feature_types: List[str] = ["CDS"],
        feature2color: Dict[str, str] = {
            "CDS": "#FFA500",
            "gene": "#0FE8E4",
            "tRNA": "#E80F0F",
            "misc_feature": "#E80FC6",
        },
        cds_limit_num: int = 500,
    ):
        """DrawGenbankFig constructor

        Args:
            gbk_list (List[Genbank]): Genbank class objects
            align_coords (List[AlignCoord], optional): AlignCoord class objects
            show_label (bool, optional): Show label or not
            show_scale (bool, optional): Show scale or not
            show_ticks (bool, optional): Show ticks or not
            label_type (str, optional): Label type
            feature_symbol (str, optional): Feature symbol
            label_angle (int, optional): Label angle
            scaleticks_interval (int, optional): Scaleticks interval
            label_fsize (int, optional): Label font size
            scaleticks_fsize (int, optional): Scaleticks font size
            fig_width (int, optional): Figure width (cm)
            fig_track_height (int, optional): Figure track height (cm)
            fig_track_size (float, optional): Figure track size
            cross_link_color (str, optional): Cross link color
            inverted_cross_link_color (str, optional): Inverted cross link color
            target_feature_types (List[str], optional): Target feature types
            feature2color (Dict[str, str], optional): Feature colors dictionary
            cds_limit_num (int, optional): CDS limit number to be drawn
        """
        self.gbk_list: List[Genbank] = gbk_list
        self.align_coords: List[AlignCoord] = align_coords
        self.show_label: bool = show_label
        self.show_scale: bool = show_scale
        self.show_ticks: bool = show_ticks
        self.label_type: str = label_type
        self.feature_symbol: str = feature_symbol
        self.label_angle: int = label_angle
        self.scaleticks_interval: int = scaleticks_interval
        self.label_fsize: int = label_fsize
        self.scaleticks_fsize: int = scaleticks_fsize
        self.fig_width: int = fig_width
        self.fig_track_height: int = fig_track_height
        self.fig_track_size: float = fig_track_size
        self.cross_link_color: str = cross_link_color
        self.inverted_cross_link_color: str = inverted_cross_link_color
        self.target_feature_types: List[str] = target_feature_types
        self.feature2color: Dict[str, str] = feature2color
        self.cds_limit_num: int = cds_limit_num

        self.gd = self._setup_genome_diagram()

    @property
    def max_range_length(self) -> int:
        """Max range length"""
        return max(gbk.range_length for gbk in self.gbk_list)

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

        for gbk in self.gbk_list:
            # Add track of one genbank
            gd_feature_set: FeatureSet = gd.new_track(
                track_level=0,
                name=gbk.name,
                greytrack=False,  # Disable greytrack
                greytrack_labels=0,
                greytrack_fontcolor=colors.black,
                start=0,
                end=gbk.range_length,
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

            range_features = gbk.extract_range_features(self.target_feature_types)
            if len(range_features) > self.cds_limit_num:
                continue

            for feature in range_features:
                # Make location fixed feature
                feature = SeqFeature(
                    location=FeatureLocation(
                        feature.location.start - gbk.min_range + 1,
                        feature.location.end - gbk.min_range + 1,
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
        cross_links = []
        for align_coord in self.align_coords:
            cross_link = align_coord.get_cross_link(
                tracks=gd.get_tracks(),
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
