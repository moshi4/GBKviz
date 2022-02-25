from __future__ import annotations

import csv
from dataclasses import astuple, dataclass
from pathlib import Path
from typing import Dict, List, Union

from Bio.Graphics.GenomeDiagram import CrossLink, Track
from reportlab.lib import colors
from reportlab.lib.colors import HexColor


@dataclass
class AlignCoord:
    """MUMmer Alignment Coordinates DataClass"""

    ref_start: int
    ref_end: int
    query_start: int
    query_end: int
    ref_length: int
    query_length: int
    identity: float
    ref_name: str
    query_name: str

    def get_cross_link(
        self,
        name2track: Dict[str, Track],
        name2start: Dict[str, int],
        normal_color: str = "#0000FF",  # Blue
        inverted_color: str = "#FF0000",  # Red
    ) -> CrossLink:
        """Get cross link object for genome comparison visualization

        Args:
            name2track (Dict[str, Track]): Name and Track dictionary
            name2start (Dict[str, int]): Name and Start(bp) dictionary
            normal_color (str): Normal cross link hexcolor (Default='#0000FF'[Blue])
            inverted_color (str): Inverted cross link hexcolor (Default='#FF0000'[Red])

        Returns:
            CrossLink: Cross link object
        """
        # Get cross link start-end of reference and query
        ref_adjust_bp = name2start[self.ref_name]
        query_adjust_bp = name2start[self.query_name]
        ref_start = min(self.ref_start, self.ref_end) - ref_adjust_bp + 1
        ref_end = max(self.ref_start, self.ref_end) - ref_adjust_bp + 1
        query_start = min(self.query_start, self.query_end) - query_adjust_bp + 1
        query_end = max(self.query_start, self.query_end) - query_adjust_bp + 1

        # GenomeDiagram cannot draw cross link color correctly in condition below
        # To resolve this drawing error, add 1 bp length to ref_start
        if self.ref_length == self.query_length and self.is_inverted:
            ref_start += 1

        # Set cross link color
        if self.is_inverted:
            cross_link_color = HexColor(inverted_color)
        else:
            cross_link_color = HexColor(normal_color)

        # Get gradient color from alignment sequence identity[%]
        gradient_cross_link_color = colors.linearlyInterpolatedColor(
            colors.white, cross_link_color, 0, 100, self.identity
        )

        return CrossLink(
            featureA=(name2track[self.ref_name], ref_start, ref_end),
            featureB=(name2track[self.query_name], query_start, query_end),
            color=gradient_cross_link_color,
            border=gradient_cross_link_color,
            flip=self.is_inverted,
        )

    @property
    def is_inverted(self) -> bool:
        """Check inverted alignment coord or not"""
        return (self.ref_end - self.ref_start) * (self.query_end - self.query_start) < 0

    @property
    def as_tsv_format(self) -> str:
        """TSV format text"""
        return "\t".join([str(v) for v in astuple(self)])

    @staticmethod
    def parse(
        coords_tsv_file: Union[str, Path],
        seqtype: str,
    ) -> List[AlignCoord]:
        """Parse MUMmer(nucmer|promer) output coords result file

        Args:
            coords_tsv_file (Union[str, Path]): MUMmer align coords file
            seqtype (str): Sequence type ('nucleotide' or 'protein')

        Returns:
            List[AlignCoord]: Align coords
        """
        align_coords = []
        with open(coords_tsv_file) as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                # Check read file contents & extract required row values
                if seqtype == "nucleotide":
                    if len(row) != 9:
                        raise ValueError(
                            f"Invalid nucmer coords file '{coords_tsv_file}'!!"
                        )
                elif seqtype == "protein":
                    if len(row) != 13:
                        raise ValueError(
                            f"Invalid promer coords file '{coords_tsv_file}'!!"
                        )
                    row = row[0:7] + row[11:13]
                else:
                    raise ValueError(f"Invalid seqtype '{seqtype}'!!")

                # Convert to correct value type
                typed_row = []
                for idx, val in enumerate(row):
                    if 0 <= idx <= 5:
                        typed_row.append(int(val))
                    elif idx == 6:
                        typed_row.append(float(val))
                    else:
                        typed_row.append(str(val))

                align_coords.append(AlignCoord(*typed_row))

        return align_coords
