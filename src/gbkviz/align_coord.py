from __future__ import annotations

import csv
from dataclasses import dataclass
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
        minus_bp: int = 0,
        normal_color: str = "#E80F13",  # Red
        inverted_color: str = "#0F0FE8",  # Blue
    ) -> CrossLink:
        # Change cross link bp
        ref_start = self.ref_start - minus_bp
        ref_end = self.ref_end - minus_bp
        query_start = self.query_start - minus_bp
        query_end = self.query_end - minus_bp

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
            # flip=flip,
        )

    @property
    def is_inverted(self) -> bool:
        """Check inverted alignment coord or not"""
        return (self.ref_end - self.ref_start) * (self.query_end - self.query_start) < 0

    @staticmethod
    def parse(
        coords_tsv_file: Union[str, Path],
        seqtype: str,
    ) -> List[AlignCoord]:
        """Parse MUMmer(nucmer|promer) output coords result file"""
        align_coords = []
        with open(coords_tsv_file) as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                # Check read file contents & extract required row values
                if seqtype == "nucleotide":
                    if len(row) != 9:
                        raise ValueError("Invalid Mummer(nucmer) coords file!!")
                elif seqtype == "protein":
                    if len(row) != 13:
                        raise ValueError("Invalid Mummer(promer) coords file!!")
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
