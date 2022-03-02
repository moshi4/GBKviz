from __future__ import annotations

import csv
from dataclasses import astuple, dataclass
from pathlib import Path
from typing import List, Union

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
        tracks: List[Track],
        normal_color: str = "#0000FF",  # Blue
        inverted_color: str = "#FF0000",  # Red
    ) -> CrossLink:
        """Get cross link object for genome comparison visualization

        Args:
            tracks (List[Track]): GenomeDiagram Track list
            normal_color (str): Normal cross link hexcolor (Default='#0000FF'[Blue])
            inverted_color (str): Inverted cross link hexcolor (Default='#FF0000'[Red])

        Returns:
            CrossLink: Cross link object
        """
        # Get cross link start-end of reference and query
        ref_start = min(self.ref_start, self.ref_end)
        ref_end = max(self.ref_start, self.ref_end)
        query_start = min(self.query_start, self.query_end)
        query_end = max(self.query_start, self.query_end)

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
        upper_limit = 100
        lower_limit = 20 if self.identity > 20 else self.identity
        gradient_cross_link_color = colors.linearlyInterpolatedColor(
            colors.white, cross_link_color, lower_limit, upper_limit, self.identity
        )

        name2track = {track.name: track for track in tracks}
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

    def add_offset(self, ref_offset: int, query_offset: int) -> AlignCoord:
        """Add offset to start-end position

        Args:
            ref_offset (int): Offset for reference genome
            query_offset (int): Offset for query genome

        Returns:
            AlignCoord: AlignCoord with offset
        """
        return AlignCoord(
            self.ref_start + ref_offset,
            self.ref_end + ref_offset,
            self.query_start + query_offset,
            self.query_end + query_offset,
            self.ref_length,
            self.query_length,
            self.identity,
            self.ref_name,
            self.query_name,
        )

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
