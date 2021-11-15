from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import List, Union


@dataclass
class AlignCoord:
    """Mummer Alignment Coordinates DataClass"""

    ref_start: int
    ref_end: int
    query_start: int
    query_end: int
    ref_length: int
    query_length: int
    identity: float
    ref_name: str
    query_name: str

    @staticmethod
    def read(
        coords_tsv_file: Union[str, Path],
        seqtype: str,
    ) -> List[AlignCoord]:
        """Read mummer(nucmer|promer) output coords result file"""
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
