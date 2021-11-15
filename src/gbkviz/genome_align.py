from __future__ import annotations

import csv
import platform
import shutil
import subprocess as sp
from dataclasses import dataclass
from pathlib import Path
from typing import List, Union

import streamlit as st


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


@dataclass
class GenomeAlign:
    """Run MUMmer Genome Alingmnent Class"""

    genome_fasta_files: List[Union[str, Path]]
    seqtype: str = "nucleotide"  # "nucleotide" or "protein"
    maptype: str = "one-to-one"  # "one-to-one" or "many-to-many"

    def __post_init__(self):
        self.genome_fasta_files = [Path(f) for f in self.genome_fasta_files]
        self.seqtype = self.seqtype.lower()
        self.maptype = self.maptype.lower()

    @st.cache(allow_output_mutation=True)
    def run(self) -> List[AlignCoord]:
        """Run MUMmer genome alignment

        Returns:
            List[AlignCoords]: Genome alignment coordinates list
        """
        align_coords: List[AlignCoord] = []
        for idx in range(0, self.genome_num - 1):
            fa_file1 = self.genome_fasta_files[idx]
            fa_file2 = self.genome_fasta_files[idx + 1]
            outdir = Path(fa_file1).parent

            # Run nucmer or promer
            prefix = outdir / f"out{idx}"
            delta_file = prefix.with_suffix(".delta")
            cmd = f"{self.align_bin} {fa_file1} {fa_file2} --prefix={prefix}"
            sp.run(cmd, shell=True)

            # Run delta-filter
            filter_delta_file = outdir / f"filter_out{idx}.delta"
            cmd = f"delta-filter {self.filter_opt} {delta_file} > {filter_delta_file}"
            sp.run(cmd, shell=True)

            # Run show-coords
            coords_file = outdir / f"coords{idx}.tsv"
            cmd = f"show-coords -H -T {filter_delta_file} > {coords_file}"
            sp.run(cmd, shell=True)

            align_coords.extend(AlignCoord.read(coords_file, self.seqtype))

        return align_coords

    @property
    def genome_num(self) -> int:
        """Input genome fasta file count"""
        return len(self.genome_fasta_files)

    @property
    def align_bin(self) -> str:
        """Get genome alignment binary name"""
        if self.seqtype == "nucleotide":
            return "nucmer"
        elif self.seqtype == "protein":
            return "promer"
        else:
            raise ValueError(f"Invalid seqtype '{self.seqtype}'")

    @property
    def filter_opt(self) -> str:
        """Get filter option"""
        if self.maptype == "one-to-one":
            return "-1"
        elif self.maptype == "many-to-many":
            return "-m"
        else:
            raise ValueError(f"Invalid maptype '{self.maptype}'")

    @staticmethod
    def check_requirements() -> bool:
        """Check requirements to run genome align

        Returns:
            bool: Check result
        """
        # MacOS or Linux only
        if platform.system() not in ("Darwin", "Linux"):
            return False
        # Mummer binary installation check
        required_bins = ["nucmer", "promer", "delta-filter", "show-coords"]
        for required_bin in required_bins:
            if not shutil.which(required_bin):
                return False
        return True
