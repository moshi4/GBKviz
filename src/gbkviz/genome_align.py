import platform
import shutil
import subprocess as sp
from dataclasses import dataclass
from pathlib import Path
from typing import List, Union

import streamlit as st

from gbkviz.align_coord import AlignCoord


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
            List[AlignCoords]: Genome alignment coordinates
        """
        align_coords: List[AlignCoord] = []
        for idx in range(0, self.genome_num - 1):
            fa_file1 = self.genome_fasta_files[idx]
            fa_file2 = self.genome_fasta_files[idx + 1]
            outdir = Path(fa_file1).parent

            # Run genome alignment using nucmer or promer
            prefix = outdir / f"out{idx}"
            delta_file = prefix.with_suffix(".delta")
            cmd = f"{self._align_bin} {fa_file1} {fa_file2} --prefix={prefix}"
            sp.run(cmd, shell=True)

            # Run delta-filter to map 'one-to-one' or 'many-to-many' relation
            filter_delta_file = outdir / f"filter_out{idx}.delta"
            cmd = f"delta-filter {self._map_opt} {delta_file} > {filter_delta_file}"
            sp.run(cmd, shell=True)

            # Run show-coords to extract alingment coords
            coords_file = outdir / f"coords{idx}.tsv"
            cmd = f"show-coords -H -T {filter_delta_file} > {coords_file}"
            sp.run(cmd, shell=True)

            align_coords.extend(AlignCoord.parse(coords_file, self.seqtype))

        return align_coords

    @property
    def genome_num(self) -> int:
        """Input genome fasta file count"""
        return len(self.genome_fasta_files)

    @property
    def _align_bin(self) -> str:
        """Get genome alignment program name ('nucmer' or 'promer')"""
        if self.seqtype == "nucleotide":
            return "nucmer"
        elif self.seqtype == "protein":
            return "promer"
        else:
            raise ValueError(f"Invalid seqtype '{self.seqtype}'")

    @property
    def _map_opt(self) -> str:
        """Get filter mapping option ('one-to-one' or 'many-to-many')"""
        if self.maptype == "one-to-one":
            return "-1"
        elif self.maptype == "many-to-many":
            return "-m"
        else:
            raise ValueError(f"Invalid maptype '{self.maptype}'")

    @staticmethod
    def check_requirements() -> bool:
        """Check requirements to run MUMmer genome alignment

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
