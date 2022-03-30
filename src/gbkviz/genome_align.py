import itertools
import multiprocessing as mp
import os
import platform
import shutil
import subprocess as sp
from pathlib import Path
from typing import List, Tuple, Union

import streamlit as st

from gbkviz.align_coord import AlignCoord


class GenomeAlign:
    """Run MUMmer Genome Alignment Class"""

    def __init__(
        self,
        genome_fasta_files: List[Path],
        outdir: Union[str, Path],
        seqtype: str = "nucleotide",
        maptype: str = "one-to-one",
    ):
        """GenomeAlign constructor

        Args:
            genome_fasta_files (List[Union[str, Path]]): Genome fasta files
            outdir (Union[str, Path]): Output directory
            seqtype (str, optional): "nucleotide" or "protein"
            maptype (str, optional): "one-to-one" or "many-to-many"
        """
        self.genome_fasta_files: List[Path] = [Path(f) for f in genome_fasta_files]
        self.outdir = Path(outdir)
        self.seqtype = seqtype.lower()
        self.maptype = maptype.lower()

    @st.cache(allow_output_mutation=True, ttl=3600)
    def run(self) -> List[AlignCoord]:
        """Run MUMmer genome alignment

        Returns:
            List[AlignCoords]: Genome alignment coordinates
        """
        # Prepare data for run MUMmer with multiprocessing
        mp_data_list: List[Tuple[Path, Path, int]] = []
        for idx in range(0, self.genome_num - 1):
            fa_file1 = self.genome_fasta_files[idx]
            fa_file2 = self.genome_fasta_files[idx + 1]
            mp_data_list.append((fa_file1, fa_file2, idx))

        # Run MUMmer with multiprocessing
        with mp.Pool(processes=self._process_num) as p:
            results = p.starmap(self._run_mummer, mp_data_list)

        return list(itertools.chain.from_iterable(results))

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

    @property
    def _process_num(self) -> int:
        cpu_num = os.cpu_count()
        return 1 if cpu_num is None or cpu_num == 1 else cpu_num - 1

    def _run_mummer(self, fa_file1: Path, fa_file2: Path, idx: int) -> List[AlignCoord]:
        """Run MUMmer function for multiprocessing

        Args:
            fa_file1 (Path): Input genome fasta 1
            fa_file2 (Path): Input genome fasta 2
            idx (int): Multiprocessing index

        Returns:
            List[AlignCoord]: AlignCoord list
        """
        # Run genome alignment using nucmer or promer
        prefix = self.outdir / f"out{idx}"
        delta_file = prefix.with_suffix(".delta")
        cmd = f"{self._align_bin} {fa_file1} {fa_file2} --prefix={prefix}"
        _ = sp.run(cmd, shell=True, capture_output=True, text=True)

        # Run delta-filter to map 'one-to-one' or 'many-to-many' relation
        filter_delta_file = self.outdir / f"filter_out{idx}.delta"
        cmd = f"delta-filter {self._map_opt} {delta_file} > {filter_delta_file}"
        _ = sp.run(cmd, shell=True, capture_output=True, text=True)

        # Run show-coords to extract alingment coords
        coords_file = self.outdir / f"coords{idx}.tsv"
        cmd = f"show-coords -H -T {filter_delta_file} > {coords_file}"
        _ = sp.run(cmd, shell=True, capture_output=True, text=True)

        align_coords = AlignCoord.parse(coords_file, self.seqtype)

        # Delete work files
        for work_file in (delta_file, filter_delta_file, coords_file):
            os.unlink(work_file)

        return align_coords

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
