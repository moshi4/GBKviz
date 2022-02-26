from __future__ import annotations

from io import StringIO
from pathlib import Path
from typing import List, Optional, Union

import streamlit as st
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from streamlit.uploaded_file_manager import UploadedFile


class Genbank:
    """Genbank Class"""

    def __init__(
        self,
        gbk_file: Union[str, StringIO, Path],
        name: str = "",
        min_range: Optional[int] = None,
        max_range: Optional[int] = None,
        reverse: bool = False,
    ):
        """Genbank constructor

        Args:
            gbk_file (Union[str, StringIO, Path]): Genbank file
            name (str, optional): Name
            min_range (Optional[int], optional): Min range
            max_range (Optional[int], optional): Max range
            reverse (bool, optional): Reverse or not
        """
        self._record: SeqRecord = list(SeqIO.parse(gbk_file, "genbank"))[0]
        self.name: str = name
        self.min_range: int = 1 if min_range is None else min_range
        self.max_range: int = len(self._record.seq) if max_range is None else max_range
        self.reverse: bool = reverse

    @property
    def full_length(self) -> int:
        """Whole genome sequence length"""
        return len(self.record.seq)

    @property
    def range_length(self) -> int:
        """Range genome sequence length"""
        return self.max_range - self.min_range + 1

    @property
    def record(self) -> SeqRecord:
        """Genbank record"""
        if self.reverse is True:
            return self._record.reverse_complement()
        else:
            return self._record

    def extract_all_features(
        self,
        feature_types: List[str] = ["CDS"],
    ) -> List[SeqFeature]:
        """Extract all features

        Args:
            feature_types (List[str]): Feature types to extract

        Returns:
            List[SeqFeature]: All features
        """
        return [f for f in self.record.features if f.type in feature_types]

    def extract_range_features(
        self,
        feature_types: List[str] = ["CDS"],
    ) -> List[SeqFeature]:
        """Extract features in range

        Args:
            feature_types (List[str]): Feature types to extract

        Returns:
            List[SeqFeature]: Features in range
        """
        range_features = []
        for feature in self.extract_all_features(feature_types):
            min_feature_loc = int(min(feature.location.start, feature.location.end))
            max_feature_loc = int(max(feature.location.start, feature.location.end))
            if (
                self.min_range <= min_feature_loc <= self.max_range
                or self.min_range <= max_feature_loc <= self.max_range
            ):
                range_features.append(feature)

        return range_features

    def write_genome_fasta(
        self,
        outfile: Union[str, Path],
        range: bool = False,
    ) -> None:
        """Write genome fasta file

        Args:
            outfile (Union[str, Path]): Output genome fasta file
            range (bool): Write range genome or full genome
        """
        if range:
            write_seq = self.record.seq[self.min_range - 1 : self.max_range - 1]
        else:
            write_seq = self.record.seq
        with open(outfile, "w") as f:
            f.write(f">{self.name}\n{write_seq}\n")

    @staticmethod
    @st.cache(allow_output_mutation=True, ttl=3600)
    def read_upload_gbk_file(upload_gbk_file: UploadedFile) -> Genbank:
        """Read uploaded genbank file from Streamlit app

        Args:
            upload_gbk_file (UploadedFile): Uploaded genbank file

        Returns:
            Genbank: Genbank class object
        """
        return Genbank(
            gbk_file=StringIO(upload_gbk_file.getvalue().decode("utf-8")),
            name=Path(upload_gbk_file.name).stem,
        )
