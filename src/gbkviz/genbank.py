from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import List, Union

import streamlit as st
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from streamlit.uploaded_file_manager import UploadedFile


@dataclass
class Genbank:
    """Genbank Class"""

    gbk_file: Union[str, StringIO, Path]
    name: str = ""

    def __post_init__(self):
        self._record: SeqRecord = list(SeqIO.parse(self.gbk_file, "genbank"))[0]

    @property
    def max_length(self) -> int:
        """Max genome sequence length"""
        return len(self._record.seq)

    def extract_features(self, target_features: List[str]) -> List[SeqFeature]:
        """Extract target features

        Args:
            target_features (List[str]): Target features type to extract

        Returns:
            List[SeqFeature]: Extract feature list

        Note:
             Target features: "CDS", "gene", "tRNA", "misc_feature"
        """
        return [f for f in self._record.features if f.type in target_features]

    @staticmethod
    @st.cache(allow_output_mutation=True)
    def read_upload_gbk_file(upload_gbk_file: UploadedFile) -> Genbank:
        """Read uploaded genbank file from Streamlit app

        Args:
            upload_gbk_file (UploadedFile): Uploaded genbank file

        Returns:
            Genbank: Genbank class object
        """
        return Genbank(
            gbk_file=StringIO(upload_gbk_file.getvalue().decode("utf-8")),
            name=upload_gbk_file.name,
        )
