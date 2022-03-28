from pathlib import Path
from typing import List

import pytest


@pytest.fixture(scope="session")
def testdata_dir() -> Path:
    """testdata directory fixture"""
    return Path(__file__).parent / "testdata"


@pytest.fixture(scope="session")
def genbank_file(testdata_dir: Path) -> Path:
    """genbank files fixture"""
    genbank_dir = testdata_dir / "genbank"
    return genbank_dir / "JX128258.1.gbk"


@pytest.fixture(scope="session")
def genbank_files(testdata_dir: Path) -> List[Path]:
    """genbank files fixture"""
    genbank_dir = testdata_dir / "genbank"
    return list(genbank_dir.glob("*.gbk"))


@pytest.fixture(scope="session")
def genome_fasta_files(testdata_dir: Path) -> List[Path]:
    """genome fasta files fixture"""
    genome_fasta_dir = testdata_dir / "genome_fasta"
    return list(genome_fasta_dir.glob("*.fa"))
