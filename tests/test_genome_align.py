from pathlib import Path
from typing import List

from gbkviz.genome_align import GenomeAlign


def test_genome_align_run_nucleotide(genome_fasta_files: List[Path], tmp_path: Path):
    """Test GenomeAlign run ('nucleotide' and 'one-to-one')"""
    genome_align = GenomeAlign(
        genome_fasta_files, tmp_path, seqtype="nucleotide", maptype="one-to-one"
    )
    align_coords = genome_align.run()
    assert len(align_coords) != 0


def test_genome_align_run_protein(genome_fasta_files: List[Path], tmp_path: Path):
    """Test GenomeAlign run ('nucleotide' and 'one-to-one')"""
    genome_align = GenomeAlign(
        genome_fasta_files, tmp_path, seqtype="protein", maptype="many-to-many"
    )
    align_coords = genome_align.run()
    assert len(align_coords) != 0
