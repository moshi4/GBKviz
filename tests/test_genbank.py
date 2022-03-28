from pathlib import Path

from gbkviz.genbank import Genbank


def test_full_length(genbank_file: Path):
    """test full length"""
    gbk = Genbank(genbank_file)
    actual_length, expected_length = gbk.full_length, 66854
    assert actual_length == expected_length


def test_range_length(genbank_file: Path):
    """test range length"""
    gbk = Genbank(genbank_file, "test", min_range=100, max_range=10099)
    actual_length, expected_length = gbk.range_length, 10000
    assert actual_length == expected_length


def test_extract_all_features(genbank_file: Path):
    """test extract all features"""
    gbk = Genbank(genbank_file)
    actual_feature_num, expected_feature_num = len(gbk.extract_all_features()), 94
    assert actual_feature_num == expected_feature_num


def test_extract_range_features(genbank_file: Path):
    """test extract range features"""
    gbk = Genbank(genbank_file, "test", min_range=1, max_range=1300)
    actual_feature_num, expected_feature_num = len(gbk.extract_range_features()), 4
    assert actual_feature_num == expected_feature_num


def test_write_genome_fasta_full(genbank_file: Path, tmp_path: Path):
    """test write genome fasta (full length)"""
    gbk = Genbank(genbank_file, "test")
    tmp_outfile = tmp_path / "tmp_genome_full.fna"
    gbk.write_genome_fasta(tmp_outfile)
    assert tmp_outfile.exists()


def test_write_genome_fasta_range(genbank_file: Path, tmp_path: Path):
    """test write genome fasta (range length)"""
    gbk = Genbank(genbank_file, "test", min_range=10000, max_range=20000)
    tmp_outfile = tmp_path / "tmp_genome_range.fna"
    gbk.write_genome_fasta(tmp_outfile, range=True)
    assert tmp_outfile.exists()
