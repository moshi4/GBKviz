from pathlib import Path
from typing import List

from gbkviz.draw_genbank_fig import DrawGenbankFig
from gbkviz.genbank import Genbank


def test_draw_genbank_fig(genbank_files: List[Path], tmp_path: Path):
    """test draw_genbank_fig"""
    gbk_list = [Genbank(gf, gf.name) for gf in genbank_files]
    gdf = DrawGenbankFig(gbk_list)
    fig_png_outfile = tmp_path / "fig.png"
    fig_svg_outfile = tmp_path / "fig.svg"
    gdf.write_figure(fig_png_outfile)
    gdf.write_figure(fig_svg_outfile)
    assert fig_png_outfile.exists() and fig_svg_outfile.exists()
