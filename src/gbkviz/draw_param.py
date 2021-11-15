from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass()
class DrawParam:
    """Draw Parameter DataClass"""

    fig_format: str = "JPG"  # "JPG", "PNG", "SVG", "PDF"
    show_label: bool = False
    show_scale: bool = True
    show_ticks: bool = True
    label_type: str = "gene"  # "gene", "protein_id", "locus_tag", "product"
    feature_symbol: str = "BIGARROW"  # "BIGARROW", "BOX", "ARROW", "OCTO"
    label_angle: int = 30  # 0, 15, 30, 45, 60, 75, 90
    scaleticks_interval: int = 10000  # 1Kbp, 5Kbp, 10Kbp, 50Kbp, 100Kbp, 1Mbp
    label_fsize: int = 10  # Range: 0 - 100 (Step: 1)
    scaleticks_fsize: int = 8  # Range: 0 - 100 (Step: 1)
    fig_width: int = 25  # Range: 10 - 100 (Step: 5)
    fig_track_height: int = 3  # Range: 1 - 10 (Step: 1)
    fig_track_size: float = 0.5  # Range: 0.1 - 1.0 (Step: 0.1)
    target_features: List[str] = field(  # "CDS", "gene", "tRNA", "misc_feature"
        default_factory=lambda: ["CDS"],
    )
    feature2color: Dict[str, str] = field(
        default_factory=lambda: {
            "CDS": "#FFA500",
            "gene": "#0FE8E4",
            "tRNA": "#E80F0F",
            "misc_future": "#E80FC6",
        }
    )
    # "Nucleotide|Protein" * "One-to-One|Multi-to-Multi"
    genome_comparison: Optional[str] = None
