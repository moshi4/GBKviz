import subprocess as sp
from pathlib import Path


def main():
    """Launch Streamlit GBKviz webapp"""
    gbkviz_dir = Path(__file__).parent.parent
    st_gbkviz_src_file = gbkviz_dir / "gbkviz_webapp.py"
    sp.run(f"streamlit run {st_gbkviz_src_file}", shell=True)


if __name__ == "__main__":
    main()
