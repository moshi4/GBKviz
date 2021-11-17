import os
import subprocess as sp
from pathlib import Path


def main():
    """Launch Streamlit GBKviz webapp"""
    # Streamlit env setting
    os.environ["STREAMLIT_THEME_BASE"] = "dark"

    # Launch Streamlit app
    gbkviz_dir = Path(__file__).parent.parent
    gbkviz_webapp_src_file = gbkviz_dir / "gbkviz_webapp.py"
    sp.run(f"streamlit run {gbkviz_webapp_src_file}", shell=True)


if __name__ == "__main__":
    main()
