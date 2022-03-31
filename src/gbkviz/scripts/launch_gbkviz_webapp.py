import argparse
import os
import subprocess as sp
from pathlib import Path

from gbkviz.__version__ import __version__


def main():
    """GBKviz main function for entrypoint"""
    args = get_args()
    port: int = args.port

    run(port)


def run(port: int = 8501):
    """Launch Streamlit GBKviz webapp

    Args:
        port (int): Port number to open web browser
    """
    # Streamlit env setting
    os.environ["STREAMLIT_THEME_BASE"] = "dark"
    os.environ["STREAMLIT_BROWSER_GATHER_USAGE_STATS"] = "false"

    # Launch Streamlit app
    gbkviz_dir = Path(__file__).parent.parent
    gbkviz_webapp_src_file = gbkviz_dir / "gbkviz_webapp.py"
    cmd = f"streamlit run {gbkviz_webapp_src_file} --server.port {port}"
    sp.run(cmd.split(" "))


def get_args():
    """Get arguments

    Returns:
        argparse.Namespace: Argument values
    """
    desc = "Simple web application to visualize and compare genomes in Genbank files"
    parser = argparse.ArgumentParser(description=desc)

    default_port = 8501
    parser.add_argument(
        "-p",
        "--port",
        type=int,
        help=f"Port number to open web server (Default: {default_port})",
        default=default_port,
        metavar="",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"v{__version__}",
        help="Print version information",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
