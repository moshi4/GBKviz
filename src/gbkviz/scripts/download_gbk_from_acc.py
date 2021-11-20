#!/usr/bin/env python3
import argparse
from pathlib import Path
from urllib.error import HTTPError, URLError

from Bio import Entrez


def main():
    """Download genbank file from accession number"""
    # Get arguments
    args = get_args()
    acc_num: str = args.acc_num
    outfile: Path = args.outfile
    email: str = args.email

    # Register email
    Entrez.email = email

    # Download genbank file from accession number
    try:
        ret_text = Entrez.efetch(
            db="nucleotide", id=acc_num, rettype="gbwithparts", retmode="text"
        )
        gbk_text = ret_text.read()
    except HTTPError:
        print(f"AccessionNumber '{acc_num}' not found!!")
        exit(1)
    except URLError:
        print("Failed to download genbank file. Please confirm network connection.")
        exit(1)

    # Write genbank file
    with open(outfile, "w") as f:
        f.write(gbk_text)
    print(f"Donwload genbank file '{outfile}' (AccessionNumber={acc_num})")


def get_args() -> argparse.Namespace:
    """Get argument values

    Returns:
        argparse.Namespace: Argument values
    """
    parser = argparse.ArgumentParser(
        description="Get genbank file from accession number"
    )

    parser.add_argument(
        "-a",
        "--acc_num",
        required=True,
        type=str,
        help="Accession number",
        metavar="",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        required=True,
        type=Path,
        help="Output genbank file",
        metavar="",
    )
    parser.add_argument(
        "-e",
        "--email",
        type=str,
        default="address@example.com",
        help="Mail address (Optional)",
        metavar="",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
