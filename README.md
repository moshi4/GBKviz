# GBKviz: Genbank Data Visualization WebApp

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/moshi4/gbkviz/main/src/gbkviz/gbkviz_webapp.py)
![Python3](https://img.shields.io/badge/Language-Python3-steelblue)
![License](https://img.shields.io/badge/License-MIT-steelblue)
[![Latest PyPI version](https://img.shields.io/pypi/v/gbkviz.svg)](https://pypi.python.org/pypi/gbkviz)  

## Table of contents

- [Overview](#overview)
- [Install](#install)
- [Command Usage](#command-usage)
- [Genome Comparison](#genome-comparison)

## Overview

GBKviz is a web-based Genbank data visualization and comparison tool developed with streamlit web framework.
GBKviz allows user to easily and flexibly draw CDSs in user-specified genomic region (PNG or SVG format is available).
It also supports drawing genome comparison results by MUMmer.
This software is developed under the strong inspiration of [EasyFig](https://mjsull.github.io/Easyfig/).

![GBKviz Demo GIF](https://raw.githubusercontent.com/moshi4/GBKviz/main/src/gbkviz/gbkviz_demo.gif)  
Click [here](https://share.streamlit.io/moshi4/gbkviz/main/src/gbkviz/gbkviz_webapp.py) to try GBKviz on Streamlit Cloud.  
>:warning: Due to the limited resources in Streamlit Cloud, it may be unstable.  

## Install

GBKviz is implemented in Python3 (Tested on Ubuntu20.04)

Install PyPI stable version with pip:

    pip install gbkviz

If you want to enable genome comparison in GBKviz, MUMmer is required.  

Install MUMmer with apt command (Ubuntu):

    sudo apt install mummer

### Dependencies

- [Streamlit](https://streamlit.io/)  
  Simple web framework for data analysis

- [BioPython](https://github.com/biopython/biopython)  
  Utility tools for computational molecular biology

- [MUMmer](https://github.com/mummer4/mummer)  
  Genome alignment tool for comparative genomics
  
## Command Usage

Launch GBKviz in web browser (<http://localhost:8501>):

    gbkviz_webapp

## Genome Comparison

In GBKviz, [MUMmer](https://github.com/mummer4/mummer) is used as genome comparison tool.  
Following four genome comparison methods are available.

- Nucleotide One-to-One Mapping
- Nucleotide Many-to-Many Mapping
- Protein One-to-One Mapping
- Protein Many-to-Many Mapping

User can download and check genome comparison results file.  
Genome comparison results file is in the following tsv format.  

| Columns      | Contents                                            |
| ------------ | --------------------------------------------------- |
| REF_START    | Reference genome alignment start position           |
| REF_END      | Reference genome alignment end position             |
| QUERY_START  | Query genome alignment start position               |
| QUERY_END    | Query genome alignment end position                 |
| REF_LENGTH   | Reference genome alignment length                   |
| QUERY_LENGTH | Query genome alignment length                       |
| IDENTITY     | Reference and query genome alignment identity (%)   |
| REF_NAME     | Reference genome name tag                           |
| QUERY_NAME   | Query genome name tag                               |
