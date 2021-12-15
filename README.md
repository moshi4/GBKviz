# GBKviz: Genbank Data Visuzalization WebApp

![Python3](https://img.shields.io/badge/Language-Python3-steelblue)
![License](https://img.shields.io/badge/License-MIT-steelblue)
[![Latest PyPI version](https://img.shields.io/pypi/v/gbkviz.svg)](https://pypi.python.org/pypi/gbkviz)  
[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/moshi4/gbkviz/main/src/gbkviz/gbkviz_webapp.py)

## Overview

GBKviz is browser-based Genbank data vizualization webapp.  
User can easily and flexibly plot specified region genome CDS figure in any format.  
In addition, GBKvis can perform comparative genome plot figure.  

Demo GIF here...

## Install

GBKviz is implemented with Python3 (Tested on Ubuntu20.04)

Install PyPI stable version with pip:

    pip install gbkviz

Install latest development version with pip:

    pip install git+git://github.com/moshi4/GBKviz.git

### Dependencies

- [Streamlit](https://streamlit.io/)  
  Web app framework for quick development

- [BioPython](https://github.com/biopython/biopython)  
  Utility tools for computational molecular biology

- [MUMmer](https://github.com/mummer4/mummer) (Optional)  
  Genome alignment tool for comparative genomics (v3 or v4)
  
## Command Usage

Launch GBKviz in web browser:

    gbkviz_webapp

## Usage
