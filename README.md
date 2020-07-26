# PyPLIF HIPPOS: A Molecular Interaction Fingerprinting Tool for Docking Results of AutoDock Vina and PLANTS

[![Build Status](https://travis-ci.com/radifar/PyPLIF-HIPPOS.svg?branch=master)](https://travis-ci.com/radifar/PyPLIF-HIPPOS)
[![Build status](https://ci.appveyor.com/api/projects/status/ai083k8ab0kxnk5w?svg=true)](https://ci.appveyor.com/project/radifar/pyplif-hippos)   
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/radifar/PyPLIF-HIPPOS.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/radifar/PyPLIF-HIPPOS/context:python)
[![codecov](https://codecov.io/gh/radifar/pyplif-hippos/branch/master/graph/badge.svg)](https://codecov.io/gh/radifar/pyplif-hippos/branch/master)
[![Documentation Status](https://readthedocs.org/projects/pyplif-hippos/badge/?version=latest&style=flat)](https://pyplif-hippos.readthedocs.io/en/latest/)
[![DOI:10.1021/acs.jcim.0c00305](https://zenodo.org/badge/DOI/10.1021/acs.jcim.0c00305.svg)](https://doi.org/10.1021/acs.jcim.0c00305)

<p align="center">
  <img alt="Icons made by Freepik from Flaticon is licensed by CC 3.0 BY" src="docs/source/hippopotamus_small.png">
</p>

<p align="center">Icons made by <a href="https://www.freepik.com/">Freepik</a> from <a href="http://www.flaticon.com">Flaticon</a> is licensed by CC 3.0 BY</p>

Welcome to PyPLIF-HIPPOS's project page. PyPLIF-HIPPOS is an upgraded version of [PyPLIF](https://github.com/radifar/pyplif/) (**Python-based Protein-Ligand Interaction Fingerprinting**), a tool for molecular docking post-analysis. It will translate the 3D coordinates of both ligand(s) (generated from docking simulation) and protein into a series of *interaction bitstring* (also known as *Interaction Fingerprint*) (see image below). **HIPPOS** (/ˌhipoʊz/) is a recursive acronym of **HIPPOS Is PyPLIF On Steroids**. From this point forward, PyPLIF-HIPPOS is simplified to HIPPOS.

Compared to PyPLIF, HIPPOS is not only faster and able to generate more customized interaction bitstring, but also supports both [PLANTS](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/pharmazie-und-biochemie/pharmazie/pharmazeutische-chemie/pd-dr-t-exner/research/plants/) & [Vina](http://vina.scripps.edu/)! More over, unlike its predecessor it is (far) more well-documented.

<p align="center">
  <img alt="Table of Content Abstract Graphic JCIM" src="docs/source/toc-abstract-graphics_small.png">
</p>

<p align="center">Reprinted with permission from https://doi.org/10.1021/acs.jcim.0c00305. Copyright 2020 American Chemical Society.</p>

<p align="center">
  <img alt="PyPLIF output from PyPLIF publication" src="docs/source/pyplif.png">
</p>

<p align="center">Illustration by Radifar et al (2013) from <a href="http://www.bioinformation.net/009/97320630009325.htm">Bioinformation.net</a> is licensed by <a href="http://creativecommons.org/licenses/by/4.0">CC 4.0 BY</a>

## Citing HIPPOS

If you are using HIPPOS please cite this paper:

Istyastono, E., Radifar, M., Yuniarti, N., Prasasty, V. and Mungkasi, S., 2020.
PyPLIF HIPPOS: A Molecular Interaction Fingerprinting Tool for Docking Results
of AutoDock Vina and PLANTS. Journal of Chemical Information and Modeling.
https://doi.org/10.1021/acs.jcim.0c00305

## Acknowledgment

This project has received funding from the [Indonesian National Research and Innovation Agency](https://international.ristekdikti.go.id/) 
under grant agreement No. 807.7/LL5/PG/2020. This project has been restructured based on the
[MOLSSI Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms)
version 1.3, and benefited greatly from [MOLSSI Python Package Development Best Practices](https://molssi.org/2020/04/20/may-webinar-series-python-package-development/)
workshop.


-----

&copy; Copyright 2020, Muhammad Radifar & Enade Perdana Istyastono
