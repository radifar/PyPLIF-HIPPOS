.. HIPPOS documentation master file, created by
   sphinx-quickstart on Fri Dec  7 13:14:20 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyPLIF HIPPOS's documentation!
=========================================

PyPLIF HIPPOS: A Molecular Interaction Fingerprinting Tool for Docking Results of AutoDock Vina and PLANTS
----------------------------------------------------------------------------------------------------------

.. image:: hippopotamus.png
	:alt: Icons made by Freepik from Flaticon is licensed by CC 3.0 BY
	:align: center
	:scale: 45%
	
.. raw:: html
	:file: attribution-hippos.html

Welcome to PyPLIF-HIPPOS's documentation. PyPLIF-HIPPOS is an upgraded version of `PyPLIF <https://github.com/radifar/pyplif/>`_ (**Python-based Protein-Ligand Interaction Fingerprinting**), a tool for molecular docking post-analysis. It will translate the 3D coordinates of both ligand(s) (generated from docking simulation) and protein into a series of *interaction bitstring* (also known as *Interaction Fingerprint*) (see image below). **HIPPOS** (/ˌhipoʊz/) is a recursive acronym of **HIPPOS Is PyPLIF On Steroids**. From this point forward, PyPLIF-HIPPOS is simplified to HIPPOS.

Compared to PyPLIF, HIPPOS is not only faster and able to generate more customized interaction bitstring, but also supports both `PLANTS <https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/pharmazie-und-biochemie/pharmazie/pharmazeutische-chemie/pd-dr-t-exner/research/plants/>`_ & `VINA <http://vina.scripps.edu/>`_! More over, unlike its predecessor it is (far) more well-documented. 

.. image:: toc-abstract-graphics_small.png
	:alt: Table of Content Abstract Graphic JCIM
	:align: center

Reprinted with permission from https://doi.org/10.1021/acs.jcim.0c00305. Copyright 2020 American Chemical Society.

.. image:: pyplif.png
	:alt: PyPLIF output from PyPLIF publication
	:align: center

.. raw:: html
	:file: attribution-pyplif.html

Citing HIPPOS
-------------

If you are using HIPPOS please cite this paper:

Istyastono, E., Radifar, M., Yuniarti, N., Prasasty, V. and Mungkasi, S., 2020.
PyPLIF HIPPOS: A Molecular Interaction Fingerprinting Tool for Docking Results
of AutoDock Vina and PLANTS. Journal of Chemical Information and Modeling.
https://doi.org/10.1021/acs.jcim.0c00305

Acknowledgment
--------------

This project has received funding from the `Indonesian National Research and Innovation Agency <https://international.ristekdikti.go.id/>`_
under grant agreement No. 807.7/LL5/PG/2020. This project has been restructured based on the
`MOLSSI Computational Molecular Science Python Cookiecutter <https://github.com/molssi/cookiecutter-cms>`_
version 1.3, and benefited greatly from `MOLSSI Python Package Development Best Practices <https://molssi.org/2020/04/20/may-webinar-series-python-package-development/>`_
workshop.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   installation
   parameter-options
   getting-started-genref
   hippos-genref-config
   getting-started-vina
   getting-started-plants
   hippos-config
   license

