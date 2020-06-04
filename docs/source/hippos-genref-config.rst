HIPPOS-genref Configuration Options
===================================

There are two kinds of options, *essential* and *optional*. When nothing declared it means the option is essential. Some options have default value. Also keep in mind that **all file names must never use space!** So use underscore instead.

You can also inserting comments at the beginning of the line or after the option-value pair by inserting the ``#`` sign before the comment. Everything inserted after the ``#`` sign will be ignored by the software.

Input Options
-------------

* **residue_name**
	value: ``list of residue_name`` , eg. ``ASP107 SER111 THR112``
	
	The list of residue_name, each residue separated by space. It is used in PLANTS post-analysis but not in VINA analysis as pdbqt don't hold the residue_name-residue_number pair field. However it is highly recommended to define it in VINA post analysis as it will be included in output file, making the results easier to interpret.
	
* **residue_number**
	value: ``list of residue_number`` , eg. ``80 84 85``
	
	The list of residue number, each number separated by space. Essential option in VINA post-analysis but optional in PLANTS post-analysis.
	
* **proteins**
	value: ``list of reference_protein``, eg. ``protein1.mol2 protein2.mol2 protein3.mol2``
	
	The list of reference_protein structure, each reference separated by space. Must be in mol2 or pdbqt format.

* **ligands**
	value: ``list of reference_ligand`` , eg. ``ligand1.mol2 ligand2.mol2 ligand3.mol2``
	
	The list of reference_ligand structure, each reference separated by space. Must be in mol2 or pdbqt format.
	Please note that the first reference_protein will be paired with first reference_ligand, the second reference_protein
	will be paired with second reference_ligand, and so on.
	
Output Options
--------------

* **outfile** (optional)
	value: ``output_file_name``, default ``genref-results.txt``

Advanced Options
----------------

* **output_mode** (optional)
	value: ``output_mode`` , default: ``full``, options: ``full``, ``full_nobb``, and ``simplified``
	
	``output_mode`` define the bistring calculation and output. ``full`` means all interaction including the 
	backbone interaction taken into account. ``full_nobb`` means backbone interactions are omitted.
	``simplified`` means irrelevant interaction omitted (see :ref:`Simplified Bitstring Rule<simplified-rule>`).
	When this option not used ``output_mode full`` will be used.
