HIPPOS Configuration Options
============================

There are two kinds of options, *essential* and *optional*. When nothing declared it means the option is essential. Some options have default value. Also keep in mind that **all file names must never use space!** So use underscore instead.

You can also inserting comments at the beginning of the line or after the option-value pair by inserting the ``#`` sign before the comment. Everything inserted after the ``#`` sign will be ignored by the software.

Basic Options
-------------

* **docking_method**
	value: ``vina`` or ``plants``

* **docking_conf**
	value: ``docking_configuration_file_name`` , eg. ``vina.conf`` or ``plants.conf``

* **direct_ifp**
	value: ``true`` or ``false``, (default: ``false``)

	If the value is set to true, it will use Direct IFP method and thus ``docking_method`` and ``docking_conf`` will be ignored.

Input Options
-------------

* **residue_name**
	value: ``list of residue_name`` , eg. ``ASP107 SER111 THR112``
	
	The list of residue_name, each residue separated by space. It is used in PLANTS post-analysis but not in VINA analysis as pdbqt don't hold the residue_name-residue_number pair field. However it is highly recommended to define it in VINA post analysis as it will be included in output file, making the results easier to interpret.
	
* **residue_number**
	value: ``list of residue_number`` , eg. ``80 84 85``
	
	The list of residue number, each number separated by space. Essential option in VINA post-analysis but optional in PLANTS post-analysis.
	
* **similarity_coef** (optional)
	value: ``list of similarity_coef`` , eg. ``tanimoto`` or ``mcconnaughey`` or ``tanimoto mcconnaughey``
	
	If the value is set then **full_ref** or **full_nobb_ref** or **simplified_ref** (depends on output_mode used)
	must be provided, so interaction fingerprint can be compared against reference.

	* **full_ref**
		value: ``bitstring1 bitstring2 ... n``
		
		At least 1 uniform bitstring must be provided

	* **full_nobb_ref**
		value: ``bitstring1 bitstring2 ... n``
		
		At least 1 uniform bitstring must be provided

	* **simplified_ref**
		value: ``bitstring1 bitstring2 ... n``
		
		At least 1 simplified bitstring must be provided
	
* **protein**
	value: ``protein file name`` , eg. ``protein.pdb``
	
	The protein file name that will be used in Direct IFP mode. Format allowed: ``pdb``, ``mol2``, ``pdbqt``.

* **ligand_files**
	value: ``list of ligand name`` , eg. ``3BO.pdb 3BR.pdb 3BS.pdb``
	
	The list of ligand name that will be paired with the protein in Direct IFP mode.

* **ligand_list**
	value: ``a file containing ligand file names`` , eg. ``ligand_list.txt``
	
	The name of the file that contain ligand file names.

* **multiple_ligand_list**
	value: ``list of file containing ligand file names`` , eg. ``ligand_list1.txt ligand_list2.txt ligand_list3.txt``
	
	The list of file that contain ligand file names.

* **complex_list**
	value: ``a filename containing list protein-ligand pair`` , eg. ``complex_list.txt``
	
	In Direct IFP mode, if this value is set then ``protein``, ``ligand_list``, and ``multiple_ligand_list`` will not be required.


Output Options
--------------

* **full_outfile** (optional)
	value: ``full_output_file_name`` , default: ``full_ifp.csv``
	
	Only used by ``full`` output_mode. It is recommended to use the csv extensions for clarity.

* **full_nobb_outfile** (optional)
	value: ``full_output_file_name`` , default: ``full_nobb_ifp.csv``
	
	Only used by ``full_nobb`` output_mode. It is recommended to use the csv extensions for clarity.

* **simplified_outfile** (optional)
	value: ``simplified_output_file_name`` , default: ``simplified_ifp.csv``
	
	Only used by ``simplified`` output_mode. It is recommended to use the csv extensions for clarity.

* **sim_outfile** (optional)
	value: ``similarity_output_file_name``, default: ``similarity.csv``
	
	Only used when the similarity coefficient is calculated.

* **logfile** (optional)
	value: ``log_file_name`` , default: ``hippos.log``
	
.. _advanced-options:
	
Advanced Options
----------------

* **output_mode** (optional)
	value: ``list of output_mode`` , default: ``full``, options: ``full``, ``full_nobb``, and ``simplified``
	
	The list of output_mode, at least one value required. When multiple value provided each value must be 
	separated by space. Multiple value can only be used in HIPPOS without reference. When nothing provided
	output_mode full will be used.
	
* **docking_score** (optional)
	value: ``yes`` or ``no`` , default: ``yes``
	
	Extract the docking score of ligand poses from docking results and attach them to output file.

* **omit_interaction** (optional)
	value: ``interaction_type`` and ``residue_name``

	where ``interaction_type`` is one of the following value:

	- ``hydrophobic`` or ``HPB``
	- ``aromatic`` or ``ARM``
	- ``h_bond`` or ``HBD``
	- ``electrostatic`` or ``ELE``
	- ``h_bond_donor`` or ``HBD_DON``
	- ``h_bond_acceptor`` or ``HBD_ACC``
	- ``electrostatic_positive`` or ``ELE_POS``
	- ``electrostatic_negative`` or ``ELE_NEG``
	- ``aromatic_facetoface`` or ``ARM_F2F``
	- ``aromatic_edgetoface`` or ``ARM_E2F``

	While ``residue_name`` specify which residue will be omitted. Usage example:

	``omit_interaction hydrophobic ARG223``

..
	* **res_weight1** (optional)
		value: ``residue_number residue_name interaction_type weight`` , eg. ``80 ASP107 electrostatic 5``
		
		Give weight to a spesific interaction on spesific amino acid residue. The example above shows that the number of electrostatic interaction bit on ASP107 will be multiplied by 5. The number 1 in **res_weight1** can be replaced with 2-5, Therefore there are 5 weight that can be applied to interaction fingerprinting.
