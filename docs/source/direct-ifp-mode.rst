Direct IFP Mode
===============


Although at first PyPLIF-HIPPOS is originally intended for analyzing the protein-ligand
interaction from docking tools, sometimes we would like to analyze protein-ligand 
interaction from various computational chemistry software. For example sometimes we just
want to analyze the protein-ligand interaction from MD trajectories, or maybe when we
would like to analyze the output of a certain docking tool. Direct IFP mode is a feature
that allow these kind of analysis.

Note that at the moment, Direct IFP mode will only accept ``mol2``, ``pdb``, and ``pdbqt``
file format.

Analyzing Protein-Ligand Interaction from one protein and multiple ligand pose
------------------------------------------------------------------------------

In this example we would use few ways in which we could analyze the interaction between
several ligands against a single protein. This example is suitable for analyzing the
docking results of a docking tool that is not supported by PyPLIF-HIPPOS.

For the protein input we simply use the protein keyword followed by the protein file name.
As for the ligand, there are basically two ways to read the ligand file names. The first one
is by directly mentioning the ligand file names, the second way is by mentioning file names
where each file is a text file that contain list of ligand file names.

First lets take a look on how we can do this via mentioning the ligand file names.
This can be done by setting three options in the configuration file: ``direct_ifp``, 
``protein``, and ``ligand_files``. The example files for this can be found in
``examples-input/09-direct_ifp/example_0/pdb`` directory. ::

    direct_ifp true
    protein 3chr_protein.pdb
    ligand_files 3chp_4BO.pdb 3chq_4BQ.pdb 3chr_4BS.pdb 3chs_4BU.pdb

    similarity_coef   tanimoto mcconnaughey

    full_ref  0000000000000000000001000000000010000000000000000000000010100000000000000000000000001000000000000000000000010000000000000000000000000

    residue_name GLN134 GLN136 ALA137 TYR267 GLY269 MET270 GLU271 TRP311 PHE314 GLU318 VAL367 LEU369 PRO374 ASP375 ALA377 TYR378 SER379 PRO382 TYR383
    residue_number 134 136 137 267 269 270 271 311 314 318 367 369 374 375 377 378 379 382 383

    full_outfile example_0_outfile.csv
    sim_outfile example_0_similarity.csv

The other way to read multiple ligand files is using expression pattern like so: ::

    direct_ifp true
    protein protein.pdb
    ligand_files 3ch*.pdb

    similarity_coef   tanimoto mcconnaughey

    full_ref  0000000000000000000001000000000010000000000000000000000010100000000000000000000000001000000000000000000000010000000000000000000000000

    residue_name GLN134 GLN136 ALA137 TYR267 GLY269 MET270 GLU271 TRP311 PHE314 GLU318 VAL367 LEU369 PRO374 ASP375 ALA377 TYR378 SER379 PRO382 TYR383
    residue_number 134 136 137 267 269 270 271 311 314 318 367 369 374 375 377 378 379 382 383

    full_outfile example_0_outfile.csv
    sim_outfile example_0_similarity.csv

The above example can be found in ``examples-input/09-direct_ifp/example_1/pdb`` directory. In above example,
the result will be the same as ``example_0`` which will generate two output file, one for the interaction bitstring
``example_0_outfile.csv``, and the other one for similarity score ``example_0_similarity.csv``. Here is the example
output for ``example_0_outfile.csv`` ::

    3chp_4BO.pdb               0000000000000000000001000000000000000000000000000100000010100000000000000000010000001000000000000000000001010000000000000000000001000
    3chq_4BQ.pdb               0000000000000000000001000100000000000000000000001000000000100000000001000000010000001000000000000000000001010000000000000000000001000
    3chr_4BS.pdb               0000000000000000000001000000000000010000000000101100000010100000000001000000010000001000000000000000000000010000000000000000000000000
    3chs_4BU.pdb               0000000000010000000001000000000000000000000000101000000010100000000001000000010000001000000000000010000001010000000000000000000001000

And here is the example output for ``example_0_similarity.csv`` :: 

    3chp_4BO.pdb     0.500 0.389
    3chq_4BQ.pdb     0.333 0.067
    3chr_4BS.pdb     0.417 0.288
    3chs_4BU.pdb     0.357 0.218

The second one is by providing the file containing a list of ligand file names. To do this we can use the
``ligand_list`` or ``multiple_ligand_list`` option instead of ``ligand_files``. ``ligand_list`` will accept a single
text file containing the ligand file names like so  ::

    direct_ifp true
    protein protein.pdb
    ligand_list ligand_input.txt

    similarity_coef   tanimoto mcconnaughey

    full_ref  0000000000000000000001000000000010000000000000000000000010100000000000000000000000001000000000000000000000010000000000000000000000000

    residue_name GLN134 GLN136 ALA137 TYR267 GLY269 MET270 GLU271 TRP311 PHE314 GLU318 VAL367 LEU369 PRO374 ASP375 ALA377 TYR378 SER379 PRO382 TYR383
    residue_number 134 136 137 267 269 270 271 311 314 318 367 369 374 375 377 378 379 382 383

    full_outfile example_0_outfile.csv
    sim_outfile example_0_similarity.csv

Where the ``ligand_input.txt`` content is ::

    3chp_4BO.pdb
    3chq_4BQ.pdb
    3chr_4BS.pdb
    3chs_4BU.pdb

While ``multiple_ligand_list`` can accept multiple file that contain ligand file names like so ::

    direct_ifp true
    protein protein.pdb
    multiple_ligand_list ligand_input1.txt ligand_input2.txt

    similarity_coef   tanimoto mcconnaughey

    full_ref  0000000000000000000001000000000010000000000000000000000010100000000000000000000000001000000000000000000000010000000000000000000000000

    residue_name GLN134 GLN136 ALA137 TYR267 GLY269 MET270 GLU271 TRP311 PHE314 GLU318 VAL367 LEU369 PRO374 ASP375 ALA377 TYR378 SER379 PRO382 TYR383
    residue_number 134 136 137 267 269 270 271 311 314 318 367 369 374 375 377 378 379 382 383

    full_outfile example_0_outfile.csv
    sim_outfile example_0_similarity.csv

The example files for ``ligand_list`` and ``multiple_ligand_list`` can be found in ``examples-input/09-direct_ifp/example_2/pdb``
and ``examples-input/09-direct_ifp/example_3/pdb``. And the output of these options will be similar to previous
examples.

Analyzing Protein-Ligand Interaction from multiple protein-ligand complex
-------------------------------------------------------------------------

In this example we will use three different protein-ligand complex, which could represent
a simple MD trajectory. Therefore this kind of method is suitable for MD trajectory analysis.

In order to analyze several pair of protein-ligand we can use ``direct_ifp`` and ``complex_list`` option.
Notice that in this example we will not using ``protein`` option since the protein already included in the
complex_list file.

First lets take a look at the ``complex_list.txt`` which contain the protein-ligand pair that will be analyzed ::

    3cho_protein.pdb 3cho_4BG.pdb
    3chr_protein.pdb 3chr_4BS.pdb
    3chs_protein.pdb 3chs_4BU.pdb

This example files can be found in ``examples-input/09-direct_ifp/example_4/pdb``. Next we can use the following
configuration file to analyze the above protein-ligand pairs ::

    direct_ifp true
    complex_list  complex_list.txt

    similarity_coef   tanimoto mcconnaughey

    full_ref  0000000000000000000001000000000010000000000000000000000010100000000000000000000000001000000000000000000000010000000000000000000000000

    residue_name GLN134 GLN136 ALA137 TYR267 GLY269 MET270 GLU271 TRP311 PHE314 GLU318 VAL367 LEU369 PRO374 ASP375 ALA377 TYR378 SER379 PRO382 TYR383
    residue_number 134 136 137 267 269 270 271 311 314 318 367 369 374 375 377 378 379 382 383

    full_outfile example_0_outfile.csv
    sim_outfile example_0_similarity.csv

Running the above example will give us the following output for ``example_0_outfile.csv`` ::

    3cho_protein.pdb_3cho_4BG.pdb
            0000000000000000000001000000000010000000000000000000000010100000000000000000000000001000000000000000000000010000000000000000000000000
    3chr_protein.pdb_3chr_4BS.pdb
            0000000000000000000001000000000000010000000000101100000010100000000001000000010000001000000000000000000000010000000000000000000000000
    3chs_protein.pdb_3chs_4BU.pdb
            0000000000010000000001000000000000010000000000101100000010100000000001000000010000001000000000000000000001010000000000000000000001000


and ``example_0_similarity.csv`` ::
    
    3cho_protein.pdb_3cho_4BG.pdb
    1.000 1.000
    3chr_protein.pdb_3chr_4BS.pdb
    0.417 0.288
    3chs_protein.pdb_3chs_4BU.pdb
    0.333 0.190
