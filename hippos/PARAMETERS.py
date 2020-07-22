'''
    Parameter for interaction distance
    Interaction exist if lower or equal 
    than these values
'''

HYDROPHOBIC = 4.5
AROMATIC = 4.0
HBOND = 3.5
ELECTROSTATIC = 4.0
'''
    Parameter for minimum H bond angle.
    Interaction exist if O --- H-D angle
    higher or equal than HBOND_ANGLE value.
'''

HBOND_ANGLE = 135
'''
    Parameter for aromatic interaction angle.
    Face to Face if:
    AROMATIC_ANGLE_LOW >= Angle between aromatic plane
    Or AROMATIC_ANGLE_HIGH <= Angle between aromatic plane
    Edge to Face if:
    AROMATIC_ANGLE_LOW < Angle between aromatic plane < 
    AROMATIC_ANGLE_HIGH
'''

AROMATIC_ANGLE_LOW = 30.0
AROMATIC_ANGLE_HIGH = 150.0
'''
    SIMP_IFP_MATRIX define which interaction of each
    amino acid residue that should be included and
    calculated. There are 6 flag (0 = OFF, 1 =  ON)
    for each residue. Each flag represent these interactions:
        1. Hydrophobic interaction
        2. Aromatic (Face to Face & Edge to Face)
        3. H-bond (Protein as Donor)
        4. H-bond (Protein as Acceptor)
        5. Electrostatic interaction (Protein +)
        6. Electrostatic interaction (Protein -)
    
    It is important to note, however, that
    only possible interaction can be given '1' flag.
    As a reference here is the matrix when all
    possible interaction flagged with 1:
        'ALA': (1,0,0,0,0,0),
        'CYS': (1,0,1,0,0,0),
        'ASP': (1,0,0,1,0,1),
        'GLU': (1,0,0,1,0,1),
        'PHE': (1,1,0,0,0,0),
        'GLY': (0,0,0,0,0,0),
        'HIS': (1,1,1,1,0,0),
        'ILE': (1,0,0,0,0,0),
        'LYS': (1,0,1,0,1,0),
        'LEU': (1,0,0,0,0,0),
        'MET': (1,0,0,0,0,0),
        'ASN': (1,0,1,1,0,0),
        'PRO': (1,0,0,0,0,0),
        'GLN': (1,0,1,1,0,0),
        'ARG': (1,0,1,0,1,0),
        'SER': (0,0,1,1,0,0),
        'THR': (1,0,1,1,0,0),
        'VAL': (1,0,0,0,0,0),
        'TRP': (1,1,1,0,0,0),
        'TYR': (1,1,1,1,0,0)
'''

SIMP_IFP_MATRIX = {
    'ALA': (1, 0, 0, 0, 0, 0),
    'CYS': (1, 0, 1, 0, 0, 0),
    'ASP': (0, 0, 0, 1, 0, 1),
    'GLU': (0, 0, 0, 1, 0, 1),
    'PHE': (0, 1, 0, 0, 0, 0),
    'GLY': (0, 0, 0, 0, 0, 0),
    'HIS': (0, 1, 1, 1, 0, 0),
    'ILE': (1, 0, 0, 0, 0, 0),
    'LYS': (0, 0, 1, 0, 1, 0),
    'LEU': (1, 0, 0, 0, 0, 0),
    'MET': (1, 0, 0, 0, 0, 0),
    'ASN': (0, 0, 1, 1, 0, 0),
    'PRO': (1, 0, 0, 0, 0, 0),
    'GLN': (0, 0, 1, 1, 0, 0),
    'ARG': (0, 0, 1, 0, 1, 0),
    'SER': (0, 0, 1, 1, 0, 0),
    'THR': (0, 0, 1, 1, 0, 0),
    'VAL': (1, 0, 0, 0, 0, 0),
    'TRP': (0, 1, 1, 0, 0, 0),
    'TYR': (0, 1, 1, 1, 0, 0),
}
