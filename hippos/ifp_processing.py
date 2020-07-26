import numpy as np
from bitarray import bitarray
from openbabel import OBMol, OBConversion, OBResidueIter, OBResidueAtomIter, OBMolAtomIter

from PARAMETERS import HYDROPHOBIC, AROMATIC, HBOND, ELECTROSTATIC, HBOND_ANGLE,\
                       AROMATIC_ANGLE_LOW, AROMATIC_ANGLE_HIGH, SIMP_IFP_MATRIX
'''
    function:
    - get_bitstring --> return residue_obj
    - get_refbitstring --> return residue_obj_list (list --> multiple reference)
    - assign_atoms --> return atom_group (dictionary),
                       value: atom indices & possible interactions
    - getRing --> return list of ring (from GetSSSR), used in assign_atoms
    - getCrossModulus --> calculate the cross and modulus for ring angle calculation
    - Residue class
'''


def get_bitstring(docking_results, hippos_config):
    protein = docking_results['protein']
    res_name = hippos_config['residue_name']
    res_num = hippos_config['residue_number']
    ligands = docking_results['docked_ligands']
    flex_protein = docking_results['docked_proteins']

    custom_settings = {
        'backbone': hippos_config['use_backbone'],
        'output_mode': hippos_config['output_mode'],
        'res_weight1': hippos_config['res_weight1'],
        'res_weight2': hippos_config['res_weight2'],
        'res_weight3': hippos_config['res_weight3'],
        'res_weight4': hippos_config['res_weight4'],
        'res_weight5': hippos_config['res_weight5'],
    }
    '''
    Assigning residue_obj
    key   => residue_name, eg. ASP107
    value => Residue object
    '''

    residue_obj = {}
    for name, num in zip(res_name, res_num):
        residue_obj[name] = Residue(protein, name, num, custom_settings)
    '''
    Calculate IFP for each Residue in residue_obj
    Calculate IFP using
    docked_ligands => list of OBMol
    ligand_atom_group => Dictionary of atom indices & possible interaction
    '''

    if hippos_config['docking_method'] == 'vina':
        ligand_atom_group = assign_atoms(ligands[0], 'vina')
        for name in res_name:
            residue_obj[name].calculateIFPVina(ligands, ligand_atom_group)
    else:
        ligand_atom_group = assign_atoms(ligands[0], 'plants')
        for name in res_name:
            residue_obj[name].calculateIFPPlants(ligands, flex_protein, ligand_atom_group)

    return residue_obj


def get_refbitstring(genref_config):
    proteins = genref_config['proteins']
    ligands = genref_config['ligands']
    res_name = genref_config['residue_name']
    res_num = genref_config['residue_number']

    custom_settings = {
        'backbone': genref_config['use_backbone'],
        'res_weight1': genref_config['res_weight1'],
        'res_weight2': genref_config['res_weight2'],
        'res_weight3': genref_config['res_weight3'],
        'res_weight4': genref_config['res_weight4'],
        'res_weight5': genref_config['res_weight5'],
        'output_mode': genref_config['output_mode']
    }

    residue_obj_list = [{} for i in range(len(ligands))]
    ref_num = 0
    for protein, ligand in zip(proteins, ligands):
        molformat = protein.split('.')[-1]

        convert = OBConversion()
        convert.SetInFormat(molformat)

        proteinmol = OBMol()
        ligandmol = OBMol()

        convert.ReadFile(proteinmol, protein)
        proteinmol.DeleteNonPolarHydrogens()

        convert.ReadFile(ligandmol, ligand)
        ligandmol.DeleteNonPolarHydrogens()
        if molformat == 'mol2':
            ligand_atom_group = assign_atoms(ligandmol, 'plants')
        if molformat == 'pdbqt':
            ligand_atom_group = assign_atoms(ligandmol, 'vina')

        for name, num in zip(res_name, res_num):
            residue_obj_list[ref_num][name] = Residue(proteinmol, name, num, custom_settings)
            residue_obj_list[ref_num][name].calculateRef(ligandmol, ligand_atom_group)

        ref_num += 1

    return residue_obj_list


def assign_atoms(ligand, docking_method):
    hydrophobic = []
    h_donor = []
    h_accept = []
    positive = []
    negative = []
    h_donorh = []

    rings = getRing(ligand)

    for atom in OBMolAtomIter(ligand):
        # Hydrophobic
        # Hydrophobic atoms SMARTS pattern retrieved from
        # http://www.molsoft.com/icm/smiles.html
        if atom.MatchesSMARTS('[C&!$(C=O)&!$(C#N),S&^3,#17,#15,#35,#53]'):
            hydrophobic.append(atom.GetId())
        # H Donor & Hydrogen bonded to H Donor
        if atom.IsHbondDonor():
            h_donor.append(atom.GetId())
            h_donorh_list = []
            for hydrogen in OBMolAtomIter(ligand):
                if atom.IsConnected(hydrogen):
                    h_donorh_list.append(hydrogen.GetId())
            h_donorh.append(h_donorh_list)
        # H Acceptor
        if atom.IsHbondAcceptor():
            h_accept.append(atom.GetId())
        # Electrostatic
        if docking_method == 'plants':
            if atom.GetPartialCharge() > 0:
                positive.append(atom.GetId())
            if atom.GetPartialCharge() < 0:
                negative.append(atom.GetId())
        if docking_method == 'vina':
            if atom.IsNitrogen() & (atom.GetPartialCharge() >= -0.235):
                positive.append(atom.GetId())
            if atom.IsOxygen() & (atom.GetPartialCharge() <= -0.648):
                negative.append(atom.GetId())

    # Create interaction matrix, similar to AAInteractionMatrix
    interactions = []

    interactions.append(1) if hydrophobic else interactions.append(0)
    interactions.append(1) if rings else interactions.append(0)
    interactions.append(1) if h_accept else interactions.append(0)
    interactions.append(1) if h_donor else interactions.append(0)
    interactions.append(1) if negative else interactions.append(0)
    interactions.append(1) if positive else interactions.append(0)

    atom_group = {
        'hydrophobic': hydrophobic,
        'h_donor': h_donor,
        'h_accept': h_accept,
        'positive': positive,
        'negative': negative,
        'rings': rings,
        'h_donorh': h_donorh,
        'interactions': interactions,
    }

    return atom_group


'''
    function getRing, used in assign_atoms
    return list of ring from GetSSSR
'''


def getRing(mol):
    ringList = []
    for ring in mol.GetSSSR():
        if ring.IsAromatic:
            path = ring._path
            ringList.append(path)

    return ringList


def getCrossModulus(path3):
    path3Coords = []
    for i in path3:
        x = i.x()
        y = i.y()
        z = i.z()
        npcoord = np.array([x, y, z])
        path3Coords.append(npcoord)

    a = path3Coords[0] - path3Coords[1]
    b = path3Coords[0] - path3Coords[2]
    crossProd = np.cross(a, b)
    crossModulus = np.sqrt((crossProd * crossProd).sum())
    return (crossProd, crossModulus)


"""
Class Residue

To contain its own atom groups and bitstring for every pose
"""


class Residue:
    bs_template = {
        'ALA': bitarray('0'),
        'CYS': bitarray('00'),
        'ASP': bitarray('000'),
        'GLU': bitarray('000'),
        'PHE': bitarray('000'),
        'GLY': bitarray('0'),
        'HIS': bitarray('000000'),
        'ILE': bitarray('0'),
        'LYS': bitarray('000'),
        'LEU': bitarray('0'),
        'MET': bitarray('0'),
        'ASN': bitarray('000'),
        'PRO': bitarray('0'),
        'GLN': bitarray('000'),
        'ARG': bitarray('000'),
        'SER': bitarray('00'),
        'THR': bitarray('000'),
        'VAL': bitarray('0'),
        'TRP': bitarray('0000'),
        'TYR': bitarray('00000'),
    }

    hydrophobic = {
        'ALA': {
            'id': [4],
            'bitarray': bitarray('1')
        },
        'CYS': {
            'id': [4],
            'bitarray': bitarray('10')
        },
        'ASP': {
            'id': [4],
            'bitarray': bitarray('100')
        },
        'GLU': {
            'id': [4, 5],
            'bitarray': bitarray('100')
        },
        'PHE': {
            'id': [4, 5, 6, 7, 8, 9, 10],
            'bitarray': bitarray('1000')
        },
        'GLY': {
            'id': [],
            'bitarray': bitarray('0')
        },
        'HIS': {
            'id': [4, 5, 7, 9],
            'bitarray': bitarray('100000')
        },
        'ILE': {
            'id': [4, 5, 6, 7],
            'bitarray': bitarray('1')
        },
        'LYS': {
            'id': [4, 5, 6],
            'bitarray': bitarray('100')
        },
        'LEU': {
            'id': [4, 5, 6, 7],
            'bitarray': bitarray('1')
        },
        'MET': {
            'id': [4, 5, 7],
            'bitarray': bitarray('1')
        },
        'ASN': {
            'id': [4],
            'bitarray': bitarray('100')
        },
        'PRO': {
            'id': [4, 5, 6],
            'bitarray': bitarray('1')
        },
        'GLN': {
            'id': [4, 5],
            'bitarray': bitarray('100')
        },
        'ARG': {
            'id': [4, 5],
            'bitarray': bitarray('100')
        },
        'THR': {
            'id': [6],
            'bitarray': bitarray('100')
        },
        'VAL': {
            'id': [4, 5, 6],
            'bitarray': bitarray('1')
        },
        'TRP': {
            'id': [4, 5, 6, 7, 9, 10, 11, 12, 13],
            'bitarray': bitarray('1000')
        },
        'TYR': {
            'id': [4, 5, 6, 7, 8, 9, 10],
            'bitarray': bitarray('10000')
        },
    }

    aromatic = {
        'PHE': {
            'id': [5, 6, 7, 8, 9, 10],
            'bitarrayF2F': bitarray('010'),
            'bitarrayE2F': bitarray('001')
        },
        'HIS': {
            'id': [5, 6, 7, 8, 9],
            'bitarrayF2F': bitarray('010000'),
            'bitarrayE2F': bitarray('001000')
        },
        'TRP': {
            'id': [5, 6, 7, 8, 9, 10, 11, 12, 13],
            'bitarrayF2F': bitarray('0100'),
            'bitarrayE2F': bitarray('0010')
        },
        'TYR': {
            'id': [5, 6, 7, 8, 9, 10],
            'bitarrayF2F': bitarray('01000'),
            'bitarrayE2F': bitarray('00100')
        },
    }

    h_donor = {
        'CYS': {
            'id': [5],
            'bitarray': bitarray('01')
        },
        'HIS': {
            'id': [6],
            'bitarray': bitarray('000100')
        },
        'LYS': {
            'id': [8],
            'bitarray': bitarray('010')
        },
        'ASN': {
            'id': [7],
            'bitarray': bitarray('010')
        },
        'GLN': {
            'id': [8],
            'bitarray': bitarray('010')
        },
        'ARG': {
            'id': [7, 9, 10],
            'bitarray': bitarray('010')
        },
        'SER': {
            'id': [5],
            'bitarray': bitarray('10')
        },
        'THR': {
            'id': [5],
            'bitarray': bitarray('010')
        },
        'TRP': {
            'id': [8],
            'bitarray': bitarray('0001')
        },
        'TYR': {
            'id': [11],
            'bitarray': bitarray('00010')
        },
    }

    h_donorh = {
        'CYS': {
            'id': [[1]]
        },
        'HIS': {
            'id': [[1]]
        },
        'LYS': {
            'id': [[1, 2, 3]]
        },
        'ASN': {
            'id': [[1, 2]]
        },
        'GLN': {
            'id': [[1, 2]]
        },
        'ARG': {
            'id': [[1], [2, 3], [4, 5]]
        },
        'SER': {
            'id': [[1]]
        },
        'THR': {
            'id': [[1]]
        },
        'TRP': {
            'id': [[1]]
        },
        'TYR': {
            'id': [[1]]
        },
    }

    h_accept = {
        'ASP': {
            'id': [6, 7],
            'bitarray': bitarray('010')
        },
        'GLU': {
            'id': [7, 8],
            'bitarray': bitarray('010')
        },
        'HIS': {
            'id': [9],
            'bitarray': bitarray('000010')
        },
        'ASN': {
            'id': [6],
            'bitarray': bitarray('001')
        },
        'GLN': {
            'id': [7],
            'bitarray': bitarray('001')
        },
        'SER': {
            'id': [5],
            'bitarray': bitarray('01')
        },
        'THR': {
            'id': [5],
            'bitarray': bitarray('001')
        },
        'TYR': {
            'id': [11],
            'bitarray': bitarray('00001')
        },
    }

    positive = {
        'LYS': {
            'id': [8],
            'bitarray': bitarray('001')
        },
        'ARG': {
            'id': [7, 9, 10],
            'bitarray': bitarray('001')
        },
        'HIS': {
            'id': [],
            'bitarray': bitarray('000001')
        }
    }

    negative = {
        'ASP': {
            'id': [6, 7],
            'bitarray': bitarray('001')
        },
        'GLU': {
            'id': [7, 8],
            'bitarray': bitarray('001')
        },
    }

    interactionList = [
        hydrophobic,
        aromatic,
        h_donor,
        h_accept,
        positive,
        negative,
    ]

    interactionNames = [
        'hydrophobic',
        'aromatic',
        'h_donor',
        'h_accept',
        'positive',
        'negative',
    ]

    AAinteractionMatrix = {
        'ALA': (1, 0, 0, 0, 0, 0),
        'CYS': (1, 0, 1, 0, 0, 0),
        'ASP': (1, 0, 0, 1, 0, 1),
        'GLU': (1, 0, 0, 1, 0, 1),
        'PHE': (1, 1, 0, 0, 0, 0),
        'GLY': (0, 0, 0, 0, 0, 0),
        'HIS': (1, 1, 1, 1, 1, 0),
        'ILE': (1, 0, 0, 0, 0, 0),
        'LYS': (1, 0, 1, 0, 1, 0),
        'LEU': (1, 0, 0, 0, 0, 0),
        'MET': (1, 0, 0, 0, 0, 0),
        'ASN': (1, 0, 1, 1, 0, 0),
        'PRO': (1, 0, 0, 0, 0, 0),
        'GLN': (1, 0, 1, 1, 0, 0),
        'ARG': (1, 0, 1, 0, 1, 0),
        'SER': (0, 0, 1, 1, 0, 0),
        'THR': (1, 0, 1, 1, 0, 0),
        'VAL': (1, 0, 0, 0, 0, 0),
        'TRP': (1, 1, 1, 0, 0, 0),
        'TYR': (1, 1, 1, 1, 0, 0),
    }

    def __init__(self, protein, res_name, res_num, custom_settings):
        self.residue_number = int(res_num) - 1
        self.residue = protein.GetResidue(self.residue_number)

        self.atoms = [atom for atom in OBResidueAtomIter(self.residue)]
        # Making atom index consistent by separating hydrogen and heavy atom
        self.heavyatoms = []
        self.hydrogens = []
        for atom in self.atoms:
            if atom.IsHydrogen():
                self.hydrogens.append(atom)
            else:
                self.heavyatoms.append(atom)

        self.output_mode = custom_settings['output_mode']
        self.simplified = False
        self.full = False
        self.full_nobb = False
        if self.output_mode['full']:
            self.full = True
        if self.output_mode['full_nobb']:
            self.full_nobb = True
        if self.output_mode['simplified']:
            self.simplified = True

        self.full_bitstring = bitarray('0000000')
        self.full_nobb_bitstring = bitarray('0000000')
        self.simp_bitstring = self.bs_template[res_name[:3]]

        self.res_name = res_name
        self.AA_name = res_name[:3]
        self.res_num = res_num
        self.interactions = self.AAinteractionMatrix[self.AA_name]

        # Classifying atoms into atom groups
        self.atomGroup = {}
        for AA_int, int_atom, int_name in \
        zip(self.interactions, self.interactionList, self.interactionNames):
            if AA_int == 1:
                atomIDs = int_atom[self.AA_name]['id']
                atomList = [self.heavyatoms[atomID] for atomID in atomIDs]
                self.atomGroup[int_name] = atomList
                if int_name == 'h_donor':
                    h_donorh_list = self.h_donorh[self.AA_name]['id']
                    h_donorh_atoms = []
                    for h_list in h_donorh_list:
                        hydrogen = [self.hydrogens[atomID] for atomID in h_list]
                        h_donorh_atoms.append(hydrogen)
                    self.atomGroup['h_donorh'] = h_donorh_atoms

        if self.AA_name in self.aromatic.keys():
            self.path3 = self.aromatic[self.AA_name]['id'][:3]
            self.path3atoms = []
            for i in self.path3:
                self.path3atoms.append(self.heavyatoms[i])
            self.cross, self.modulus = getCrossModulus(self.path3atoms)

        self.res_weight1 = custom_settings['res_weight1']
        self.res_weight2 = custom_settings['res_weight2']
        self.res_weight3 = custom_settings['res_weight3']
        self.res_weight4 = custom_settings['res_weight4']
        self.res_weight5 = custom_settings['res_weight5']

    def calculateIFPVina(self, ligands, ligand_atom_group):
        pose_num = len(ligands)
        self.full_bits_list = [self.full_bitstring.copy() for i in range(pose_num)]
        self.simp_bits_list = [self.simp_bitstring.copy() for i in range(pose_num)]
        self.full_nobb_list = [self.full_nobb_bitstring.copy() for i in range(pose_num)]

        possible_interactions = []
        for x, y in zip(self.interactions, ligand_atom_group['interactions']):
            if x & y:
                possible_interactions.append(1)
            else:
                possible_interactions.append(0)

        for ligand, full, full_nobb, simp in \
        zip(ligands, self.full_bits_list, self.full_nobb_list, self.simp_bits_list):
            interaction_flags = [0, 0, 0, 0, 0, 0, 0]

            # hydrophobic
            if possible_interactions[0]:
                for ligand_id in ligand_atom_group['hydrophobic']:
                    ligand_atom = ligand.GetAtomById(ligand_id)
                    for atom in self.atomGroup['hydrophobic']:
                        distance = ligand_atom.GetDistance(atom)
                        if distance <= HYDROPHOBIC:
                            if self.simplified:
                                simp |= self.hydrophobic[self.AA_name]['bitarray']
                            if self.full:
                                full |= bitarray('1000000')
                            if self.full_nobb:
                                full_nobb |= bitarray('1000000')
                            interaction_flags[0] = 1
                            break
                    if interaction_flags[0]:
                        break

            # aromatic EF & FF
            if possible_interactions[1]:
                for path in ligand_atom_group['rings']:
                    path3 = path[:3]
                    path3atoms = []
                    for i in path3:
                        path3atoms.append(ligand.GetAtom(i))
                    for ligand_id in path:
                        for atom in self.atomGroup['aromatic']:
                            distance = ligand.GetAtom(ligand_id).GetDistance(atom)
                            if distance <= AROMATIC:
                                ligCross, ligModulus = getCrossModulus(path3atoms)

                                dot = np.dot(self.cross, ligCross)
                                cos_angle = dot / self.modulus / ligModulus
                                ring_angle = np.arccos(cos_angle) * 180 / np.pi

                                if (AROMATIC_ANGLE_LOW >= ring_angle) or (AROMATIC_ANGLE_HIGH <= ring_angle):
                                    if self.simplified:
                                        simp |= self.aromatic[self.AA_name]['bitarrayF2F']
                                    if self.full:
                                        full |= bitarray('0100000')
                                    if self.full_nobb:
                                        full_nobb |= bitarray('0100000')
                                    interaction_flags[1] = 1
                                else:
                                    if self.simplified:
                                        simp |= self.aromatic[self.AA_name]['bitarrayE2F']
                                    if self.full:
                                        full |= bitarray('0010000')
                                    if self.full_nobb:
                                        full_nobb |= bitarray('0010000')
                                    interaction_flags[2] = 1

                                break
                        if interaction_flags[1] | interaction_flags[2]:
                            break
                    if interaction_flags[1] & interaction_flags[2]:
                        break

            # H acceptor
            if possible_interactions[2]:
                for ligand_id in ligand_atom_group['h_accept']:
                    ligand_atom = ligand.GetAtomById(ligand_id)
                    for atom, h_list in zip(self.atomGroup['h_donor'], self.atomGroup['h_donorh']):
                        distance = ligand_atom.GetDistance(atom)
                        if distance <= HBOND:
                            angle_flag = 0
                            for hydrogen in h_list:
                                angle = atom.GetAngle(hydrogen, ligand_atom)
                                if angle > HBOND_ANGLE:
                                    angle_flag = 1
                                    break

                            if angle_flag:
                                if self.simplified:
                                    simp |= self.h_donor[self.AA_name]['bitarray']
                                if self.full:
                                    full |= bitarray('0001000')
                                if self.full_nobb:
                                    full_nobb |= bitarray('0001000')
                                interaction_flags[3] = 1
                                break
                    if interaction_flags[3]:
                        break

            # H donor
            if possible_interactions[3]:
                for ligand_id, h_list in zip(ligand_atom_group['h_donor'], ligand_atom_group['h_donorh']):
                    ligand_atom = ligand.GetAtomById(ligand_id)
                    for atom in self.atomGroup['h_accept']:
                        distance = ligand_atom.GetDistance(atom)
                        if distance <= HBOND:
                            angle_flag = 0
                            for h in h_list:
                                hydrogen = ligand.GetAtomById(h)
                                angle = ligand_atom.GetAngle(hydrogen, atom)
                                if angle > HBOND_ANGLE:
                                    angle_flag = 1
                                    break
                            if angle_flag:
                                if self.simplified:
                                    simp |= self.h_accept[self.AA_name]['bitarray']
                                if self.full:
                                    full |= bitarray('0000100')
                                if self.full_nobb:
                                    full_nobb |= bitarray('0000100')
                                interaction_flags[4] = 1
                        if interaction_flags[4]:
                            break

            # ligand negative - protein positive
            if possible_interactions[4]:
                for ligand_id in ligand_atom_group['negative']:
                    ligand_atom = ligand.GetAtomById(ligand_id)
                    for atom in self.atomGroup['positive']:
                        distance = ligand_atom.GetDistance(atom)
                        if distance <= ELECTROSTATIC:
                            if self.simplified:
                                simp |= self.positive[self.AA_name]['bitarray']
                            if self.full:
                                full |= bitarray('0000010')
                            if self.full_nobb:
                                full_nobb |= bitarray('0000010')
                            interaction_flags[5] = 1
                            break
                    if interaction_flags[5]:
                        break

            # ligand positive - protein negative
            if possible_interactions[5]:
                for ligand_id in ligand_atom_group['positive']:
                    ligand_atom = ligand.GetAtomById(ligand_id)
                    for atom in self.atomGroup['negative']:
                        distance = ligand_atom.GetDistance(atom)
                        if distance <= ELECTROSTATIC:
                            if self.simplified:
                                simp |= self.negative[self.AA_name]['bitarray']
                            if self.full:
                                full |= bitarray('0000001')
                            if self.full_nobb:
                                full_nobb |= bitarray('0000001')
                            interaction_flags[6] = 1
                            break
                    if interaction_flags[6]:
                        break

            if self.full:
                if ligand_atom_group['interactions'][0]:
                    for ligand_id in ligand_atom_group['hydrophobic']:
                        ligand_atom = ligand.GetAtomById(ligand_id)
                        bb_atom = self.heavyatoms[1]
                        distance = ligand_atom.GetDistance(bb_atom)
                        if distance <= HYDROPHOBIC:
                            full |= bitarray('1000000')
                            break
                notPRO = False if self.AA_name == 'PRO' else True
                if ligand_atom_group['interactions'][3] & notPRO:
                    for ligand_id in ligand_atom_group['h_accept']:
                        bb_atom = self.heavyatoms[0]
                        ligand_atom = ligand.GetAtomById(ligand_id)
                        distance = ligand_atom.GetDistance(bb_atom)
                        if distance <= HBOND:
                            hydrogen = self.hydrogens[0]
                            angle = bb_atom.GetAngle(hydrogen, ligand_atom)
                            if angle > HBOND_ANGLE:
                                full |= bitarray('0001000')
                                break
                if ligand_atom_group['interactions'][2]:
                    for ligand_id, h_list in zip(ligand_atom_group['h_donor'], ligand_atom_group['h_donorh']):
                        ligand_atom = ligand.GetAtomById(ligand_id)
                        bb_atom = self.heavyatoms[3]
                        distance = ligand_atom.GetDistance(bb_atom)
                        if distance <= HBOND:
                            angle_flag = 0
                            for h in h_list:
                                hydrogen = ligand.GetAtomById(h)
                                angle = ligand_atom.GetAngle(hydrogen, bb_atom)
                                if angle > HBOND_ANGLE:
                                    angle_flag = 1
                                    break
                            if angle_flag:
                                full |= bitarray('0000100')
                                break

    def calculateIFPPlants(self, ligands, flex_proteins, ligand_atom_group):
        pose_num = len(ligands)
        self.flex_residues = [residue.GetName() for residue in OBResidueIter(flex_proteins[0])]

        self.full_bits_list = [self.full_bitstring.copy() for i in range(pose_num)]
        self.simp_bits_list = [self.simp_bitstring.copy() for i in range(pose_num)]
        self.full_nobb_list = [self.full_nobb_bitstring.copy() for i in range(pose_num)]

        possible_interactions = []
        for x, y in zip(self.interactions, ligand_atom_group['interactions']):
            possible_interactions.append(1) if x & y else possible_interactions.append(0)
        for ligand, flex, full, full_nobb, simp in \
        zip(ligands, flex_proteins, self.full_bits_list, self.full_nobb_list, self.simp_bits_list):
            interaction_flags = [0, 0, 0, 0, 0, 0, 0]

            if possible_interactions[0]:
                for ligand_id in ligand_atom_group['hydrophobic']:
                    ligand_atom = ligand.GetAtomById(ligand_id)
                    for atom in self.atomGroup['hydrophobic']:
                        distance = ligand_atom.GetDistance(atom)
                        if distance <= HYDROPHOBIC:
                            if self.simplified:
                                simp |= self.hydrophobic[self.AA_name]['bitarray']
                            if self.full:
                                full |= bitarray('1000000')
                            if self.full_nobb:
                                full_nobb |= bitarray('1000000')
                            interaction_flags[0] = 1
                            break
                    if interaction_flags[0]:
                        break

            if possible_interactions[1]:
                for path in ligand_atom_group['rings']:
                    path3 = path[:3]
                    path3atoms = []
                    for i in path3:
                        path3atoms.append(ligand.GetAtom(i))
                    for ligand_id in path:
                        for atom in self.atomGroup['aromatic']:
                            distance = ligand.GetAtom(ligand_id).GetDistance(atom)
                            if distance <= AROMATIC:
                                ligCross, ligModulus = getCrossModulus(path3atoms)

                                dot = np.dot(self.cross, ligCross)
                                cos_angle = dot / self.modulus / ligModulus
                                ring_angle = np.arccos(cos_angle) * 180 / np.pi

                                if (AROMATIC_ANGLE_LOW >= ring_angle) or (AROMATIC_ANGLE_HIGH <= ring_angle):
                                    if self.simplified:
                                        simp |= self.aromatic[self.AA_name]['bitarrayF2F']
                                    if self.full:
                                        full |= bitarray('0100000')
                                    if self.full_nobb:
                                        full_nobb |= bitarray('0100000')
                                    interaction_flags[1] = 1
                                else:
                                    if self.simplified:
                                        simp |= self.aromatic[self.AA_name]['bitarrayE2F']
                                    if self.full:
                                        full |= bitarray('0010000')
                                    if self.full_nobb:
                                        full_nobb |= bitarray('0010000')
                                    interaction_flags[2] = 1

                                break
                        if interaction_flags[1] | interaction_flags[2]:
                            break
                    if interaction_flags[1] & interaction_flags[2]:
                        break

            if possible_interactions[2]:
                for ligand_id in ligand_atom_group['h_accept']:
                    ligand_atom = ligand.GetAtomById(ligand_id)
                    for atom, h_list in zip(self.atomGroup['h_donor'], self.atomGroup['h_donorh']):
                        distance = ligand_atom.GetDistance(atom)
                        if distance <= HBOND:
                            angle_flag = 0
                            if self.res_name in self.flex_residues:
                                for flex_atom in OBMolAtomIter(flex):
                                    if flex_atom.IsHydrogen() & (flex_atom.GetResidue().GetName() == self.res_name):
                                        angle = atom.GetAngle(flex_atom, ligand_atom)
                                        if angle > HBOND_ANGLE:
                                            angle_flag = 1
                                            break

                            else:
                                for hydrogen in h_list:
                                    angle = atom.GetAngle(hydrogen, ligand_atom)
                                    if angle > HBOND_ANGLE:
                                        angle_flag = 1
                                        break

                            if angle_flag:
                                if self.simplified:
                                    simp |= self.h_donor[self.AA_name]['bitarray']
                                if self.full:
                                    full |= bitarray('0001000')
                                if self.full_nobb:
                                    full_nobb |= bitarray('0001000')
                                interaction_flags[3] = 1
                                break
                    if interaction_flags[3]:
                        break

            if possible_interactions[3]:
                for ligand_id, h_list in zip(ligand_atom_group['h_donor'], ligand_atom_group['h_donorh']):
                    ligand_atom = ligand.GetAtomById(ligand_id)
                    for atom in self.atomGroup['h_accept']:
                        distance = ligand_atom.GetDistance(atom)
                        if distance <= HBOND:
                            angle_flag = 0
                            for h in h_list:
                                hydrogen = ligand.GetAtomById(h)
                                angle = ligand_atom.GetAngle(hydrogen, atom)
                                if angle > HBOND_ANGLE:
                                    angle_flag = 1
                                    break
                            if angle_flag:
                                if self.simplified:
                                    simp |= self.h_accept[self.AA_name]['bitarray']
                                if self.full:
                                    full |= bitarray('0000100')
                                if self.full_nobb:
                                    full_nobb |= bitarray('0000100')
                                interaction_flags[4] = 1
                        if interaction_flags[4]:
                            break

            if possible_interactions[4]:
                for ligand_id in ligand_atom_group['negative']:
                    ligand_atom = ligand.GetAtomById(ligand_id)
                    for atom in self.atomGroup['positive']:
                        distance = ligand_atom.GetDistance(atom)
                        if distance <= ELECTROSTATIC:
                            if self.simplified:
                                simp |= self.positive[self.AA_name]['bitarray']
                            if self.full:
                                full |= bitarray('0000010')
                            if self.full_nobb:
                                full_nobb |= bitarray('0000010')
                            interaction_flags[5] = 1
                            break
                    if interaction_flags[5]:
                        break

            if possible_interactions[5]:
                for ligand_id in ligand_atom_group['positive']:
                    ligand_atom = ligand.GetAtomById(ligand_id)
                    for atom in self.atomGroup['negative']:
                        distance = ligand_atom.GetDistance(atom)
                        if distance <= ELECTROSTATIC:
                            if self.simplified:
                                simp |= self.negative[self.AA_name]['bitarray']
                            if self.full:
                                full |= bitarray('0000001')
                            if self.full_nobb:
                                full_nobb |= bitarray('0000001')
                            interaction_flags[6] = 1
                            break
                    if interaction_flags[6]:
                        break

            if self.full:
                if ligand_atom_group['interactions'][0]:
                    for ligand_id in ligand_atom_group['hydrophobic']:
                        ligand_atom = ligand.GetAtomById(ligand_id)
                        bb_atom = self.heavyatoms[1]
                        distance = ligand_atom.GetDistance(bb_atom)
                        if distance <= HYDROPHOBIC:
                            full |= bitarray('1000000')
                            break
                notPRO = False if self.AA_name == 'PRO' else True
                if ligand_atom_group['interactions'][3] & notPRO:
                    for ligand_id in ligand_atom_group['h_accept']:
                        bb_atom = self.heavyatoms[0]
                        ligand_atom = ligand.GetAtomById(ligand_id)
                        distance = ligand_atom.GetDistance(bb_atom)
                        if distance <= HBOND:
                            hydrogen = self.hydrogens[0]
                            angle = bb_atom.GetAngle(hydrogen, ligand_atom)
                            if angle > HBOND_ANGLE:
                                full |= bitarray('0001000')
                                break
                if ligand_atom_group['interactions'][2]:
                    for ligand_id, h_list in zip(ligand_atom_group['h_donor'], ligand_atom_group['h_donorh']):
                        ligand_atom = ligand.GetAtomById(ligand_id)
                        bb_atom = self.heavyatoms[3]
                        distance = ligand_atom.GetDistance(bb_atom)
                        if distance <= HBOND:
                            angle_flag = 0
                            for h in h_list:
                                hydrogen = ligand.GetAtomById(h)
                                angle = ligand_atom.GetAngle(hydrogen, bb_atom)
                                if angle > HBOND_ANGLE:
                                    angle_flag = 1
                                    break
                            if angle_flag:
                                full |= bitarray('0000100')
                                break

    def calculateRef(self, ligand, ligand_atom_group):
        possible_interactions = []
        for x, y in zip(self.interactions, ligand_atom_group['interactions']):
            possible_interactions.append(1) if x & y else possible_interactions.append(0)

        self.full_bit = self.full_bitstring.copy()
        self.simp = self.simp_bitstring.copy()
        self.nobb_bit = self.full_nobb_bitstring.copy()
        interaction_flags = [0, 0, 0, 0, 0, 0, 0]

        if possible_interactions[0]:
            for ligand_id in ligand_atom_group['hydrophobic']:
                ligand_atom = ligand.GetAtomById(ligand_id)
                for atom in self.atomGroup['hydrophobic']:
                    distance = ligand_atom.GetDistance(atom)
                    if distance <= HYDROPHOBIC:
                        if self.simplified:
                            self.simp |= self.hydrophobic[self.AA_name]['bitarray']
                        if self.full:
                            self.full_bit |= bitarray('1000000')
                        if self.full_nobb:
                            self.nobb_bit |= bitarray('1000000')
                        interaction_flags[0] = 1
                        break
                if interaction_flags[0]:
                    break

        if possible_interactions[1]:
            for path in ligand_atom_group['rings']:
                path3 = path[:3]
                path3atoms = []
                for i in path3:
                    path3atoms.append(ligand.GetAtom(i))
                for ligand_id in path:
                    for atom in self.atomGroup['aromatic']:
                        distance = ligand.GetAtom(ligand_id).GetDistance(atom)
                        if distance <= AROMATIC:
                            ligCross, ligModulus = getCrossModulus(path3atoms)

                            dot = np.dot(self.cross, ligCross)
                            cos_angle = dot / self.modulus / ligModulus
                            ring_angle = np.arccos(cos_angle) * 180 / np.pi

                            if (AROMATIC_ANGLE_LOW >= ring_angle) or (AROMATIC_ANGLE_HIGH <= ring_angle):
                                if self.simplified:
                                    self.simp |= self.aromatic[self.AA_name]['bitarrayF2F']
                                if self.full:
                                    self.full_bit |= bitarray('0100000')
                                if self.full_nobb:
                                    self.nobb_bit |= bitarray('0100000')
                                interaction_flags[1] = 1
                            else:
                                if self.simplified:
                                    self.simp |= self.aromatic[self.AA_name]['bitarrayE2F']
                                if self.full:
                                    self.full_bit |= bitarray('0010000')
                                if self.full_nobb:
                                    self.nobb_bit |= bitarray('0010000')
                                interaction_flags[2] = 1

                            break
                    if interaction_flags[1] | interaction_flags[2]:
                        break
                if interaction_flags[1] & interaction_flags[2]:
                    break

        if possible_interactions[2]:
            for ligand_id in ligand_atom_group['h_accept']:
                ligand_atom = ligand.GetAtomById(ligand_id)
                for atom, h_list in zip(self.atomGroup['h_donor'], self.atomGroup['h_donorh']):
                    distance = ligand_atom.GetDistance(atom)
                    if distance <= HBOND:
                        angle_flag = 0
                        for hydrogen in h_list:
                            angle = atom.GetAngle(hydrogen, ligand_atom)
                            if angle > HBOND_ANGLE:
                                angle_flag = 1
                                break

                        if angle_flag:
                            if self.simplified:
                                self.simp |= self.h_donor[self.AA_name]['bitarray']
                            if self.full:
                                self.full_bit |= bitarray('0001000')
                            if self.full_nobb:
                                self.nobb_bit |= bitarray('0001000')
                            interaction_flags[3] = 1
                            break
                if interaction_flags[3]:
                    break

        if possible_interactions[3]:
            for ligand_id, h_list in zip(ligand_atom_group['h_donor'], ligand_atom_group['h_donorh']):
                ligand_atom = ligand.GetAtomById(ligand_id)
                for atom in self.atomGroup['h_accept']:
                    distance = ligand_atom.GetDistance(atom)
                    if distance <= HBOND:
                        angle_flag = 0
                        for h in h_list:
                            hydrogen = ligand.GetAtomById(h)
                            angle = ligand_atom.GetAngle(hydrogen, atom)
                            if angle > HBOND_ANGLE:
                                angle_flag = 1
                                break
                        if angle_flag:
                            if self.simplified:
                                self.simp |= self.h_accept[self.AA_name]['bitarray']
                            if self.full:
                                self.full_bit |= bitarray('0000100')
                            if self.full_nobb:
                                self.nobb_bit |= bitarray('0000100')
                            interaction_flags[4] = 1
                    if interaction_flags[4]:
                        break

        if possible_interactions[4]:
            for ligand_id in ligand_atom_group['negative']:
                ligand_atom = ligand.GetAtomById(ligand_id)
                for atom in self.atomGroup['positive']:
                    distance = ligand_atom.GetDistance(atom)
                    if distance <= ELECTROSTATIC:
                        if self.simplified:
                            self.simp |= self.positive[self.AA_name]['bitarray']
                        if self.full:
                            self.full_bit |= bitarray('0000010')
                        if self.full_nobb:
                            self.nobb_bit |= bitarray('0000010')
                        interaction_flags[5] = 1
                        break
                if interaction_flags[5]:
                    break

        if possible_interactions[5]:
            for ligand_id in ligand_atom_group['positive']:
                ligand_atom = ligand.GetAtomById(ligand_id)
                for atom in self.atomGroup['negative']:
                    distance = ligand_atom.GetDistance(atom)
                    if distance <= ELECTROSTATIC:
                        if self.simplified:
                            self.simp |= self.negative[self.AA_name]['bitarray']
                        if self.full:
                            self.full_bit |= bitarray('0000001')
                        if self.full_nobb:
                            self.nobb_bit |= bitarray('0000001')
                        interaction_flags[6] = 1
                        break
                if interaction_flags[6]:
                    break

        if self.full:
            if ligand_atom_group['interactions'][0]:
                for ligand_id in ligand_atom_group['hydrophobic']:
                    ligand_atom = ligand.GetAtomById(ligand_id)
                    bb_atom = self.heavyatoms[1]
                    distance = ligand_atom.GetDistance(bb_atom)
                    if distance <= HYDROPHOBIC:
                        self.full_bit |= bitarray('1000000')
                        break
            notPRO = False if self.AA_name == 'PRO' else True
            if ligand_atom_group['interactions'][3] & notPRO:
                for ligand_id in ligand_atom_group['h_accept']:
                    bb_atom = self.heavyatoms[0]
                    ligand_atom = ligand.GetAtomById(ligand_id)
                    distance = ligand_atom.GetDistance(bb_atom)
                    if distance <= HBOND:
                        hydrogen = self.hydrogens[0]
                        angle = bb_atom.GetAngle(hydrogen, ligand_atom)
                        if angle > HBOND_ANGLE:
                            self.full_bit |= bitarray('0001000')
                            break
            if ligand_atom_group['interactions'][2]:
                for ligand_id, h_list in zip(ligand_atom_group['h_donor'], ligand_atom_group['h_donorh']):
                    ligand_atom = ligand.GetAtomById(ligand_id)
                    bb_atom = self.heavyatoms[3]
                    distance = ligand_atom.GetDistance(bb_atom)
                    if distance <= HBOND:
                        angle_flag = 0
                        for h in h_list:
                            hydrogen = ligand.GetAtomById(h)
                            angle = ligand_atom.GetAngle(hydrogen, bb_atom)
                            if angle > HBOND_ANGLE:
                                angle_flag = 1
                                break
                        if angle_flag:
                            self.full_bit |= bitarray('0000100')
                            break
