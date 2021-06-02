"""
=======================================================================
  Rosetta stone for residue translation
=======================================================================

  Data
  ----
    rosetta_residues
    rosetta_atoms

"""

# aa translations ( www.bmrb.wisc.edu/ref_info/atom_nom.tbl and GROMACS documentation $GMXPREFIX/share/gromacs/top/atom_nom.tbl )
rosetta_residues = {
    'amber': ['HIS', 'HID', 'HIE', 'HIP', 'CYS', 'CYM', 'CYX', 'LYS', 'LYN', 'ASP', 'ASH', 'GLU', 'GLH'],
    'charmm': ['HIS', 'HSD', 'HSE', 'HSP', 'CYS', 'CYS', 'CYS', 'LYS', 'LSN', 'ASP', 'ASP', 'GLU', 'GLU'],
    'opls': ['HIS', 'HIS', 'HIS', 'HIS', 'CYS', 'CYS', 'CYS', 'LYS', 'LYS', 'ASP', 'ASP', 'GLU', 'GLU']
    }

rosetta_atoms = {
    'dynamo': {
        'ALA': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'ARG': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'HD1', 'HD2', 'NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'ASN': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'OD1', 'ND2', 'HD21', 'HD22', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'ASP': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'OD1', 'OD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'CYS': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'SG', 'HG', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'GLN': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'OE1', 'NE2', 'HE21', 'HE22', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'GLU': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'OE1', 'OE2', 'HE2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'GLY': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA1', 'HA2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HIS': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'ILE': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23', 'CG1', 'HG11', 'HG12', 'CD', 'HD1', 'HD2', 'HD3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'LEU': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG', 'CD1', 'HD11', 'HD12', 'HD13', 'CD2', 'HD21', 'HD22', 'HD23', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'LYS': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'HD1', 'HD2', 'CE', 'HE1', 'HE2', 'NZ', 'HZ1', 'HZ2', 'HZ3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'MET': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'SD', 'CE', 'HE1', 'HE2', 'HE3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'PHE': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'HZ', 'CE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'PRO': ['N', 'H1', 'H2', 'H3', 'CD', 'HD1', 'HD2', 'CG', 'HG1', 'HG2', 'CB', 'HB1', 'HB2', 'CA', 'HA', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'SER': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'OG', 'HG', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'THR': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23', 'OG1', 'HG1', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'TRP': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'CD1', 'HD1', 'NE1', 'HE1', 'CE2', 'CZ2', 'HZ2', 'CH2', 'HH2', 'CZ3', 'HZ3', 'CE3', 'HE3', 'CD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'TYR': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'OH', 'HH', 'CE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'VAL': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG1', 'HG11', 'HG12', 'HG13', 'CG2', 'HG21', 'HG22', 'HG23', 'C', 'O', 'OC1', 'OC2', 'OXT'],

        'HID': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HIE': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HIP': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HSD': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HSE': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HSP': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'CYX': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'SG', 'C', 'O', 'OC1', 'OC2', 'OXT']
        },
    'gmx': {
        'ALA': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'ARG': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'HD1', 'HD2', 'NE', 'HE', 'CZ', 'NH1', '1HH1', '2HH1', 'NH2', '1HH2', '2HH2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'ASN': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'OD1', 'ND2', '1HD2', '2HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'ASP': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'OD1', 'OD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'CYS': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'SG', 'HG', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'GLN': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'OE1', 'NE2', '1HE2', '2HE2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'GLU': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'OE1', 'OE2', 'HE2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'GLY': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA1', 'HA2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HIS': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'ILE': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG2', '1HG2', '2HG2', '3HG2', 'CG1', '1HG1', '2HG1', 'CD', 'HD1', 'HD2', 'HD3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'LEU': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG', 'CD1', '1HD1', '2HD1', '3HD1', 'CD2', '1HD2', '2HD2', '3HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'LYS': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'HD1', 'HD2', 'CE', 'HE1', 'HE2', 'NZ', 'HZ1', 'HZ2', 'HZ3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'MET': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'SD', 'CE', 'HE1', 'HE2', 'HE3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'PHE': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'HZ', 'CE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'PRO': ['N', 'H1', 'H2', 'H3', 'CD', 'HD1', 'HD2', 'CG', 'HG1', 'HG2', 'CB', 'HB1', 'HB2', 'CA', 'HA', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'SER': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'OG', 'HG', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'THR': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG2', '1HG2', '2HG2', '3HG2', 'OG1', 'HG1', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'TRP': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'CD1', 'HD1', 'NE1', 'HE1', 'CE2', 'CZ2', 'HZ2', 'CH2', 'HH2', 'CZ3', 'HZ3', 'CE3', 'HE3', 'CD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'TYR': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'OH', 'HH', 'CE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'VAL': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG1', '1HG1', '2HG1', '3HG1', 'CG2', '1HG2', '2HG2', '3HG2', 'C', 'O', 'OC1', 'OC2', 'OXT'],

        'HID': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HIE': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HIP': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HSD': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HSE': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HSP': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'CYX': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'SG', 'C', 'O', 'OC1', 'OC2', 'OXT']
        },
    'maestro': {
        'ALA': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'ARG': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'ASN': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'OD1', 'ND2', 'HD21', 'HD22', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'ASP': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'OD1', 'OD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'CYS': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'SG', 'HG', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'GLN': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'OE1', 'NE2', 'HE21', 'HE22', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'GLU': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'OE1', 'OE2', 'HE2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'GLY': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA2', 'HA3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HIS': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'ILE': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23', 'CG1', 'HG12', 'HG13', 'CD1', 'HD11', 'HD12', 'HD13', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'LEU': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG', 'CD1', 'HD11', 'HD12', 'HD13', 'CD2', 'HD21', 'HD22', 'HD23', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'LYS': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'NZ', 'HZ1', 'HZ2', 'HZ3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'MET': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'SD', 'CE', 'HE1', 'HE2', 'HE3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'PHE': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'HZ', 'CE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'PRO': ['N', 'H1', 'H2', 'H3', 'CD', 'HD2', 'HD3', 'CG', 'HG2', 'HG3', 'CB', 'HB2', 'HB3', 'CA', 'HA', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'SER': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'OG', 'HG', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'THR': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23', 'OG1', 'HG1', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'TRP': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'NE1', 'HE1', 'CE2', 'CZ2', 'HZ2', 'CH2', 'HH2', 'CZ3', 'HZ3', 'CE3', 'HE3', 'CD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'TYR': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'OH', 'HH', 'CE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'VAL': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG1', 'HG11', 'HG12', 'HG13', 'CG2', 'HG21', 'HG22', 'HG23', 'C', 'O', 'OC1', 'OC2', 'OXT'],

        'HID': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HIE': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HIP': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HSD': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HSE': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'HSP': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
        'CYX': ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'SG', 'C', 'O', 'OC1', 'OC2', 'OXT']
        }
    }
