#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# File: pdb4all.py
# Description : Protein conversion between common pdb formats and name conventions
# Version : 0.3.4
# Last update : 02-12-2020
# Author : Sergio Boneta

#######################################################################
##                                                                   ##
##                              pdb4all                              ##
##                                                                   ##
#######################################################################
#
#
# Copyright (C) 2020, Sergio Boneta
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see https://www.gnu.org/licenses/
#
#
# Check for updates:  https://github.com/boneta/pdb4all


"""
pdb4all
=======

Protein conversion between common pdb formats and name conventions


Classes
-------

    pdb

Functions
---------

    standard_preparation
"""

#######################################################################
##  DEPENDENCIES                                                     ##
#######################################################################
## Support: Python >2.7 and >3.7

import os
import sys
import math as m
import argparse

__version__ = '0.3.5'

#######################################################################
##  PARSER                                                           ##
#######################################################################
## OpenBabel inspired

def __parserbuilder():
    # supported formats
    input_formats = ['maestro', 'gmx', 'dynamo']
    output_formats = ['fasta', 'gmx', 'dynamo', 'intseq']
    ff_formats = ['amber', 'charmm', 'opls']

    input_formats_print = "\n ".join(input_formats)
    output_formats_print = "\n ".join(output_formats)
    ff_formats_print = "\n ".join(ff_formats)

    # parser building
    parser = argparse.ArgumentParser(prog='pdb4all',
                                    description=' -- Convert between common pdb formats and names --\n\n'+
                                    'input formats:\n '+input_formats_print + "\n\n" +
                                    'output formats:\n '+output_formats_print + "\n\n" +
                                    'FF formats:\n '+ff_formats_print + "\n",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-v',
                        '--version',
                        action='version',
                        version='pdb4all  v{}\nby Sergio Boneta / GPL'.format(__version__))
    parser.add_argument('I',
                        metavar='.pdb',
                        type=str,
                        help='input pdb file')
    parser.add_argument('-i',
                        required=True,
                        metavar='<format>',
                        type=str,
                        choices=input_formats,
                        help='input pdb format (software)')
    parser.add_argument('-o',
                        required=True,
                        metavar='<format>',
                        type=str,
                        choices=output_formats,
                        help='output pdb format (software)')
    parser.add_argument('-O',
                        metavar='.pdb',
                        type=str,
                        help='output pdb file')
    parser.add_argument('-iff',
                        metavar='ff',
                        type=str,
                        choices=ff_formats,
                        default='amber',
                        help='input Force Field format (def: amber)')
    parser.add_argument('-off',
                        metavar='ff',
                        type=str,
                        choices=ff_formats,
                        default='amber',
                        help='output Force Field format (def: amber)')
    parser.add_argument('--dynamize',
                        action='store_true',
                        help='extra fDynamo files (.crd, seq, lig_opls)')
    parser.add_argument('--center',
                        action='store_true',
                        help='move system to center (0,0,0)')
    parser.add_argument('--simple',
                        action='store_false',
                        help='read simplified pdb based on common columns\n'+
                        '  ATOM/HETATM  #Atom  Atom  Res  [Chain]  #Res  X  Y  Z  Occup  TempFac  [Segment]')
    return parser


#######################################################################
##  INITIALIZATION                                                   ##
#######################################################################

## Periodic Table ·····················································
# www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses
# m: mass of most common isotope,  m_std: standard atomic weight
Ptable = {  'H'  : { 'N':   1, 'm':   1.00782503223, 'm_std':   (1.00784 + 1.00811)/2 },
            'He' : { 'N':   2, 'm':   4.00260325413, 'm_std':    4.002602 },
            'Li' : { 'N':   3, 'm':   7.0160034366 , 'm_std':   (6.938 + 6.997)/2 },
            'Be' : { 'N':   4, 'm':   9.012183065  , 'm_std':    9.0121831 },
            'B'  : { 'N':   5, 'm':  11.00930536   , 'm_std':  (10.806 + 10.821)/2 },
            'C'  : { 'N':   6, 'm':  12.00000000000, 'm_std':  (12.0096 + 12.0116)/2 },
            'N'  : { 'N':   7, 'm':  14.00307400443, 'm_std':  (14.00643 + 14.00728)/2 },
            'O'  : { 'N':   8, 'm':  15.99491461957, 'm_std':  (15.99903 + 15.99977)/2 },
            'F'  : { 'N':   9, 'm':  18.99840316273, 'm_std':   18.998403163 },
            'Ne' : { 'N':  10, 'm':  19.9924401762 , 'm_std':   20.1797 },
            'Na' : { 'N':  11, 'm':  22.9897692820 , 'm_std':   22.98976928 },
            'Mg' : { 'N':  12, 'm':  23.985041697  , 'm_std':  (24.304 + 24.307)/2 },
            'Al' : { 'N':  13, 'm':  26.98153853   , 'm_std':   26.9815385 },
            'Si' : { 'N':  14, 'm':  27.97692653465, 'm_std':  (28.084 + 28.086)/2 },
            'P'  : { 'N':  15, 'm':  30.97376199842, 'm_std':   30.973761998 },
            'S'  : { 'N':  16, 'm':  31.9720711744 , 'm_std':  (32.059 + 32.076)/2 },
            'Cl' : { 'N':  17, 'm':  34.968852682  , 'm_std':  (35.446 + 35.457)/2 },
            'Ar' : { 'N':  18, 'm':  39.9623831237 , 'm_std':   39.948 },
            'K'  : { 'N':  19, 'm':  38.9637064864 , 'm_std':   39.0983 },
            'Ca' : { 'N':  20, 'm':  39.962590863  , 'm_std':   40.078 },
            'Sc' : { 'N':  21, 'm':  44.95590828   , 'm_std':   44.955908 },
            'Ti' : { 'N':  22, 'm':  47.94794198   , 'm_std':   47.867 },
            'V'  : { 'N':  23, 'm':  50.94395704   , 'm_std':   50.9415 },
            'Cr' : { 'N':  24, 'm':  51.94050623   , 'm_std':   51.9961 },
            'Mn' : { 'N':  25, 'm':  54.93804391   , 'm_std':   54.938044 },
            'Fe' : { 'N':  26, 'm':  55.93493633   , 'm_std':   55.845 },
            'Co' : { 'N':  27, 'm':  58.93319429   , 'm_std':   58.933194 },
            'Ni' : { 'N':  28, 'm':  57.93534241   , 'm_std':   58.6934 },
            'Cu' : { 'N':  29, 'm':  62.92959772   , 'm_std':   63.546 },
            'Zn' : { 'N':  30, 'm':  63.92914201   , 'm_std':   65.38 },
            'Ga' : { 'N':  31, 'm':  68.9255735    , 'm_std':   69.723 },
            'Ge' : { 'N':  32, 'm':  73.921177761  , 'm_std':   72.630 },
            'As' : { 'N':  33, 'm':  74.92159457   , 'm_std':   74.921595 },
            'Se' : { 'N':  34, 'm':  79.9165218    , 'm_std':   78.971 },
            'Br' : { 'N':  35, 'm':  78.9183376    , 'm_std':  (79.901 + 79.907)/2 },
            'Kr' : { 'N':  36, 'm':  83.9114977282 , 'm_std':   83.798 },
            'Rb' : { 'N':  37, 'm':  84.9117897379 , 'm_std':   85.4678 },
            'Sr' : { 'N':  38, 'm':  87.9056125    , 'm_std':   87.62 },
            'Y'  : { 'N':  39, 'm':  88.9058403    , 'm_std':   88.90584 },
            'Zr' : { 'N':  40, 'm':  89.9046977    , 'm_std':   91.224 },
            'Nb' : { 'N':  41, 'm':  92.9063730    , 'm_std':   92.90637 },
            'Mo' : { 'N':  42, 'm':  97.90540482   , 'm_std':   95.95 },
            'Tc' : { 'N':  43, 'm':  98.           , 'm_std':   98. },
            'Ru' : { 'N':  44, 'm': 101.9043441    , 'm_std':  101.07 },
            'Rh' : { 'N':  45, 'm': 102.9054980    , 'm_std':  102.90550 },
            'Pd' : { 'N':  46, 'm': 105.9034804    , 'm_std':  106.42 },
            'Ag' : { 'N':  47, 'm': 106.9050916    , 'm_std':  107.8682 },
            'Cd' : { 'N':  48, 'm': 113.90336509   , 'm_std':  112.414 },
            'In' : { 'N':  49, 'm': 114.903878776  , 'm_std':  114.818 },
            'Sn' : { 'N':  50, 'm': 119.90220163   , 'm_std':  118.710 },
            'Sb' : { 'N':  51, 'm': 120.9038120    , 'm_std':  121.760 },
            'Te' : { 'N':  52, 'm': 129.906222748  , 'm_std':  127.60 },
            'I'  : { 'N':  53, 'm': 126.9044719    , 'm_std':  126.90447 },
            'Xe' : { 'N':  54, 'm': 131.9041550856 , 'm_std':  131.293 },
            'Cs' : { 'N':  55, 'm': 132.9054519610 , 'm_std':  132.90545196 },
            'Ba' : { 'N':  56, 'm': 137.90524700   , 'm_std':  137.327 },
            'La' : { 'N':  57, 'm': 138.9063563    , 'm_std':  138.90547 },
            'Ce' : { 'N':  58, 'm': 139.9054431    , 'm_std':  140.116 },
            'Pr' : { 'N':  59, 'm': 140.9076576    , 'm_std':  140.90766 },
            'Nd' : { 'N':  60, 'm': 141.9077290    , 'm_std':  144.242 },
            'Pm' : { 'N':  61, 'm': 145.           , 'm_std':  145. },
            'Sm' : { 'N':  62, 'm': 151.9197397    , 'm_std':  150.36 },
            'Eu' : { 'N':  63, 'm': 152.9212380    , 'm_std':  151.964 },
            'Gd' : { 'N':  64, 'm': 157.9241123    , 'm_std':  157.25 },
            'Tb' : { 'N':  65, 'm': 158.9253547    , 'm_std':  158.92535 },
            'Dy' : { 'N':  66, 'm': 163.9291819    , 'm_std':  162.500 },
            'Ho' : { 'N':  67, 'm': 164.9303288    , 'm_std':  164.93033 },
            'Er' : { 'N':  68, 'm': 165.9302995    , 'm_std':  167.259 },
            'Tm' : { 'N':  69, 'm': 168.9342179    , 'm_std':  168.93422 },
            'Yb' : { 'N':  70, 'm': 173.9388664    , 'm_std':  173.054 },
            'Lu' : { 'N':  71, 'm': 174.9407752    , 'm_std':  174.9668 },
            'Hf' : { 'N':  72, 'm': 179.9465570    , 'm_std':  178.49 },
            'Ta' : { 'N':  73, 'm': 180.9479958    , 'm_std':  180.94788 },
            'W'  : { 'N':  74, 'm': 183.95093092   , 'm_std':  183.84 },
            'Re' : { 'N':  75, 'm': 186.9557501    , 'm_std':  186.207 },
            'Os' : { 'N':  76, 'm': 191.9614770    , 'm_std':  190.23 },
            'Ir' : { 'N':  77, 'm': 192.9629216    , 'm_std':  192.217 },
            'Pt' : { 'N':  78, 'm': 194.9647917    , 'm_std':  195.084 },
            'Au' : { 'N':  79, 'm': 196.96656879   , 'm_std':  196.966569 },
            'Hg' : { 'N':  80, 'm': 201.97064340   , 'm_std':  200.592 },
            'Tl' : { 'N':  81, 'm': 204.9744278    , 'm_std':  (204.382 + 204.385)/2 },
            'Pb' : { 'N':  82, 'm': 207.9766525    , 'm_std':  207.2 },
            'Bi' : { 'N':  83, 'm': 208.9803991    , 'm_std':  208.98040 },
            'Po' : { 'N':  84, 'm': 209.           , 'm_std':  209. },
            'At' : { 'N':  85, 'm': 210.           , 'm_std':  210. },
            'Rn' : { 'N':  86, 'm': 222.           , 'm_std':  222. },
            'Fr' : { 'N':  87, 'm': 223.0197360    , 'm_std':  223. },
            'Ra' : { 'N':  88, 'm': 226.           , 'm_std':  226. },
            'Ac' : { 'N':  89, 'm': 227.           , 'm_std':  227. },
            'Th' : { 'N':  90, 'm': 232.0380558    , 'm_std':  232.0377 },
            'Pa' : { 'N':  91, 'm': 231.0358842    , 'm_std':  231.03588 },
            'U'  : { 'N':  92, 'm': 238.0507884    , 'm_std':  238.02891 },
            'Np' : { 'N':  93, 'm': 237.           , 'm_std':  237. },
            'Pu' : { 'N':  94, 'm': 244.           , 'm_std':  244. },
            'Am' : { 'N':  95, 'm': 243.           , 'm_std':  243. },
            'Cm' : { 'N':  96, 'm': 247.           , 'm_std':  247. },
            'Bk' : { 'N':  97, 'm': 247.           , 'm_std':  247. },
            'Cf' : { 'N':  98, 'm': 251.           , 'm_std':  251. },
            'Es' : { 'N':  99, 'm': 252.           , 'm_std':  252. },
            'Fm' : { 'N': 100, 'm': 257.           , 'm_std':  257. },
            'Md' : { 'N': 101, 'm': 258.           , 'm_std':  258. },
            'No' : { 'N': 102, 'm': 259.           , 'm_std':  259. },
            'Lr' : { 'N': 103, 'm': 262.           , 'm_std':  262. },
            'Rf' : { 'N': 104, 'm': 267.           , 'm_std':  267. },
            'Db' : { 'N': 105, 'm': 268.           , 'm_std':  268. },
            'Sg' : { 'N': 106, 'm': 271.           , 'm_std':  271. },
            'Bh' : { 'N': 107, 'm': 272.           , 'm_std':  272. },
            'Hs' : { 'N': 108, 'm': 270.           , 'm_std':  270. },
            'Mt' : { 'N': 109, 'm': 276.           , 'm_std':  276. },
            'Ds' : { 'N': 110, 'm': 281.           , 'm_std':  281. },
            'Rg' : { 'N': 111, 'm': 280.           , 'm_std':  280. },
            'Cn' : { 'N': 112, 'm': 285.           , 'm_std':  285. },
            'Nh' : { 'N': 113, 'm': 284.           , 'm_std':  284. },
            'Fl' : { 'N': 114, 'm': 289.           , 'm_std':  289. },
            'Mc' : { 'N': 115, 'm': 288.           , 'm_std':  288. },
            'Lv' : { 'N': 116, 'm': 293.           , 'm_std':  293. },
            'Ts' : { 'N': 117, 'm': 292.           , 'm_std':  292. },
            'Og' : { 'N': 118, 'm': 294.           , 'm_std':  294. } }
Ptable.update({ 'D'  : { 'N':   1, 'm':   2.01410177812, 'm_std':    2.01410177812 },
                'T'  : { 'N':   1, 'm':   3.0160492779 , 'm_std':    3.0160492779  } })
Ptable_upper = { i.upper():j for i,j in Ptable.items() }



## AminoAcids ·························································
aa_letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I',
              'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
aa_letters.update({'HID':'H', 'HIE':'H', 'HIP':'H', 'CYX':'C', 'CYM':'C', 'LYN':'K', 'ASH':'D', 'GLH':'E'})
aa_letters.update({'HSD':'H', 'HSE':'H', 'HSP':'H'})
aa_letters.update({'SEC':'U'})

# aa set (more efficiency)
aa = { i for i in aa_letters.keys() }

# aa translations ( www.bmrb.wisc.edu/ref_info/atom_nom.tbl and GROMACS documentation $GMXPREFIX/share/gromacs/top/atom_nom.tbl )
rosetta_residues = { 'amber'  : [ 'HIS', 'HID', 'HIE', 'HIP', 'CYS', 'CYM', 'CYX', 'LYS', 'LYN', 'ASP', 'ASH', 'GLU', 'GLH' ],
                     'charmm' : [ 'HIS', 'HSD', 'HSE', 'HSP', 'CYS', 'CYS', 'CYS', 'LYS', 'LSN', 'ASP', 'ASP', 'GLU', 'GLU' ],
                     'opls'   : [ 'HIS', 'HIS', 'HIS', 'HIS', 'CYS', 'CYS', 'CYS', 'LYS', 'LYS', 'ASP', 'ASP', 'GLU', 'GLU' ]  }

rosetta_atoms = {
                  'dynamo' : { 'ALA' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'ARG' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'HD1', 'HD2', 'NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'ASN' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'OD1', 'ND2', 'HD21', 'HD22', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'ASP' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'OD1', 'OD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'CYS' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'SG', 'HG', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'GLN' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'OE1', 'NE2', 'HE21', 'HE22', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'GLU' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'OE1', 'OE2', 'HE2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'GLY' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA1', 'HA2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HIS' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'ILE' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23', 'CG1', 'HG11', 'HG12', 'CD', 'HD1', 'HD2', 'HD3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'LEU' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG', 'CD1', 'HD11', 'HD12', 'HD13', 'CD2', 'HD21', 'HD22', 'HD23', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'LYS' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'HD1', 'HD2', 'CE', 'HE1', 'HE2', 'NZ', 'HZ1', 'HZ2', 'HZ3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'MET' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'SD', 'CE', 'HE1', 'HE2', 'HE3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'PHE' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'HZ', 'CE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'PRO' : ['N', 'H1', 'H2', 'H3', 'CD', 'HD1', 'HD2', 'CG', 'HG1', 'HG2', 'CB', 'HB1', 'HB2', 'CA', 'HA', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'SER' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'OG', 'HG', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'THR' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23', 'OG1', 'HG1', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'TRP' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'CD1', 'HD1', 'NE1', 'HE1', 'CE2', 'CZ2', 'HZ2', 'CH2', 'HH2', 'CZ3', 'HZ3', 'CE3', 'HE3', 'CD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'TYR' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'OH', 'HH', 'CE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'VAL' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG1', 'HG11', 'HG12', 'HG13', 'CG2', 'HG21', 'HG22', 'HG23', 'C', 'O', 'OC1', 'OC2', 'OXT'],

                               'HID' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HIE' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HIP' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HSD' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HSE' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HSP' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'CYX' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'SG', 'C', 'O', 'OC1', 'OC2', 'OXT']
                             },
                  'gmx'    : { 'ALA' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'ARG' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'HD1', 'HD2', 'NE', 'HE', 'CZ', 'NH1', '1HH1', '2HH1', 'NH2', '1HH2', '2HH2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'ASN' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'OD1', 'ND2', '1HD2', '2HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'ASP' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'OD1', 'OD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'CYS' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'SG', 'HG', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'GLN' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'OE1', 'NE2', '1HE2', '2HE2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'GLU' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'OE1', 'OE2', 'HE2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'GLY' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA1', 'HA2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HIS' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'ILE' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG2', '1HG2', '2HG2', '3HG2', 'CG1', '1HG1', '2HG1', 'CD', 'HD1', 'HD2', 'HD3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'LEU' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG', 'CD1', '1HD1', '2HD1', '3HD1', 'CD2', '1HD2', '2HD2', '3HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'LYS' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'HD1', 'HD2', 'CE', 'HE1', 'HE2', 'NZ', 'HZ1', 'HZ2', 'HZ3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'MET' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'SD', 'CE', 'HE1', 'HE2', 'HE3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'PHE' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'HZ', 'CE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'PRO' : ['N', 'H1', 'H2', 'H3', 'CD', 'HD1', 'HD2', 'CG', 'HG1', 'HG2', 'CB', 'HB1', 'HB2', 'CA', 'HA', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'SER' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'OG', 'HG', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'THR' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG2', '1HG2', '2HG2', '3HG2', 'OG1', 'HG1', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'TRP' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'CD1', 'HD1', 'NE1', 'HE1', 'CE2', 'CZ2', 'HZ2', 'CH2', 'HH2', 'CZ3', 'HZ3', 'CE3', 'HE3', 'CD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'TYR' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'OH', 'HH', 'CE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'VAL' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG1', '1HG1', '2HG1', '3HG1', 'CG2', '1HG2', '2HG2', '3HG2', 'C', 'O', 'OC1', 'OC2', 'OXT'],

                               'HID' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HIE' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HIP' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HSD' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HSE' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HSP' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'CYX' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'SG', 'C', 'O', 'OC1', 'OC2', 'OXT']
                             },
                  'maestro': { 'ALA' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'ARG' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'ASN' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'OD1', 'ND2', 'HD21', 'HD22', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'ASP' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'OD1', 'OD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'CYS' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'SG', 'HG', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'GLN' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'OE1', 'NE2', 'HE21', 'HE22', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'GLU' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'OE1', 'OE2', 'HE2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'GLY' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA2', 'HA3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HIS' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'ILE' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23', 'CG1', 'HG12', 'HG13', 'CD1', 'HD11', 'HD12', 'HD13', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'LEU' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG', 'CD1', 'HD11', 'HD12', 'HD13', 'CD2', 'HD21', 'HD22', 'HD23', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'LYS' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'NZ', 'HZ1', 'HZ2', 'HZ3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'MET' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'SD', 'CE', 'HE1', 'HE2', 'HE3', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'PHE' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'HZ', 'CE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'PRO' : ['N', 'H1', 'H2', 'H3', 'CD', 'HD2', 'HD3', 'CG', 'HG2', 'HG3', 'CB', 'HB2', 'HB3', 'CA', 'HA', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'SER' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'OG', 'HG', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'THR' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23', 'OG1', 'HG1', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'TRP' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'NE1', 'HE1', 'CE2', 'CZ2', 'HZ2', 'CH2', 'HH2', 'CZ3', 'HZ3', 'CE3', 'HE3', 'CD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'TYR' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'OH', 'HH', 'CE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'VAL' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG1', 'HG11', 'HG12', 'HG13', 'CG2', 'HG21', 'HG22', 'HG23', 'C', 'O', 'OC1', 'OC2', 'OXT'],

                               'HID' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HIE' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HIP' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HSD' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HSE' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'HSP' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O', 'OC1', 'OC2', 'OXT'],
                               'CYX' : ['N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'SG', 'C', 'O', 'OC1', 'OC2', 'OXT']
                             }
                }

## Ions ·······························································
# resName : segment
ions = { 'CL' : 'CL',
         'NA' : 'NA',
         'SOD' : 'SOD',
        #  'CA' : 'CA',
         'MG' : 'MG',
         'K'  : 'K' ,
         'RB' : 'RB',
         'CS' : 'CS',
         'LI' : 'LI',
         'ZN' : 'ZN' }

## Solvent ····························································
# resName : segment
solvent = { 'SOL' : 'BOX' }

## Segments ···························································
segment_letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']


#######################################################################
##  DEFINITIONS                                                      ##
#######################################################################

# PDB Strict formatting
# ATOM/HETATM  serial  name   altLoc  resName  chainID  resSeq  iCode  x       y     z      occupancy  tempFactor  segment  element  charge
# 1-6          7-11    13-16  17      18-20    22       23-26   27     31-38  39-46  47-54  55-60      61-66       73-76    77-78    79-80

# PDB Column formatting (simple)
# ATOM/HETATM  #Atom  Atom  Res  [Chain]  #Res  X  Y  Z  Occup  TempFac  [Segment]

##  Main PDB class  ###################################################
class pdb:
    '''
        Main PDB class

        Attributes
        ----------
        title : str
        pdb : list of dic

        Properties
        ----------
        natoms : int
            number of atoms
        nres : int
            number of protein residues
        nres_tot : int
            number of total residues
        nseg : int
            number of segments
        segments : dic
            dictionary of segmenets: {seg:#res}
        sequence : list
            list of protein residues: [[#,XXX,X]]
        ligands : list
            list of ligands
        ssbons : list
            list of disulphide bonds: [[#resSeq,#resSeq]]
        protein_weight : float
            molecular average protein weight (Da)

        Methods
        -------
        read(file,strict)
            read pdb from file
        write(file,title,remark4,renum_atoms,onlyProtein,notProtein)
            write pdb to file
        write_fasta(file,gaps)
            write fasta to file
        write_intseq(file)
            write sequence in fDynamo's interaction format to file
        write_crd(file)
            write fDynamo's crd to file
        write_seq(file,variants,ssbonds)
            write fDynamo's seq to file
        write_ligand(ligand)
            write fDynamo's topology of ligand to file
        substitute(field,origin,destination,protectProtein,onlyProtein)
            change matching values of field
        clean_field(field,value)
            change all values of field to zero/empty/value
        remove(field,value)
            remove atoms that match field
        all2ATOM()
            change all ATOM/HETATM to ATOM
        translate_residues(destination,origin)
            translate resName between standards
        translate_names(destination,origin)
            translate atom names between standards
        canonical_order(canon)
            reorder protein atoms inside residues
        renum_atoms()
            renumerate all atoms from 1
        renum_res(continuous,protectProtein,guess_segments)
            renumerate residues, each group from 1: protein/ligands/solvent/ions
        guess_elements(keepknown)
            guess elements of all atoms
        guess_segments(keepknown,useChains)
            guess segments, one for each group: protein/ligands/solvent/ions
        weight(guess_elements,onlyProtein,monoisotopic)
            molecular weight (Da)
        center(guess_elements,center,monoisotopic)
            move system's COM to center
        cys2cyx()
            change resName of CYS to CYX if corresponds
        guess_his()
            change HIS to corresponding HID/HIE/HIP
        guess_glh()
            change GLU to GLH if corresponds
        guess_ash()
            change ASP to ASH if corresponds
        guess_lyn()
            change LYS to LYN if corresponds
        guess_protonres()
            change resName based on protonation for HIS/GLU/ASP/LYS
            shortcut to call guess_his()/guess_glh()/guess_ash()/guess_lyn()
    '''

    def __init__( self ):
        self.title = ''
        self.pdb = []

    def __del__( self ):
        del self.title
        del self.pdb

    ## read pdb from file ---------------------------------------------
    def read( self, file, strict=True ):
        '''Read PDB file'''

        # initialize object (in case of reuse)
        self.__init__()

        # open file
        with open( file, 'rt' ) as inp:
            inpdb = inp.readlines()
            inpdb = map(str.strip, inpdb)

        # create pdb list of dictionaries
        for line in inpdb:

            # finish if 'END'
            if line.startswith('END'): break

            # get title
            if line.startswith('TITLE'):
                self.title += line[6:].strip() + ' '
                continue

            # ignore no ATOM/HETATM lines
            if not line.startswith('ATOM') and not line.startswith('HETATM'): continue

            if strict:      # strict mode
                newline = { 'ATOM'       : str(line[:6]).strip(),
                            'serial'     : int(line[6:11]),
                            'name'       : str(line[12:16]).strip(),
                            'altLoc'     : str(line[16]).strip(),
                            'resName'    : str(line[17:20]).strip(),
                            'chainID'    : str(line[21]).strip(),
                            'resSeq'     : int(line[22:26]),
                            'iCode'      : str(line[26]).strip(),
                            'x'          : float(line[30:38]),
                            'y'          : float(line[38:46]),
                            'z'          : float(line[46:54]),
                            'occupancy'  : float(line[54:60]),
                            'tempFactor' : float(line[60:66]),
                            'segment'    : str(line[72:76]).strip(),
                            'element'    : str(line[76:78]).strip(),
                            'charge'     : str(line[78:80]).strip() }
            else:           # column mode
                line = line.split()
                # check if chainID column present
                try: int(line[4]); c = 0
                except: c = 1
                newline = { 'ATOM'       : str(line[0]),
                            'serial'     : int(line[1]),
                            'name'       : str(line[2]),
                            'resName'    : str(line[3]),
                            'resSeq'     : int(line[4+c]),
                            'x'          : float(line[5+c]),
                            'y'          : float(line[6+c]),
                            'z'          : float(line[7+c]),
                            'occupancy'  : float(line[8+c]),
                            'tempFactor' : float(line[9+c]),
                            'altLoc' : str(''), 'chainID' : str(''), 'iCode' : str(''),
                            'element' : str(''), 'charge' : str('') }
                # check if segment column present
                try: newline.update({ 'segment': str(line[10+c]) })
                except: newline.update({ 'segment' : str('') })

            self.pdb.append(newline)

    ## write pdb to file ----------------------------------------------
    def write( self, file, title=False, remark4=False, renum_atoms=True, onlyProtein=False, notProtein=False ):
        '''Write PDB file'''
        # renumerate atoms
        if renum_atoms: self.renum_atoms()
        # open file
        with open( file, 'wt' ) as outp:
            if remark4: outp.write("REMARK   4      COMPLIES WITH FORMAT V. 3.3, 21-NOV-2012\n")
            if title: outp.write("TITLE     {:70s}\n".format(self.title))
            for n in range(self.natoms):
                line = self.pdb[n].copy()
                # discart if not aminoacid and onlyProtein
                if onlyProtein and line['resName'] not in aa: continue
                # discart if aminoacid and notProtein
                if notProtein and line['resName'] in aa: continue
                # correct alignment of atom name
                if len(line['name']) == 3: line['name'] = ' ' + line['name']
                # format pdb
                # formatted_line = "{:<6s}{:>5d} {:^4s}{:1s}{:>3s} {:1s}{:>4d}{:1s}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}      {:<4s}{:>2s}{:<2s}" \
                # .format( line['ATOM'], line['serial'], line['name'], line['altLoc'], line['resName'], line['chainID'], line['resSeq'], line['iCode'], \
                # line['x'], line['y'], line['z'], line['occupancy'], line['tempFactor'], line['segment'], line['element'], line['charge'] )
                formatted_line = "{:<6s}{:>5d} {:^4s}{:1s}{:>3s} {:1s}{:>4.4}{:1s}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}      {:<4s}{:>2s}{:<2s}" \
                .format( line['ATOM'], line['serial'], line['name'], line['altLoc'], line['resName'], line['chainID'], str(line['resSeq']), line['iCode'], \
                line['x'], line['y'], line['z'], line['occupancy'], line['tempFactor'], line['segment'], line['element'], line['charge'] )
                # write
                outp.write( formatted_line + "\n" )
            outp.write("END\n")  # final 'END'

    ## write fasta to file --------------------------------------------
    def write_fasta( self, file, gaps=False ):
        '''Write fasta file'''
        if not gaps:
            seq = [ i[2] for i in self.sequence ]
        else:
            seq = []; n = 0
            for i in range(self.sequence[0][0], self.sequence[-1][0]+1):
                if self.sequence[n][0] == i:
                    seq.append(self.sequence[n][2])
                    n += 1
                else: seq.append('-')
        # open file
        with open( file, 'wt' ) as outp:
            # title
            if self.title: outp.write("> {:<s} \n".format(self.title))
            else: outp.write("> {:<s} \n".format("UNKNOWN PROTEIN"))
            # sequence
            nline = int(m.ceil(self.nres / 70))
            for n in range(nline):
                formatted_line = "{:<s}".format(''.join(seq[0+n*70:70+n*70]))
                outp.write( formatted_line + "\n" )

    ## write interaction sequence to file -----------------------------
    def write_intseq( self, file ):
        '''Write interaction sequence to file'''
        seq = self.sequence
        # open file
        with open( file, 'wt' ) as outp:
            # sequence format: X#
            for res in seq:
                formatted_line = "{:<s}{:<d}".format(res[2],res[0])
                outp.write( formatted_line + "\n" )

    ## write dynamo crd to file ---------------------------------------
    def write_crd( self, file ):
        '''Write dynamo crd file'''
        # open file
        with open(file, 'wt') as outp:
            # crd total properties
            outp.write("!===============================================================================\n")
            outp.write("{:<d} {:<d}  {:<d}\n".format(self.natoms, self.nres_tot, self.nseg))
            outp.write("!===============================================================================\n")
            outp.write("!Symmetry\t1\n!CUBIC \t\tLattice_(Angstroms)\n")
            outp.write("!===============================================================================\n")
            # initializations
            natoms = self.natoms
            segments = self.segments
            n = 0
            subsys = ''  # aka 'segment'
            nsubsys = 0
            # main crd loop
            while n < natoms:
                # subsystem statement
                if self.pdb[n]['segment'] != subsys:
                    subsys = self.pdb[n]['segment']
                    nsubsys += 1
                    outp.write("Subsystem {:>5d}  {:s}\n".format(nsubsys, subsys))
                    outp.write("{:>7d}\n".format(segments[subsys]))
                # residue statement
                resName = self.pdb[n]['resName']
                resSeq = self.pdb[n]['resSeq']
                outp.write("Residue {:>5d}  {:s}\n".format(resSeq, resName))
                # atom statement (count number of atoms in residue and print all)
                nres = n
                while nres < natoms and self.pdb[nres]['resName'] == resName and self.pdb[nres]['resSeq'] == resSeq : nres += 1
                outp.write("{:7d}\n".format(nres - n))
                while n < nres:
                    outp.write("{:>7d} {:<4s} {:>10d} {:>18.10f}{:>18.10f}{:>18.10f}\n"
                        .format(self.pdb[n]['serial'], self.pdb[n]['name'],
                        Ptable[self.pdb[n]['element']]['N'], self.pdb[n]['x'],
                        self.pdb[n]['y'], self.pdb[n]['z']))
                    n += 1

    ## write dynamo seq to file ---------------------------------------
    def write_seq( self, file, variants=True, ssbonds=True ):
        '''Write dynamo seq file'''
        # open file
        with open(file, 'wt') as outp:
            outp.write("Sequence\n{:>5d}\n".format(self.nseg))
            # initializations
            natoms = self.natoms
            segments = self.segments
            sequence = self.sequence
            n = 0
            subsys = ''  # aka 'segment'
            # main crd loop
            while n < natoms:
                # subsystem statement
                if self.pdb[n]['segment'] != subsys:
                    nline = 0
                    subsys = self.pdb[n]['segment']
                    outp.write("\nSubsystem  {:s}\n{:>5d}\n".format(subsys,segments[subsys]))
                # residues
                if self.pdb[n]['resName'] != self.pdb[n-1]['resName'] or self.pdb[n]['resSeq'] != self.pdb[n-1]['resSeq']:
                    nline += 1
                    outp.write("{:s} ; ".format(self.pdb[n]['resName']))
                    if nline % 13 == 0:
                        nline = 0
                        outp.write("\n")
                if n+1 < natoms and self.pdb[n+1]['segment'] != subsys:
                    if variants and subsys == 'A':
                        outp.write("\nVARIANT          N_TERMINAL  {:<s}  {:<d}".format(sequence[0][1],sequence[0][0]))
                        outp.write("\nVARIANT          C_TERMINAL  {:<s}  {:<d}".format(sequence[-1][1],sequence[-1][0]))
                        outp.write("\nEnd\n")
                    else:
                        outp.write("\nEnd\n")
                n += 1
            outp.write("\nEnd\n")
            if ssbonds:
                for ssbond in self.ssbonds:
                    outp.write("\nLink DISULPHIDE_BRIDGE A CYX {:>3d} A CYX {:>3d}".format(*ssbond))
            outp.write("\nEnd")

    ## write dynamo ligand opls to file -------------------------------
    def write_ligand( self, ligand ):
        '''Write dynamo ligand topology file'''
        file="ligand_"+ligand
        # get ligand lines index
        lig_index = [i for i, n in enumerate(self.pdb) if n['resName']==ligand]
        # bonds: calculate upper distance matrix, consider bond if below threshold
        bond_thr = 1.7  # bond threshold
        bonds = []
        for i in lig_index:
            for j in range(i+1, lig_index[-1]+1):
                dist = m.sqrt( (self.pdb[i]['x']-self.pdb[j]['x'])**2 + (self.pdb[i]['y']-self.pdb[j]['y'])**2 + (self.pdb[i]['z']-self.pdb[j]['z'])**2 )
                if dist < bond_thr and (self.pdb[i]['element'] or self.pdb[j]['element']) != 'H':
                    bonds.append([self.pdb[i]['name'],self.pdb[j]['name']])
        # improper dihedrals (central in third position, rest alphabetical)
        impropers = []
        for i in lig_index:
            name = self.pdb[i]['name']
            # atoms bonded to atom_i and alphabetically sorted
            bonds_i = list({k for j in bonds for k in j if name in j and name != k})
            # if it has 3 bonds, assumed to be sp2
            if len(bonds_i) == 3:
                bonds_sorted = sorted(bonds_i)
                bonds_sorted.insert(2,name)
                impropers.append(bonds_sorted)
        # open file
        with open(file, 'wt') as outp:
            # header
            outp.write("!-------------------------------------------------------------------------------\n")
            outp.write("Residue {:s}\n".format(ligand))
            outp.write("!-------------------------------------------------------------------------------\n")
            outp.write("! # Atoms, bonds and impropers.\n")
            outp.write(" {}   {}   {}\n""".format(len(lig_index),len(bonds),len(impropers)))
            # atoms
            for n in lig_index:
                if self.pdb[n]['element'] != '':
                    outp.write("{:<4s}  {:<2s}    0.0\n".format(self.pdb[n]['name'], self.pdb[n]['element']))
                else:
                    outp.write("{:<4s}  {:<2s}    0.0\n".format(self.pdb[n]['name'], list(self.pdb[n]['name'])[0]))
            # bonds
            outp.write("\n")
            nline = 0
            for b in bonds:
                nline += 1
                outp.write("{:<4s} {:<4s} ; ".format(b[0],b[1]))
                if nline % 6 == 0:
                    nline = 0
                    outp.write("\n")
            # impropers
            outp.write("\n\n")
            nline = 0
            for b in impropers:
                nline += 1
                outp.write("{:<4s} {:<4s} {:<4s} {:<4s} ; ".format(b[0],b[1],b[2],b[3]))
                if nline % 3 == 0:
                    nline = 0
                    outp.write("\n")

    ## substitute field -----------------------------------------------
    def substitute( self, field, origin, destination, protectProtein=False, onlyProtein=False ):
        '''Change field values that match'''
        if protectProtein and onlyProtein:
            return
        elif protectProtein:
            for n in self.pdb:
                if not n['resName'] in aa and n[field] == origin: n[field] = destination
        elif onlyProtein:
            for n in self.pdb:
                if n['resName'] in aa and n[field] == origin: n[field] = destination
        else:
            for n in self.pdb:
                if n[field] == origin: n[field] = destination

    ## empty/zero or value field --------------------------------------
    def clean_field( self, field, value=None ):
        '''Makes field zero/empty or value'''

        # check field
        if field not in {'ATOM', 'name', 'altLoc', 'resName', 'chainID', 'iCode', 'segment', 'element', 'charge', 'serial', 'resSeq', 'x', 'y', 'z', 'occupancy', 'tempFactor' }:
            raise ValueError("Not valid field")

        # check input value or empty
        if value != None:
            pass
        elif field in { 'ATOM', 'name', 'altLoc', 'resName', 'chainID', 'iCode', 'segment', 'element', 'charge' }:
            value = ''
        elif field in { 'serial', 'resSeq' }:
            value = 0
        elif field in { 'x', 'y', 'z', 'occupancy', 'tempFactor' }:
            value = 0.0

        # substitute field value
        for n in self.pdb: n[field] = value

    ## remove atom based on field -------------------------------------
    def remove( self, field, value ):
        '''Remove atoms based on match on a field'''
        n = 0
        while n < self.natoms:
            if self.pdb[n][field] == value:
                del self.pdb[n]
            else:
                n += 1

    ## all to ATOM ----------------------------------------------------
    def all2ATOM( self ):
        '''Change all ATOM/HETATM to ATOM'''
        for n in self.pdb:
            n['ATOM'] = 'ATOM'

    ## transate resNames between ff -----------------------------------
    def translate_residues( self, destination, origin=None ):
        '''Translate resNames between ff'''
        if origin is None:
            formats = rosetta_residues.keys()
        else:
            formats = [origin]
        for j in formats:
            if j == destination: continue
            for n in self.pdb:
                resName = n['resName']
                if resName in rosetta_residues[j]:
                    n['resName'] = rosetta_residues[destination][rosetta_residues[j].index(resName)]

    ## transate names between formats ---------------------------------
    def translate_names( self, destination, origin ):
        '''Translate atom names between formats'''
        for n in self.pdb:
            resName = n['resName']
            name = n['name']
            if resName in rosetta_atoms[origin] and name in rosetta_atoms[origin][resName]:
                n['name'] = rosetta_atoms[destination][resName][rosetta_atoms[origin][resName].index(name)]

    ## reorder protein atoms ------------------------------------------
    def canonical_order( self, canon='dynamo' ):
        '''Reorder protein atoms'''
        nres = self.nres
        # find protein
        for linen in range(self.natoms):
            if self.pdb[linen]['resName'] in aa: break
        res = 1
        while res <= nres:
            # residue number and name
            resSeq = self.pdb[linen]['resSeq']
            resName = self.pdb[linen]['resName']
            # find last residue line
            linen_end = linen
            while True and linen_end < self.natoms:
                if self.pdb[linen_end]['resSeq'] == resSeq and self.pdb[linen_end]['resName'] == resName: linen_end += 1
                else: break
            # reorder atoms (modified bubble sort algorithm)
            for j in range(linen, linen_end-1):
                swapped = False
                for i in range(linen, linen_end-j+linen-1):
                    try:
                        if rosetta_atoms[canon][resName].index(self.pdb[i]['name']) > rosetta_atoms[canon][resName].index(self.pdb[i+1]['name']):
                            self.pdb[i], self.pdb[i+1] = self.pdb[i+1], self.pdb[i]
                            swapped = True
                    except ValueError:
                        self.pdb[i], self.pdb[i+1] = self.pdb[i+1], self.pdb[i]
                        swapped = True
                if not swapped: break
            linen = linen_end +1
            res += 1

    ## renumerate atoms -----------------------------------------------
    def renum_atoms( self ):
        '''Renumerate all atoms from 1'''
        for n in range(self.natoms):
            self.pdb[n]['serial'] = n + 1

    ## renumerate residues --------------------------------------------
    def renum_res( self, continuous=False, protectProtein=False, guess_segments=False ):
        '''Renumerate residues, each group from 1: Protein / Ligands / Solvent / Ions'''
        if guess_segments: self.guess_segments()
        resSeq = 0
        resSeq_prev = 0
        resName_prev = ''
        segment_prev = ''
        if continuous:  # all residues in global list, ignore protectProtein
            for n in self.pdb:
                if n['resName'] != resName_prev or n['resSeq'] != resSeq_prev:
                    resName_prev = n['resName']
                    resSeq_prev = n['resSeq']
                    resSeq += 1
                n['resSeq'] = resSeq
        else:  # renumerate with a new list for each segment
            for n in self.pdb:
                if protectProtein and n['resName'] in aa: continue
                if n['segment'] != segment_prev:
                    segment_prev = n['segment']
                    resSeq = 0
                if n['resName'] != resName_prev or n['resSeq'] != resSeq_prev:
                    resName_prev = n['resName']
                    resSeq_prev = n['resSeq']
                    resSeq += 1
                n['resSeq'] = resSeq

    ## guess elements -------------------------------------------------
    def guess_elements( self, keepknown=True ):
        '''Guess elements for all atoms'''
        for n in self.pdb:
            if keepknown and n['element'] != '': continue
            name = ''.join([i for i in n['name'] if i.isalpha()])
            if n['resName'] in aa:  # protein
                if name[0] in Ptable:
                    n['element'] = name[0]
                elif name[1] in Ptable:
                    n['element'] = name[1]
            elif n['resName'] in ions:  # ions
                n['element'] = name[0] + name[1].lower()
            else:  # rest
                n['element'] = name[0]

    ## guess segments -------------------------------------------------
    def guess_segments( self, keepknown=True, useChains=False ):
        '''Guess segments for each group: Protein / Ligands / Solvent / Ions'''
        ligands = self.ligands
        for n in self.pdb:
            if not keepknown: n['segment'] = ''
            if useChains and n['chainID'] != '':
                n['segment'] = n['chainID']
                continue
            if keepknown and n['segment'] != '': continue
            resName = n['resName']
            if resName in aa:
                n['segment'] = segment_letters[0]
            elif resName in solvent:
                n['segment'] = solvent[resName]
            elif resName in ions:
                n['segment'] = "ION"
            elif resName in ligands:
                n['segment'] = segment_letters[ligands.index(resName) + 1]

    ## molecular weight -----------------------------------------------
    def weight( self, guess_elements=True, onlyProtein=True, monoisotopic=False ):
        '''Molecular weight (Da)'''
        if guess_elements: self.guess_elements()
        mass = 0.0
        if monoisotopic: atomic_mass = 'm'
        else: atomic_mass = 'm_std'
        for n in self.pdb:
            if not onlyProtein or n['resName'] in aa:
                mass += Ptable[ n['element'] ][atomic_mass]
        return mass

    ## move com to center ---------------------------------------------
    def center( self, guess_elements=True, center=[0., 0., 0.], monoisotopic=False ):
        '''Move system com to the center'''
        if guess_elements: self.guess_elements()
        if monoisotopic: atomic_mass = 'm'
        else: atomic_mass = 'm_std'
        m_total = self.weight(guess_elements=False, onlyProtein=False, monoisotopic=monoisotopic)
        com = [0.0, 0.0, 0.0]
        # calculate com
        for n in self.pdb:
            for i, j in zip( [0,1,2], ['x','y','z'] ):
                com[i] += n[j] * Ptable[n['element']][atomic_mass]
        com = [ i / m_total for i in com ]
        # move com
        for n in self.pdb:
            for i, j in zip( [0,1,2], ['x','y','z'] ):
                n[j] -= com[i] + center[i]

    ## change names CYS for CYX ---------------------------------------
    def cys2cyx( self ):
        '''Change resName of CYS for CYX if not HG'''
        # get resSeq for SG and HG atoms
        sg = { line['resSeq'] for n, line in enumerate(self.pdb) if line['name']=='SG' and line['resName'] in aa }
        hg = { line['resSeq'] for n, line in enumerate(self.pdb) if line['name']=='HG' and line['resName'] in aa }
        # change names
        for n in self.pdb:
            if n['resName'] == 'CYS' and n['resSeq'] in sg.difference(hg): n['resName'] = 'CYX'
            elif n['resName'] == 'CYX' and n['resSeq'] in sg.intersection(hg): n['resName'] = 'CYS'

    ## change HIS for HID/HIE/HIP -------------------------------------
    def guess_his( self ):
        '''Change HIS for corresponding HID/HIE/HIP'''
        # get resSeq for N HIS
        nhis = [ n for n, line in enumerate(self.pdb) if line['name']=='N' and line['resName']=='HIS' ]
        for linen in nhis:
            # residue number and name
            resSeq = self.pdb[linen]['resSeq']
            resName = self.pdb[linen]['resName']
            # find last residue line
            linen_end = linen
            while linen_end < self.natoms:
                if self.pdb[linen_end]['resSeq'] == resSeq and self.pdb[linen_end]['resName'] == resName: linen_end += 1
                else: break
            # determine HID/HIE/HIP
            keyres = [ line['name'] for line in self.pdb[linen:linen_end] if line['name']=='HD1' or line['name']=='HE2' ]
            if len(keyres)==2:
                histype = 'HIP'
            elif len(keyres)==1:
                if keyres[0]=='HD1': histype = 'HID'
                elif keyres[0]=='HE2': histype = 'HIE'
            else:
                raise ValueError('Strange HIS found')
            # change resName
            for n in range(linen,linen_end):
                self.pdb[n]['resName'] = histype

    ## change GLU for GLH ---------------------------------------------
    def guess_glh( self ):
        '''Change GLU for corresponding GLH'''
        # get resSeq for N HIS
        nglu = [ n for n, line in enumerate(self.pdb) if line['name']=='N' and line['resName']=='GLU' ]
        for linen in nglu:
            # residue number and name
            resSeq = self.pdb[linen]['resSeq']
            resName = self.pdb[linen]['resName']
            # find last residue line
            linen_end = linen
            while linen_end < self.natoms:
                if self.pdb[linen_end]['resSeq'] == resSeq and self.pdb[linen_end]['resName'] == resName: linen_end += 1
                else: break
            # determine GLH
            if 'HE2' in { line['name'] for line in self.pdb[linen:linen_end] }:
                # change resName
                for n in range(linen,linen_end):
                    self.pdb[n]['resName'] = 'GLH'

    ## change ASP for ASH ---------------------------------------------
    def guess_ash( self ):
        '''Change ASP for corresponding ASH'''
        # get resSeq for N HIS
        nasp = [ n for n, line in enumerate(self.pdb) if line['name']=='N' and line['resName']=='ASP' ]
        for linen in nasp:
            # residue number and name
            resSeq = self.pdb[linen]['resSeq']
            resName = self.pdb[linen]['resName']
            # find last residue line
            linen_end = linen
            while linen_end < self.natoms:
                if self.pdb[linen_end]['resSeq'] == resSeq and self.pdb[linen_end]['resName'] == resName: linen_end += 1
                else: break
            # determine ASH
            if 'HD2' in { line['name'] for line in self.pdb[linen:linen_end] }:
                # change resName
                for n in range(linen,linen_end):
                    self.pdb[n]['resName'] = 'ASH'

    ## change LYS for LYN ---------------------------------------------
    def guess_lyn( self ):
        '''Change LYS for corresponding LYN'''
        # get resSeq for N HIS
        nasp = [ n for n, line in enumerate(self.pdb) if line['name']=='N' and line['resName']=='LYS' ]
        for linen in nasp:
            # residue number and name
            resSeq = self.pdb[linen]['resSeq']
            resName = self.pdb[linen]['resName']
            # find last residue line
            linen_end = linen
            while linen_end < self.natoms:
                if self.pdb[linen_end]['resSeq'] == resSeq and self.pdb[linen_end]['resName'] == resName: linen_end += 1
                else: break
            # determine ASH
            if 'HZ3' not in { line['name'] for line in self.pdb[linen:linen_end] }:
                # change resName
                for n in range(linen,linen_end):
                    self.pdb[n]['resName'] = 'LYN'

    ## guess protonation resName of key residues ----------------------
    def guess_protonres( self ):
        '''Change residues based on their protonation (HIS,GLU,ASP,LYS)'''
        self.guess_his()
        self.guess_glh()
        self.guess_ash()
        self.guess_lyn()

    ## number of atoms ------------------------------------------------
    @property
    def natoms( self ):
        '''Number of atoms'''
        return len(self.pdb)

    ## number of protein residues -------------------------------------
    @property
    def nres( self ):
        '''Number of protein residues'''
        return len(self.sequence)

    ## number of total residues ---------------------------------------
    @property
    def nres_tot( self ):
        '''Number of total residues'''
        nres_tot = 1
        for n in range(1, self.natoms):
            if self.pdb[n]['resName'] != self.pdb[n-1]['resName'] or self.pdb[n]['resSeq'] != self.pdb[n-1]['resSeq']:
                nres_tot += 1
        return nres_tot

    ## number of segments ---------------------------------------------
    @property
    def nseg( self ):
        '''Number of segments'''
        return len(self.segments.keys())

    ## dictionary of segments -----------------------------------------
    @property
    def segments( self ):
        '''Dictionary of segments: { seg : #res }'''
        seg = {}
        for n in range(self.natoms):
            if self.pdb[n]['segment'] == '': raise ValueError("Not all segments assigned")
            if self.pdb[n]['segment'] not in seg:
                seg.update({ self.pdb[n]['segment'] : 1 })
            elif self.pdb[n]['resSeq'] != self.pdb[n-1]['resSeq'] or self.pdb[n]['resName'] != self.pdb[n-1]['resName']:
                seg[self.pdb[n]['segment']] += 1
        return seg

    ## residue sequence -----------------------------------------------
    @property
    def sequence( self ):
        '''List of protein residues: [[ #, XXX, X ]]'''
        letter_seq = []
        for n in self.pdb:
            if n['resName'] in aa:
             res = n['resSeq']; break
        for n in self.pdb:
            if n['resName'] in aa and n['resSeq'] >= res:
                res = n['resSeq'] + 1
                letter_seq.append([ n['resSeq'], n['resName'], aa_letters[n['resName']] ])
        return letter_seq

    ## ligands list ---------------------------------------------------
    @property
    def ligands( self ):
        '''List of ligands: [ XXX ]'''
        ligands = []
        for n in self.pdb:
            resName = n['resName']
            if resName not in aa and resName not in solvent \
               and resName not in ions and resName not in ligands:
                ligands.append(n['resName'])
        return ligands

    ## dishuphide bridges list ----------------------------------------
    @property
    def ssbonds( self ):
        '''S-S bonds list: [[ #resSeq, #resSeq ]]'''
        ss_bond_thr = 3.0  # maximum distance to consider S-S bond (A)
        # get lines from SG atoms from CYX
        sg = [ line for n, line in enumerate(self.pdb) if line['name']=='SG' and line['resName']=='CYX' ]
        # calculate upper distance matrix, consider bond if below threshold
        ss_bonds = []
        for i in range(len(sg)):
            for j in range(i+1, len(sg)):
                dist = m.sqrt( (sg[i]['x']-sg[j]['x'])**2 + (sg[i]['y']-sg[j]['y'])**2 + (sg[i]['z']-sg[j]['z'])**2 )
                if dist < ss_bond_thr:
                    ss_bonds.append([sg[i]['resSeq'],sg[j]['resSeq']])
        if len(ss_bonds) != len(sg)/2: raise ValueError("Issues with SS bonds")
        return ss_bonds

    ## protein weight -------------------------------------------------
    @property
    def protein_weight( self ):
        '''Molecular average protein weight (Da)'''
        return self.weight( guess_elements=True, onlyProtein=True, monoisotopic=False )


##  Standard protein preparation  #####################################
def standard_preparation( molec, inpformat, outformat, ff='amber' ):
    '''Standard protein preparation'''
    molec.all2ATOM()
    molec.translate_residues(ff, origin=None)
    molec.translate_names(outformat, origin=inpformat)
    molec.canonical_order(outformat)
    molec.guess_elements(keepknown=False)


#######################################################################
##  MAIN                                                             ##
#######################################################################
if __name__ == '__main__':

    # parser
    parser = __parserbuilder()
    args = parser.parse_args()

    # names of files
    infile    = args.I
    basein    = os.path.splitext(infile)[0]
    inpformat = args.i
    outformat = args.o
    inpff     = args.iff
    outff     = args.off
    if args.O is not None:
        outfile = args.O
        baseout = os.path.splitext(outfile)[0]
    else:
        if outformat == 'fasta':
            baseout = basein
            outfile = baseout + ".fasta"
        elif outformat == 'intseq':
            baseout = basein
            outfile = baseout + ".resid"
        else:
            baseout = basein + "_" + outformat
            outfile = baseout + ".pdb"

    # read pdb
    my_pdb = pdb()
    my_pdb.read(infile, strict=args.simple)

    # main selection
    if outformat == 'fasta':
        my_pdb.write_fasta(outfile, gaps=False)
    elif outformat == 'intseq':
        my_pdb.write_intseq(outfile)
    else:
        if outformat == 'gmx':
            my_pdb.substitute('name', 'HXT', 'OXT')
            standard_preparation(my_pdb, inpformat, outformat, outff)
            my_pdb.guess_protonres()
            if args.center: my_pdb.center(guess_elements=False, center=[0., 0., 0.], monoisotopic=False)
            my_pdb.write(outfile, title=True, remark4=True, renum_atoms=True, onlyProtein=True)

        elif outformat == 'dynamo':
            standard_preparation(my_pdb, inpformat, outformat, outff)
            my_pdb.cys2cyx()
            my_pdb.remove('name', 'OC2')
            my_pdb.substitute('name', 'OC1', 'O')
            my_pdb.substitute('resName', 'NA', 'SOD')
            my_pdb.guess_his()
            my_pdb.renum_res(continuous=False, protectProtein=False, guess_segments=True)
            my_pdb.renum_atoms()
            if args.center: my_pdb.center(guess_elements=False, center=[0., 0., 0.], monoisotopic=False)
            if args.dynamize:
                my_pdb.write_crd(baseout + '.crd')
                my_pdb.write_seq(baseout + '.seq', variants=True, ssbonds=True)
                for ligand in my_pdb.ligands: my_pdb.write_ligand(ligand)
            for n in ['occupancy', 'tempFactor', 'element', 'charge', 'chainID']: my_pdb.clean_field(n)
            my_pdb.write(outfile, title=True, remark4=True, renum_atoms=True, onlyProtein=False)

    sys.stderr.write("{:s} --> {:s} \n".format(infile, outfile))
