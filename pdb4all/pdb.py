"""
=======================================================================
  Main PDB class
=======================================================================

  Classes
  -------
    pdb

  Functions
  ---------
    standard_preparation

"""

import math as m
from copy import deepcopy

from pdb4all.constants import *
from pdb4all.rosetta import *


# PDB Strict formatting
# ATOM/HETATM  serial  name   altLoc  resName  chainID  resSeq  iCode  x       y     z      occupancy  tempFactor  segment  element  charge
# 1-6          7-11    13-16  17      18-20    22       23-26   27     31-38  39-46  47-54  55-60      61-66       73-76    77-78    79-80

# PDB Column formatting (simple)
# ATOM/HETATM  #Atom  Atom  Res  [Chain]  #Res  X  Y  Z  Occup  TempFac  [Segment]


##  Main PDB class  ###################################################
class PDB:
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
        __init__(file)
            initalization and optional pdb read
        read(file,strict)
            read pdb from file
        read_crd(file)
            read dynamo's crd from file
        read_xyz(file)
            read XYZ from file
        write(file,title,remark4,renum_atoms,onlyProtein,notProtein)
            write pdb to file
        write_fasta(file,gaps)
            write fasta to file
        write_xyz(file)
            write XYZ to file
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

    atom_empty = {
        'ATOM'       : "ATOM",
        'serial'     : 0,
        'name'       : "",
        'altLoc'     : "",
        'resName'    : "",
        'chainID'    : "",
        'resSeq'     : 0,
        'iCode'      : "",
        'x'          : 0.0,
        'y'          : 0.0,
        'z'          : 0.0,
        'occupancy'  : 0.0,
        'tempFactor' : 0.0,
        'segment'    : "",
        'element'    : "",
        'charge'     : ""
        }

    def __init__( self, file=None ):
        '''Class initialization'''
        self.title = ''
        self.pdb = []
        if file is not None:
            self.read(file)

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
            inpdb = [line.strip() for line in inpdb if line.strip()]

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
                            'altLoc' : '', 'chainID' : '', 'iCode' : '',
                            'element' : '', 'charge' : '' }
                # check if segment column present
                try: newline.update({ 'segment': str(line[10+c]) })
                except: newline.update({ 'segment' : str('') })

            self.pdb.append(newline)

    ## read crd from file ---------------------------------------------
    def read_crd( self, file, ):
        '''Read dynamo's crd from file'''

        # initialize object (in case of reuse)
        self.__init__()

        # open file
        with open( file, 'rt' ) as inp:
            incrd = inp.readlines()
            incrd = [line.strip() for line in incrd if line.strip() and not line.startswith('!')]

        # convert to pdb list of dictionaries
        a = deepcopy(self.atom_empty)
        for line in incrd:
            line = line.split("!")[0].split()   # to list and remove trailing comments
            if line[0].lower() == "subsystem":
                a['segment'] = str(line[2])
            elif line[0].lower() == "residue":
                a['resSeq']  = int(line[1])
                a['resName'] = str(line[2])
            elif len(line) != 6:
                continue
            else:
                a['serial']  = int(line[0])
                a['name']    = str(line[1])
                a['element'] = Ptable_inv_n[int(line[2])]
                a['x']       = float(line[3])
                a['y']       = float(line[4])
                a['z']       = float(line[5])
                self.pdb.append(deepcopy(a))

    ## read XYZ from file ---------------------------------------------
    def read_xyz( self, file, ):
        '''Read XYZ from file'''

        # initialize object (in case of reuse)
        self.__init__()

        # open file
        with open( file, 'rt' ) as inp:
            inxyz = inp.readlines()
            inxyz = [line.strip() for line in inxyz if line.strip()]

        # number of atoms and comment line
        natoms = int(inxyz.pop(0))
        self.title = str(inxyz.pop(0))

        # convert to pdb list of dictionaries
        a = deepcopy(self.atom_empty)
        for n in range(natoms):
            line = inxyz[n].split()
            a['element'] = str(line.pop(0))
            a['x'], a['y'], a['z'] = map(float, line)
            self.pdb.append(deepcopy(a))

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

    ## write XYZ to file ----------------------------------------------
    def write_xyz( self, file ):
        '''Write XYZ file'''
        self.guess_elements()
        with open( file, 'wt' ) as outp:
            outp.write(" {}\n".format(self.natoms))
            outp.write(" {}\n".format(self.title))
            for atom in self.pdb:
                outp.write(" {:<6s}  {:>18.10f} {:>18.10f} {:>18.10f}\n"
                           .format(atom['element'], atom['x'], atom['y'], atom['z']))

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
        file=ligand+".ff"
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
