"""
=======================================================================
  Main
=======================================================================

"""

import os
import sys
import argparse

from pdb4all import __version__
from pdb4all.pdb import PDB, standard_preparation


##  PARSER  ###########################################################
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
                                     description=' -- Convert between common pdb formats and names --\n\n' +
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
                        help='read simplified pdb based on common columns\n' +
                        '  ATOM/HETATM  #Atom  Atom  Res  [Chain]  #Res  X  Y  Z  Occup  TempFac  [Segment]')
    return parser


##  MAIN  #############################################################

def main():
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

    # read crd or pdb
    my_pdb = PDB()
    if infile.split('.')[-1] == 'crd':
        my_pdb.read_crd(infile)
    else:
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
            my_pdb.substitute('name', 'OC1', 'O')
            my_pdb.substitute('name', 'OC2', 'OXT')
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


if __name__ == "__main__":
    main()
