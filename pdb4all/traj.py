"""
=======================================================================
  Trajectory class
=======================================================================

  Classes
  -------
    Traj

"""

import math as m
import sys
from copy import deepcopy
from struct import error, pack, unpack

from pdb4all import PDB


##  Trajectory class  #################################################
class Traj:
    '''
        Trajectory class

        Attributes
        ----------
        top : PDB
        frames_xyz : list
        nfixed : int
        nfree : int
        qcrys : int
        sele : list

        Properties
        ----------
        natoms : int
            number of atoms
        nframes : int
            number of frames

        Methods
        -------
        read_dcd(file, top, stride, top_format)
            load trajectory from DCD file
        nframe(nframe)
            return PDB object of nframe
        slice(start, stop, step)
            slice the trajectory
        join(traj)
            join trajectory
        rmsd(ndx, ndx)
            calculate RMSD respect the first frame of trajectory
    '''

    def __init__(self, top=None) -> None:
        self.top = top if type(top) == PDB else PDB()
        self.frames_xyz = []
        self.nfixed = 0
        self.nfree = 0
        self.qcrys = 0
        self.sele = []
        self._end = ""
        self._natoms_dcd = 0
        self._nframes_dcd = 0

    def __repr__(self) -> str:
        return '<Traj: {:d} frames :: {:d} atoms >' \
            .format(self.nframes, self.natoms)

    def __bool__(self) -> bool:
        return bool(self.nframes)

    def __len__(self) -> int:
        return self.nframes

    def __iter__(self) -> 'PDB':
        top = deepcopy(self.top)
        for frame_xyz in self.frames_xyz:
            for natom in range(self.natoms):
                top.pdb[natom]['x'], top.pdb[natom]['y'], top.pdb[natom]['z'] = frame_xyz[natom]
            yield top

    def deepcopy(self) -> 'Traj':
        return deepcopy(self)

    ## read dcd -------------------------------------------------------
    def read_dcd(self, file, top=None, stride=1, top_format=None) -> None:
        '''Load trajectory from DCD file'''
        # topology
        if not (top or self.top):
            raise Exception("Topology is required")
        elif type(top) == PDB:
            self.__init__(deepcopy(top))
        elif type(top) == str:
            self.__init__()
            self.top.read(top, format=top_format)
        elif self.top:
            self.__init__(self.top)
        # trajectory
        with open(file, 'rb') as f:
            # header
            self._end = "<" if unpack("i", f.read(4))[0] == 84 else ">"
            endi = f"{self._end}i"
            f.read(4)
            self._nframes_dcd = unpack(endi, f.read(4))[0]
            f.read(28)
            self.nfixed = unpack(endi, f.read(4))[0]
            f.read(4)
            self.qcrys = unpack(endi, f.read(4))[0]
            f.read(40)
            f.read(unpack(endi, f.read(4))[0] + 8)
            self._natoms_dcd = unpack(endi, f.read(4))[0]
            if self._natoms_dcd != self.top.natoms:
                raise Exception("Number of atoms in topology and trajectory do not match")
            self.nfree = self._natoms_dcd - self.nfixed
            f.read(4)
            if self.nfixed > 0:
                f.read(4)
                for natom in unpack(f"{self._end}{self.nfree:d}i", f.read(4 * self.nfree)):
                    self.sele.append(natom-1)
                f.read(4)
            # frames
            frames_xyz = self.top.xyz.copy()
            try:
                for nframe in range(self._nframes_dcd):
                    f.read(4)
                    if self.qcrys:
                        f.read(56)
                    if self.nfixed > 0 and nframe > 0:
                        natoms = self.nfree
                        natom_list = [self.sele[natom] for natom in range(self.nfree)]
                    else:
                        natoms = self._natoms_dcd
                        natom_list = range(self._natoms_dcd)
                    x_coord = unpack(f"{self._end}{natoms:d}f", f.read(4 * natoms))
                    f.read(8)
                    y_coord = unpack(f"{self._end}{natoms:d}f", f.read(4 * natoms))
                    f.read(8)
                    z_coord = unpack(f"{self._end}{natoms:d}f", f.read(4 * natoms))
                    f.read(4)
                    if nframe % stride == 0:
                        for natom in range(natoms):
                            frames_xyz[natom_list[natom]] = x_coord[natom], y_coord[natom], z_coord[natom]
                        self.frames_xyz.append(frames_xyz.copy())
            except error:
                sys.stderr.write(f"WARNING: {self.nframes} frames read instead of {self._nframes_dcd} claimed by header\n ")

    ## get nframe -----------------------------------------------------
    def nframe(self, nframe) -> 'PDB':
        '''Return PDB object of nframe'''
        if nframe >= self.nframes:
            raise ValueError(f"Frame requested ({nframe}) greater than number of frames ({self.nframes})")
        top = deepcopy(self.top)
        for natom in range(self.natoms):
            top.pdb[natom]['x'], top.pdb[natom]['y'], top.pdb[natom]['z'] = self.frames_xyz[nframe][natom]
        return top

    ## slice ----------------------------------------------------------
    def slice(self, start, stop, step=1) -> 'Traj':
        '''Slice the trajectory'''
        traj_slice = deepcopy(self)
        traj_slice.frames_xyz = traj_slice.frames_xyz[start:stop:step]
        return traj_slice

    ## join -----------------------------------------------------------
    def join(self, traj) -> 'Traj':
        '''Join trajectory'''
        if self.top != traj.top:
            raise Exception("Topologies do not match")
        traj_join = deepcopy(self)
        traj_join.frames_xyz.extend(traj.frames_xyz)
        return traj_join

    ## calculate RMSD -------------------------------------------------
    def rmsd(self, ndx=None) -> float:
        '''Calculate RMSD respect the first frame of trajectory'''
        atoms_ndx = ndx or range(self.natoms)
        return m.sqrt(sum(
            [(self.frames_xyz[0][natom][coord] - frame_xyz[natom][coord])**2
            for coord in range(3)
            for natom in atoms_ndx
            for frame_xyz in self.frames_xyz])
            / (self.nframes * len(atoms_ndx)))

    ## number of atoms ------------------------------------------------
    @property
    def natoms(self) -> int:
        return self.top.natoms

    ## number of frames -----------------------------------------------
    @property
    def nframes(self) -> int:
        return len(self.frames_xyz)
