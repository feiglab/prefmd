#!/usr/bin/env python

import os
import sys
import argparse

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

class CHARMMparam(dict):
    def __init__(self):
        self["xpar"] = None
        self["xtop"] = None
        self["dyntstep"] = 0.002
        self["dynsteps"] = 0
        self["dynoutfrq"] = 5000
        self["dyntemp"] = 298.0
        self["lang"] = True
        self["langfbeta"] = 0.01
        self["boxsize"] = [None, None, None]
    def __getitem__(self, key):
        if key in self:
            return dict.__getitem__(self, key)
        else:
            return None
    def status(self):
        if self['xpar'] is None or self['xtop'] is None:
            return False
        elif None in self['boxsize']:
            return False
        else:
            return True
    def setParameter(self, key, p):
        if key in ['dyntstep', 'langfbeta', 'dyntemp']:
            self[key] = float(p)
        elif key in ['dynsteps', 'dynoutfrq']:
            self[key] = int(p)
        elif key in ['lang']:
            self[key] = (int(p) == 1)
        elif key in ['boxsize']:
            self[key] = (float(p)/10.0, float(p)/10.0, float(p)/10.0) # in nm
        elif key.startswith("box") and key[3] in ['x','y','z']:
            self['boxsize'][['x','y','z'].index(key[3])] = float(p)/10.0 # in nm 
        else:
            self[key] = p
    @classmethod
    def parse(cls, line):
        param = cls()
        #
        pars = line.strip().split(",")
        for par in pars:
            if '=' in par:
                key,p = par.split("=",1)
            else:
                key = par
                p = True
            if p == '':
                sys.stderr.write("A parameter key (%s) is presented, but corresponding value is not given\n"%key)
                continue
            if key in param or key in ['boxx','boxy','boxz']:
                param.setParameter(key, p)
        return param

class CHARMMcons:
    def __init__(self, par):
        if par[0].lower() in ['ca','cb', 'heavy']:
            self.atoms = [par[0].upper()]
        elif par[0].lower() == 'cab':
            self.atoms = ['CA', 'CB']
        elif par[0].lower() == 'cabp':
            self.atoms = [''] # TODO
        self.ref_pdb = par[1]
        #
        x = par[2].split("_")
        r = x[0].split(":")
        self.residues = range(int(r[0]), int(r[1])+1)
        self.force_const = float(x[1])

def run(arg, par, cons=None):
    ff = CharmmParameterSet(par['xtop'], par['xpar'])
    psf = CharmmPsfFile(arg.input_psf[0])
    psf.setBox(par['boxsize'][0], par['boxsize'][1], par['boxsize'][2])
    pdb = CharmmCrdFile(arg.input_psf[1])
    #pdb = PDBFile(arg.input_pdb)
    if cons is not None:
        if cons.ref_pdb == 'self':
            ref = pdb
            ref.atom = [ref.resno, ref.attype]
        else:
            ref = PDBFile(cons.ref_pdb)
            ref.atom = [[int(atom.residue.id) for atom in ref.topology.atoms()],\
                        [atom.name for atom in ref.topology.atoms()]]
    #
    system = psf.createSystem(ff, nonbondedMethod=PME, \
                              switchDistance=0.8*nanometers,\
                              nonbondedCutoff=1.0*nanometers,\
                              constraints=HBonds)

    if par['lang']:
        integrator = LangevinIntegrator(par['dyntemp']*kelvin,\
                                        par['langfbeta']/picosecond,\
                                        par['dyntstep']*picosecond)
    else:
        integrator = None
    #
    if cons is not None:
        restraint = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        restraint.addGlobalParameter("k", cons.force_const*kilocalories_per_mole/angstroms**2)
        restraint.addPerParticleParameter("x0")
        restraint.addPerParticleParameter("y0")
        restraint.addPerParticleParameter("z0")
        n = 0
        for i,atom_crd in enumerate(ref.positions):
            if ref.atom[0][i] in cons.residues and \
                    (ref.atom[1][i] in cons.atoms or \
                     (cons.atoms[0] == 'HEAVY' and ref.atom[1][i][0] in ['C','N','O','S'])):
                restraint.addParticle(i, atom_crd.value_in_unit(nanometers))
                n += 1
        system.addForce(restraint)
    #
    simulation = Simulation(psf.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    if arg.restart is not None:
        with open(arg.restart, 'rb') as fp:
            simulation.context.loadCheckpoint(fp.read())

    simulation.context.setVelocitiesToTemperature(par['dyntemp']*kelvin)
    if arg.trajout is not None:
        simulation.reporters.append(DCDReporter(arg.trajout, par['dynoutfrq']))
    if arg.final is not None:
        simulation.reporters.append(PDBReporter(arg.final, par['dynsteps']))

    if arg.log is None:
        logfile = sys.stdout
    else:
        logfile = arg.log
    simulation.reporters.append(StateDataReporter(logfile, par['dynoutfrq'], step=True, \
        time=True, kineticEnergy=True, potentialEnergy=True, temperature=True, progress=True, \
        remainingTime=True, speed=True, totalSteps=par['dynsteps'], separator='\t'))

    simulation.step(par['dynsteps'])

    with open(arg.restout, 'wb') as fout:
        fout.write(simulation.context.createCheckpoint())

def main():
    arg = argparse.ArgumentParser(prog='mdOpenMM')
    arg.add_argument(dest='input_pdb', metavar='PDBfile',\
            help='Initial structure in PDB format')
    arg.add_argument("-par", dest='par', default=None,\
            help='CHARMMparams')
    arg.add_argument("-cons", dest='cons', default=None, nargs='*', \
            help='[ca|cb|cab|cabp|heavy] refpdb|self min:max[_force][=...]')
    arg.add_argument('-psf', dest='input_psf', default=None, required=True, nargs=2)
    arg.add_argument("-restart", dest='restart', default=None)
    arg.add_argument("-restout", dest='restout', default=None)
    arg.add_argument("-trajout", dest='trajout', default=None)
    arg.add_argument("-log",     dest='log',     default=None)
    arg.add_argument("-final",   dest='final',   default=None)
    if len(sys.argv) == 1:
        arg.print_help()
        return
    #
    arg = arg.parse_args()
    #
    par = CHARMMparam.parse(arg.par)
    if not par.status():
        sys.stderr.write("Error: incomplete parameters\n")
        return
    #
    if arg.cons is not None:
        cons = CHARMMcons(arg.cons)
    else:
        cons = None
    #
    run(arg, par, cons=cons)

if __name__=='__main__':
    main()
