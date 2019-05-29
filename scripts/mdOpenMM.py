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
        self["pstream"] = None
        self["dyntstep"] = 0.002
        self["dynsteps"] = 0
        self["dynoutfrq"] = 5000
        self["dyntemp"] = 298.0
        self["lang"] = True
        self["langfbeta"] = 0.01
        self["boxsize"] = [None, None, None]
    def toppar(self):
        toppar = []
        toppar.extend(self['xpar'].split(":"))
        toppar.extend(self['xtop'].split(":"))
        if self['pstream'] is not None:
            toppar.extend(self['pstream'].split(":"))
        return tuple(toppar)
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
        if par[0].lower() != 'ca':
            sys.exit("Error: it is not supported cons type %s\n"%par[0])
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
        if len(x) >= 3:
            self.bottom_radius = float(x[2])
        else:
            self.bottom_radius = 0.0

def construct_restraint(psf, cons, ref):
    if cons.bottom_radius == 0.0:
        restraint = CustomExternalForce("k0*dsq ; dsq=((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    else:
        restraint = CustomExternalForce("k0*(max(d-d0, 0.0))^2 ; d=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        restraint.addGlobalParameter("d0", cons.bottom_radius*angstroms)

    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")
    restraint.addPerParticleParameter('k0')

    calphaIndex = []
    for i,atom in enumerate(psf.topology.atoms()):
        if atom.name == 'CA':
            calphaIndex.append(i)
    #
    atom_s = [atom for atom in ref.topology.atoms()]
    #
    k = -1
    for i,atom_crd in enumerate(ref.positions):
        if ref.atom[1][i] not in cons.atoms:
            continue
        #
        k += 1
        if ref.atom[0][i] in cons.residues:
            mass = atom_s[i].element.mass
            force_const = cons.force_const*mass*kilocalories_per_mole/angstroms**2
            parameters = list(atom_crd.value_in_unit(nanometers))
            parameters.append(force_const)
            restraint.addParticle(calphaIndex[k], parameters)
    return restraint

def run(arg, par, cons=None):
    #
    ff = CharmmParameterSet(*par.toppar())
    psf = CharmmPsfFile(arg.input_psf[0])
    psf.setBox(par['boxsize'][0], par['boxsize'][1], par['boxsize'][2])
    if arg.input_psf[1].endswith("crd"):
        pdb = CharmmCrdFile(arg.input_psf[1])
    elif arg.input_psf[1].endswith("pdb"):
        pdb = PDBFile(arg.input_psf[1])
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
        system.addForce(construct_restraint(psf, cons, ref))
    #
    simulation = Simulation(psf.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    if arg.restart is not None:
        with open(arg.restart, 'rb') as fp:
            simulation.context.loadCheckpoint(fp.read())
    else:
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
    #

if __name__=='__main__':
    main()
