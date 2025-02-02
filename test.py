from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdForceFieldHelpers

from rdkit.Chem import PyMol
import copy
import numpy as np

from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

inpcrd = AmberInpcrdFile('diala.inpcrd')
prmtop = AmberPrmtopFile('diala.prmtop', periodicBoxVectors=inpcrd.boxVectors)
system = prmtop.createSystem(nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer,
        constraints=HBonds)
integrator = VerletIntegrator(0.001*picoseconds)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(10000)