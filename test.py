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
restraint = PeriodicTorsionForce()
system.addForce(restraint)
restraint.addTorsion(t1, t2, t3, t4, 1, restraintAngle * radians, 100*kilojoules_per_mole) 
integrator = VerletIntegrator(0.001*picoseconds)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(10000)

for bond in rotatableBonds:
    for i in range(36):
        addRestraint()
        runSimulation()
        conformerList.append(simulationResult)
finalConformations = filterResults(conformerList)

restraintEnergy = (restraint_angle-actual_angle)*force_constant
restraintEnergy = absolute_value((30-100)*100) = 0
restraintEnergy = absolute_value((30-100)*100) = 7000