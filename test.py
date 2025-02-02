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

# 90 degree dihedral angle
point1 = [0, 0, 0]
point2 = [1, 0, 0]
point3 = [0, 1, 0]
point4 = [0, 0, 1]

angle = angle_of_dihydral(point1, point2, point3, point4)
print(f"The dihedral angle is {angle:.2f} degrees")

# 45 degree dihedral angle
p1 = [0, 0, 0]
p2 = [1, 0, 0]
p3 = [0, 1, 1]
p4 = [1, 1, 0]

angle = angle_of_dihedral(p1, p2, p3, p4)
print(f"The dihedral angle is {angle:.2f} degrees")

# 120 degree dihedral angle
pt1 = [0, 0, 0]
pt2 = [1, 0, 0]
pt3 = [0, 2, 1]
pt4 = [2, 2, 0]

angle = angle_of_dihedral(pt1, pt2, pt3, pt4)
print(f"The dihedral angle is {angle:.2f} degrees")