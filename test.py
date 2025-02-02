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


# def vector_from_points(point1, point2):
  #  return np.array(point2) - np.array(point1)
  #  print(point1)
  #  print(point1.shape)
  #  print(point2)
  #  print(point2.shape)

# def cross_product(vector1, vector2):
  #  return np.cross(vector1, vector2)
  #  print(vector1)
  #  print(vector1.shape)
  #  print(vector2)
  #  print(vector2.shape)

# def dot_product(vector1, vector2):
  #  return np.dot(vector1, vector2)
  #  print(vector1)
  #  print(vector1.shape)

# def dihedral_angle(point1, point2, point3, point4):
  #  bond1 = vector_from_points(point1, point2)
  #  bond2 = vector_from_points(point2, point3)
  #  bond3 = vector_from_points(point3, point4)

  #  normal1 = cross_product(bond1, bond2)
  #  normal2 = cross_product(bond2, bond3)

  #  normal1 /= np.linalg.norm(normal1)
  #  normal2 /= np.linalg.norm(normal2)

# 90 degree dihdral angle
point1 = [0, 0, 0,]
point2 = [1, 0, 0]
point3 = [0, 1, 0]
point4 = [0, 0, 1]

angle = angle_of_dihedral(point1, point2, point3, point4)
print(f"The dihedral angle is {angle:.2f} degrees")