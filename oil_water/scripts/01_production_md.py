from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import sys

# ===============================
# Load AMBER topology and restart
# ===============================
prmtop = AmberPrmtopFile('complex.prmtop')
rst    = AmberInpcrdFile('equil4.rst')

# ===============================
# Build system
# ===============================
system = prmtop.createSystem(
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometer,
    constraints=HBonds
)

# NPT barostat: 1 atm, 300 K, every 25 steps
system.addForce(MonteCarloBarostat(1*atmosphere, 300*kelvin, 25))

# Langevin integrator: 300 K, gamma=1 ps^-1, dt=2 fs
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# ===============================
# Platform (GPU)
# ===============================
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('DeviceIndex', '2')

# ===============================
# Simulation
# ===============================
simulation = Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(rst.positions)
if rst.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*rst.boxVectors)
simulation.context.setVelocitiesToTemperature(300*kelvin)

# ===============================
# Reporters
# ===============================
nsteps   = 25_000_000   # 50 ns  (dt=2 fs → 25M steps)
dcd_freq = 25_000       # save every 50 ps → 1000 frames total

simulation.reporters.append(DCDReporter('traj_50ns.dcd', dcd_freq))
simulation.reporters.append(StateDataReporter(
    'prod50ns.log', dcd_freq,
    step=True, time=True, potentialEnergy=True,
    kineticEnergy=True, temperature=True, volume=True, density=True
))

print("Starting 50 ns production MD...")
simulation.step(nsteps)
print("Done. Output: traj_50ns.dcd, prod50ns.log")
