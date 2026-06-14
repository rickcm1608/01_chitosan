from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import sys

# -------------------------
# Input files
# -------------------------
prmtop = AmberPrmtopFile('complex.prmtop')
rst = AmberInpcrdFile('equil4.rst')

# -------------------------
# System setup
# -------------------------
system = prmtop.createSystem(nonbondedMethod=PME,
                             nonbondedCutoff=1.0*nanometer,
                             constraints=HBonds)

# NPT barostat (300 K, 1 atm)
system.addForce(MonteCarloBarostat(1*atmosphere, 300*kelvin, 25))

# Langevin integrator (2 fs timestep)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# GPU platform
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('DeviceIndex', '2')

# -------------------------
# Simulation
# -------------------------
simulation = Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(rst.positions)
if rst.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*rst.boxVectors)

# -------------------------
# Reporters
# -------------------------
nsteps   = 25_000_000   # 50 ns at 2 fs timestep
dcd_freq = 25_000       # save every 50 ps → 1000 frames total
log_freq = 25_000

simulation.reporters.append(DCDReporter('traj_50ns.dcd', dcd_freq))
simulation.reporters.append(StateDataReporter(sys.stdout, log_freq, step=True,
    time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
    temperature=True, density=True, progress=True, remainingTime=True,
    speed=True, totalSteps=nsteps, separator='\t'))

simulation.reporters.append(
    StateDataReporter('prod50ns.log', log_freq, step=True, time=True,
        potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
        temperature=True, density=True, progress=True, remainingTime=True,
        speed=True, totalSteps=nsteps, separator='\t'))

# -------------------------
# Run production MD: 50 ns
# -------------------------
print("Starting 50 ns production simulation (chitosan + TIP4PEW water)...")
simulation.step(nsteps)
print("Simulation complete.")

# Save final snapshot
positions = simulation.context.getState(getPositions=True).getPositions()
with open('prod50ns.pdb', 'w') as f:
    PDBFile.writeFile(simulation.topology, positions, f)

simulation.saveState('prod50ns.state')
simulation.saveCheckpoint('prod50ns.chk')
