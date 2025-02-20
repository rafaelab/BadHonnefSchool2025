"""
This notebook exemplifies how to simulate the propagation of UHECRs over cosmological distances. 
The simulation is done assuming this spectrum to be of the form E^(-1).
In the post-processing stage, it will be reweighted to the desired spectrum.

(This simulations may take a couple of minutes.)
"""

from crpropa import *



# general options
nEvents = 30000
energyRange = (1 * EeV, 1000 * EeV) # energy range of emission
distanceRange = (0., 1000. * Mpc) # UHECR sources uniformly distributed in this range
A, Z = 1, 1
cmb = CMB()
ebl = IRB_Gilmore12()
outputFile = 'simulations/sim1D-A_%02i_Z_%02i.txt' % (A, Z)


# source distribution: uniform with power-law spectrum
position = SourceUniform1D(*distanceRange)
direction = SourceDirection(Vector3d(-1, 0, 0)) # emit in the -x direction (1D simulation)
redshift = SourceRedshift1D() # takes the positions and assign the corresponding redshifts
energySpectrum = SourcePowerLawSpectrum(*energyRange, -1)
particleType = SourceParticleType(nucleusId(A, Z))
source = Source()
source.add(position)
source.add(redshift) 
source.add(direction)
source.add(energySpectrum)
source.add(particleType)

# output
outputType = Output.Event1D
output = TextOutput(outputFile, outputType)
output.disable(output.CandidateTagColumn)
output.setEnergyScale(eV)
output.setLengthScale(Mpc)

# observer 
observerType = Observer1D()
observer = Observer()
observer.add(observerType)
observer.onDetection(output)

# interactions and energy losses
pdCMB = PhotoDisintegration(cmb) # not needed for protons
pdEBL = PhotoDisintegration(ebl) # not needed for protons
pppCMB = PhotoPionProduction(cmb)
pppEBL = PhotoPionProduction(ebl)
eppCMB = ElectronPairProduction(cmb)
eppEBL = ElectronPairProduction(ebl)
nd = NuclearDecay()
z = Redshift()
processes = [z, pdCMB, pdEBL, pppCMB, pppEBL, eppCMB, eppEBL, nd]

# propagator: one-dimensional with minimum step 0.1 kpc and maximum step 1 Mpc
propagator = SimplePropagation(0.1 * kpc, 1 * Mpc)

# break conditions
breakEnergy = MinimumEnergy(1 * EeV)

# assemble simulation components
sim = ModuleList()
sim.add(propagator)
for interaction in processes:
	sim.add(interaction)
sim.add(observer)
sim.add(breakEnergy)
sim.setShowProgress(True)
sim.run(source, nEvents, True)

output.close()
