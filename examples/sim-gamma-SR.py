"""
This notebook exemplifies how to simulate the propagation of gamma rays over cosmological distances. 
The simulation is done assuming this spectrum to be of the form E^(-1).
In the post-processing stage, it will be reweighted to the desired spectrum.
This case is the standard (Lorentz invariant approach), using only CRPropa.

(This simulations may take a couple of minutes.)
"""

from crpropa import *

# general options
nEvents = 10000
energyMinimum = 1 * GeV
energyMaximum = 400 * TeV
redshift = 0.14
distance = redshift2ComovingDistance(redshift)
electrons = photons = True
thinning = 1.
cmb = CMB()
ebl = IRB_Gilmore12()
outputFile = 'simulations/sim1D-gamma-SR.txt'

# source distribution: uniform with power-law spectrum
position = SourcePosition(Vector3d(distance, 0, 0))
direction = SourceDirection(Vector3d(-1, 0, 0)) # emit in the -x direction (1D simulation)
redshifts = SourceRedshift1D() # takes the positions and assign the corresponding redshifts
energySpectrum = SourcePowerLawSpectrum(energyMinimum, energyMaximum, -1)
particleType = SourceParticleType(22) # we are interested in gamma rays
source = Source()
source.add(position)
source.add(redshifts) 
source.add(direction)
source.add(energySpectrum)
source.add(particleType)

# output
outputType = Output.Event1D
output = TextOutput(outputFile, outputType)
output.disable(output.CandidateTagColumn)
output.enable(output.WeightColumn) # since we are using thinning
output.setEnergyScale(eV)
output.setLengthScale(Mpc)

# observer 
observerType = Observer1D()
observer = Observer()
observer.add(observerType)
observer.onDetection(output)

# interactions
ppCMB = EMPairProduction(cmb, electrons, thinning)
ppEBL = EMPairProduction(ebl, electrons, thinning)
icsCMB = EMInverseComptonScattering(cmb, photons, thinning)
icsEBL = EMInverseComptonScattering(ebl, photons, thinning)
z = Redshift()
processes = [z, ppCMB, ppEBL, icsCMB, icsEBL]

# propagator: one-dimensional
propagator = SimplePropagation(0.1 * kpc, 100 * kpc)

# break conditions
breakEnergy = MinimumEnergy(1 * GeV)

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