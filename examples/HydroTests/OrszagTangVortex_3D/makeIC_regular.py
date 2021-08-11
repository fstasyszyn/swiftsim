################################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
#               2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#               2021 Federico Stasyszyn (fstasyszyn@unc.edu.ar)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

import h5py
import numpy as np

# Generates a swift IC file for the Kelvin-Helmholtz vortex in a periodic box

# Parameters
gamma = 5./3.     # Gas adiabatic index
B0    = 1./(np.sqrt(4.*np.pi)) # Bfield
P0    = gamma * pow(B0,2)   # Pressure
rho0  = gamma * P0             # Density
N1D   = 64        # Number of particles in one dimension

fileOutputName = "OrszagTangVortex.hdf5"

#---------------------------------------------------

N = N1D ** 3
x = np.linspace(0., 1., N1D + 1)
x = 0.5 * (x[1:] + x[:-1])
y = x
z = x
xx, yy, zz = np.meshgrid(x, y, z)
pos = np.zeros((N, 3))
pos[:, 0] = xx.reshape((N))
pos[:, 1] = yy.reshape((N))
pos[:, 2] = zz.reshape((N))
h = np.ones(N) * 2. / N1D

vol = 1.

# Generate extra arrays
v = np.zeros((N, 3))
b = np.zeros((N, 3))
epa = np.zeros(N)
epb = np.zeros(N)
ids = np.linspace(1, N, N)
m = np.ones(N) * rho0 * vol / N
u = np.ones(N) * P0 / (rho0 * (gamma - 1.))


v[:, 0] =  -np.sin(2. * np.pi * pos[:, 1]) 
v[:, 1] =   np.sin(2. * np.pi * pos[:, 0]) 
v[:, 2] =  0. 
b[:, 0] =  -B0 * np.sin(2. * np.pi * pos[:, 1]) 
b[:, 1] =   B0 * np.sin(4. * np.pi * pos[:, 0]) 
b[:, 2] =  0. 
epa[:]  = B0/(2.*np.pi)*(np.cos(2.*np.pi*pos[:,1]) + np.cos(4.*np.pi*pos[:,0]))
epb[:]  = pos[:,2]
            
#File
fileOutput = h5py.File(fileOutputName, 'w')

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [1., 1., 1.]
grp.attrs["NumPart_Total"] =  [N, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

#Units
grp = fileOutput.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.
grp.attrs["Unit mass in cgs (U_M)"] = 1.
grp.attrs["Unit time in cgs (U_t)"] = 1.
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

#Particle group
grp = fileOutput.create_group("/PartType0")
grp.create_dataset("Coordinates", data = pos, dtype = 'd')
grp.create_dataset("Velocities", data = v, dtype = 'f')
grp.create_dataset("Masses", data = m, dtype = 'f')
grp.create_dataset("SmoothingLength", data = h, dtype = 'f')
grp.create_dataset("InternalEnergy", data = u, dtype = 'f')
grp.create_dataset("ParticleIDs", data = ids, dtype = 'L')
grp.create_dataset("Bfield", data = b, dtype = 'f')
grp.create_dataset("EPalpha", data = epa, dtype = 'f')
grp.create_dataset("EPbeta" , data = epb, dtype = 'f')

fileOutput.close()
