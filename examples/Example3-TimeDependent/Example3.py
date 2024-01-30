import h5py
import matplotlib.pyplot as plt
import numpy as np

hid = h5py.File('Example3.out.h5','r')




# plot beam parameters - current profile and energy
s = hid['Global']['s'][()]/3e8*1e15
current = hid['Beam']['current'][()][0,:]*1e-3
energy = hid['Beam']['energy'][()][0,:]*0.511*1e-3


fig,ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel(r'$s$ (fs)')
ax1.set_ylabel(r'$I$ (kA)',color = color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.plot(s,current,color = color)

ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel(r'$E$ (GeV)',color = color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.plot(s,energy,color = color)
plt.show()





# plot global parameters
z = hid['Lattice']['zplot'][()]
energy = hid['Field']['Global']['energy'][()]*1e6
energy0=energy[220]
farfield = hid['Field']['Global']['intensity-farfield'][()]*1e-21

fig,ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel(r'$z$ (m)')
ax1.set_ylabel(r'$E$ ($\mu$J)',color = color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.semilogy(z,energy,color = color)
ax1.set_ylim([1e-3,1e3])

ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel(r'$<I(\theta=0)>$ (GW/$\mu$rad$^2$)',color = color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.semilogy(z,farfield,color = color)
ax2.set_ylim([1e-2,1e5])
plt.show()

# power profile
power1 = hid['Field']['power'][()][50,:]*1e-9
farfield1 = hid['Field']['intensity-farfield'][()][50,:]*1e-24

power2 = hid['Field']['power'][()][220,:]*1e-9
farfield2 = hid['Field']['intensity-farfield'][()][220,:]*1e-24

fig,ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel(r'$t$ (fs)')
ax1.set_ylabel(r'$P$ (MW)',color = color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.plot(s,power1*1e3,color = color)
#ax1.set_ylim([1e-3,1e3])

ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel(r'$I(\theta=0)$ (GW/$\mu$rad$^2$)',color = color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.plot(s,farfield1*1e3,color = color)
#ax2.set_ylim([1e-2,1e5])
plt.show()

fig,ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel(r'$t$ fs)')
ax1.set_ylabel(r'$P$ (GW)',color = color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.plot(s,power2,color = color)
#ax1.set_ylim([1e-3,1e3])

ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel(r'$I(\theta=0)$ (TW/$\mu$rad$^2$)',color = color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.plot(s,farfield2,color = color)
#ax2.set_ylim([1e-2,1e5])
plt.show()

fig,ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel(r'$t$ fs)')
ax1.set_ylabel(r'$P$ (GW)',color = color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.plot(s,power2,color = color)
ax1.set_xlim([37,42])
#ax1.set_ylim([1e-3,1e3])

ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel(r'$I(\theta=0)$ (TW/$\mu$rad$^2$)',color = color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.plot(s,farfield2,color = color)
ax2.set_xlim([37,42])
#ax2.set_ylim([1e-2,1e5])
plt.show()

power =  hid['Field']['power'][()]
pmean= np.mean(power, axis = 1)
zmean =np.mean(power, axis=1)
for i in range(len(pmean)):
    if pmean[i] == 0:
        pmean[i]=1.
    power[i,:]*=1./pmean[i]
plt.imshow(np.flipud(power), aspect='auto', interpolation='none', extent=(np.min(s),np.max(s),np.min(z),np.max(z)))
plt.xlabel(r'$t$ (fs)')
plt.ylabel(r'$z$ (m)')
plt.show()

# plot spectrum

# plot sample spectrum

sig = hid['Field']['intensity-farfield'][()][220,:]
phi = hid['Field']['phase-farfield'][()][220,:] 
freq = hid['Global']['frequency'][()]

signal = np.sqrt(sig)*np.exp(1j*phi)
spec = np.abs(np.fft.fftshift(np.fft.fft(signal)))**2
norm = energy0/np.sum(spec)/(freq[1]-freq[0])
spec =norm*spec

plt.plot(freq*1e-3,spec)
plt.xlim([12.400,12.800])
plt.xlabel(r'$E_{ph}$ (keV)')
plt.ylabel(r'$P(E_{ph})$ ($\mu$J/eV)')
plt.show()

exit()
                   


hid.close()
