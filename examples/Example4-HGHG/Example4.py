import h5py
import matplotlib.pyplot as plt
import numpy as np


#sig = 0.1
#dE=7*sig
#lam0=53e-9
#lamr=265e-9
#kr = 2*np.pi/lamr
#N=100000
#theta=np.random.uniform(0,2*np.pi,N)
#gamma=np.random.normal(0,sig,N)

#gamma+= dE*np.sin(theta)


#Nscan=100
#gamma0=1500
#R56=np.linspace(0,3e-4,num=Nscan)
#b = np.zeros((Nscan))
#for i,r in enumerate(R56):
#    t0 = theta + kr*r*gamma/gamma0
#    b[i]=np.abs(np.mean(np.exp(-1j*5*t0)))
#plt.plot(R56,b)
#plt.show()
#exit()

def plotPS(filename):
    hid = h5py.File(filename,'r')
    idx0=450
    dI = 10
    for i in range(dI):
        slice='slice%6.6d' % (idx0+i)
        theta = np.mod(hid[slice]['theta'][()][::5],2*np.pi)+i*2*np.pi
        gamma = hid[slice]['gamma'][()][::5]-1500.
        plt.scatter(theta,gamma,s=0.5)
    plt.xlabel(r'$k_r\cdot \Delta s$')
    plt.ylabel(r'$\Delta\gamma$')
    hid.close()
    plt.show()

plotPS('dump_a.par.h5')
plotPS('dump_b.par.h5')

hida = h5py.File('Example4_a.out.h5','r')
hidb = h5py.File('Example4_b.out.h5','r')
buncha = hida['Beam']['bunching'][()][0,:]
bunchb = hidb['Beam']['bunching'][()][0,:]
bphia = hida['Beam']['bunchingphase'][()][0,:]
bphib = hidb['Beam']['bunchingphase'][()][0,:]
cura = hida['Beam']['current'][()][0,:]
curb = hidb['Beam']['current'][()][0,:]
kr = 2*np.pi/59e-9
sa = hida['Global']['s'][()]*kr-451*2*np.pi
sb = hidb['Global']['s'][()]*kr-451*2*np.pi

plt.plot(sa,buncha,label='Beamlets')
plt.plot(sb,bunchb,label='One-4-One')
plt.legend()
plt.xlabel(r'$k_r\cdot \Delta s$')
plt.ylabel(r'|<exp(i$\theta$)>|')
plt.xlim((0,20*np.pi))
plt.show()

plt.plot(sa,cura,label='Beamlets')
plt.plot(sb,curb,label='One-4-One')
plt.legend()
plt.xlabel(r'$k_r\cdot \Delta s$')
plt.ylabel(r'$I$ (A)')
plt.xlim((0,20*np.pi))
plt.show()


za = hida['Lattice']['zplot'][()]
zb = hidb['Lattice']['zplot'][()]
powera = np.mean(hida['Field']['power'][()],axis=1)*1e-6
powerb = np.mean(hidb['Field']['power'][()],axis=1)*1e-6


ns = len(sa)
f0 = 1240./59.
f = np.linspace(0,f0,num=ns)+0.5*f0


plt.plot(za,powera,label='Beamlets')
plt.plot(zb,powerb,label='One-4-One')
plt.legend()
plt.xlabel(r'$z$ (m)')
plt.ylabel(r'<P> (MW)')
plt.show()


siga=np.abs(np.fft.fftshift(np.fft.fft(buncha*cura*np.exp(1j*bphia))))
sigb=np.abs(np.fft.fftshift(np.fft.fft(bunchb*curb*np.exp(1j*bphib))))

plt.plot(f,siga*1e-4,label='Beamlets')
plt.plot(f,sigb*1e-4,label='One-4-One')
plt.legend()
plt.xlabel(r'$E_{ph}$ (eV)')
plt.ylabel(r'$P_b(E_{ph})$ (arb. units)')
plt.show()

nz=220




powera =  hida['Field']['power'][()]
powerb =  hidb['Field']['power'][()]

plt.plot(sa,powera[nz,:]*1e-6,label='Beamlets')
plt.plot(sb,powerb[nz,:]*1e-6,label='One-4-One')
plt.legend()
plt.xlabel(r'$k_r\cdot \Delta s$')
plt.ylabel(r'$P$ (MW)')
plt.show()



pmean= np.max(powera, axis = 1)
for i in range(len(pmean)):
    if pmean[i] == 0:
        pmean[i]=1.
    powera[i,:]*=1./pmean[i]
plt.imshow(np.flipud(powera), aspect='auto', interpolation='none', extent=(np.min(sa),np.max(sa),np.min(za),np.max(za)))
plt.xlabel(r'$k_r\cdot \Delta s$')
plt.ylabel(r'$z$ (m)')
plt.show()




powerb =  hidb['Field']['power'][()]
pmean= np.max(powerb, axis = 1)
for i in range(len(pmean)):
    if pmean[i] == 0:
        pmean[i]=1.
    powerb[i,:]*=1./pmean[i]
plt.imshow(np.flipud(powerb), aspect='auto', interpolation='none', extent=(np.min(sa),np.max(sa),np.min(za),np.max(za)))
plt.xlabel(r'$k_r\cdot \Delta s$')
plt.ylabel(r'$z$ (m)')
plt.show()

