import h5py
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np


def getWF(filename,slice=1):
    hfl = h5py.File(filename,'r')
    slc = 'slice%6.6d' % slice
    ng = hfl['gridpoints'][()][0]
    dg = hfl['gridsize'][()][0]
    fre = hfl[slc]['field-real'][()]
    fim = hfl[slc]['field-imag'][()]
    inten = np.reshape(fre*fre+fim*fim, (ng,ng))
    return inten,dg*(ng-1)*0.5*1e3


def getPS(filename,slice=1):
    hfl = h5py.File(filename,'r')
    slc = 'slice%6.6d' % slice
    gamma = hfl[slc]['gamma'][()]*0.511e-3
    theta = np.mod(hfl[slc]['theta'][()]-np.pi*0.5,2*np.pi)
    return theta, gamma

def getTS(filename,slice=1):
    hfl = h5py.File(filename,'r')
    slc = 'slice%6.6d' % slice
    x = hfl[slc]['x'][()]*1e6
    theta = np.mod(hfl[slc]['theta'][()]-np.pi*0.5,2*np.pi)
    return theta, x

# plot wavefront
istep = 184
figs, axs = plt.subplots(2, 2)
color = 'yellow'
for i1 in range(2):
    for i2 in range(2):
        i = (i2*2+i1 +1)*istep
        inten, dg = getWF('Example2.%d.fld.h5' % i, 1)
        axs[i2, i1].imshow(inten,extent=(-dg,dg,-dg,dg))
        txt = 'z = %3.1f m' % (9.5*(i2*2+i1+1))
        axs[i2, i1].text(-0.15, 0.15, txt ,color = color)

axs[1,0].set_xlabel(r'$x$ (mm)')
axs[1,1].set_xlabel(r'$x$ (mm)')
axs[0,0].set_ylabel(r'$y$ (mm)')
axs[1,0].set_ylabel(r'$y$ (mm)')

plt.show()

# get range for phase space plots
hid = h5py.File('Example2.out.h5','r')
emin = np.min(hid['Beam']['emin'][()])*0.511e-3
emax = np.max(hid['Beam']['emax'][()])*0.511e-3
xmin = np.min(hid['Beam']['xmin'][()])*1e6
xmax = np.max(hid['Beam']['xmax'][()])*1e6
hid.close()


# plot final phase space
t,g = getPS('Example2.700.par.h5',1)
plt.scatter(t,g,s=0.2)
plt.xlabel(r'$\theta$ (rad)')
plt.ylabel(r'$E$ (GeV)')
plt.xlim([0,2*np.pi])
plt.ylim([emin,emax])
plt.show()


# Animate phase space

fig = plt.figure()
ax = plt.axes(xlim=(0,2*np.pi),ylim=(emin,emax))
ax.set_xlabel(r'$\theta$ (rad)')
ax.set_ylabel(r'$E$ (GeV)')
scat = ax.scatter([],[],s=0.2)


def init():
    scat.set_offsets([])
    return scat,

def animate(i):
    file = 'Example2.%d.par.h5' % (2*i)
    print(file)
    t, g = getPS(file,1)
    scat.set_offsets(np.hstack((t[:,np.newaxis],g[:,np.newaxis])))
    return scat,

anim = animation.FuncAnimation(fig,animate,init_func=init,blit = False,interval=20,frames = 500)
anim.save('Animation1.mp4')

fig = plt.figure()
ax = plt.axes(xlim=(0,2*np.pi),ylim=(xmin,xmax))
ax.set_xlabel(r'$\theta$ (rad)')
ax.set_ylabel(r'$x$ ($\mu$m)')
scat = ax.scatter([],[],s=0.2)


def init():
    scat.set_offsets([])
    return scat,

def animate(i):
    file = 'Example2.%d.par.h5' % (2*i)
    print(file)
    t, g = getTS(file,1)
    scat.set_offsets(np.hstack((t[:,np.newaxis],g[:,np.newaxis])))
    return scat,

anim = animation.FuncAnimation(fig,animate,init_func=init,blit = False,interval=20,frames = 500)
anim.save('Animation2.mp4')