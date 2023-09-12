import h5py
import matplotlib.pyplot as plt

hid = h5py.File('Example1.out.h5','r')

# plot the lattice

z = hid['Lattice']['z'][()]
aw = hid['Lattice']['aw'][()]
qf = hid['Lattice']['qf'][()]

fig,ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel(r'$z$ (m)')
ax1.set_ylabel(r'$a_w$',color = color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.step(z,aw,color = color, where='post')

ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel(r'$k_1$ (m$^{-2}$)',color = color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.step(z,qf,color = color,where='post')
plt.show()

# plot the beam sizes
z = hid['Lattice']['zplot'][()]
bx = hid['Beam']['xsize'][()]
by = hid['Beam']['ysize'][()]
fx = hid['Field']['xsize'][()]
fy = hid['Field']['ysize'][()]
plt.plot(z,bx*1e6,label=r'Beam: $\sigma_x$')
plt.plot(z,by*1e6,label=r'Beam: $\sigma_y$')
plt.plot(z,fx*1e6,label=r'Field: $\sigma_x$')
plt.plot(z,fy*1e6,label=r'Field: $\sigma_y$')
plt.legend()
plt.xlabel(r'$z$ (m)')
plt.ylabel(r'$\sigma_{x,y}$ ($\mu$m)')
plt.ylim([0,60])
plt.show()

# plot power and bunching
z = hid['Lattice']['zplot'][()]
b = hid['Beam']['bunching'][()]
p = hid['Field']['power'][()]

fig,ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel(r'$z$ (m)')
ax1.set_ylabel(r'$P$ (W)',color = color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.semilogy(z,p,color = color)

ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel(r'$<\exp(i\theta)>$',color = color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.semilogy(z,b,color = color)
ax2.set_ylim([1e-3,0.5])
plt.show()

hid.close()