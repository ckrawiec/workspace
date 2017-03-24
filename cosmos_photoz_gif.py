import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib.animation as anim

class Zcolors:

    def __init__(self, photoz_file, zbin_size):
        self.ztable = Table.read(photoz_file)
        
        self.g = self.ztable['MAG_DETMODEL_G']
        self.r = self.ztable['MAG_DETMODEL_R']
        self.i = self.ztable['MAG_DETMODEL_I']
        self.z = self.ztable['MAG_DETMODEL_Z']
        
        self.photoz = self.ztable['zminchi2']
        
        self.zmin = 0.0
        self.zmax = self.zmin+zbin_size

        self.zmask = (self.photoz > self.zmin) & (self.photoz < self.zmax)

    def colors(self):
        ri = self.r[self.zmask]-self.i[self.zmask]
        gr = self.g[self.zmask]-self.r[self.zmask]
        
        return (ri, gr)

    def step(self, dz):
        self.zmask = (self.photoz > self.zmin) & (self.photoz < self.zmax)
        self.zmin += dz
        self.zmax += dz

filename = '/Users/Christina/DES/data/sva1_gold_detmodel_MC1_good_regions_cosmos.fits'
save_file = '/Users/Christina/DES/data/sva1_gold_detmodel_MC1_good_regions_cosmos_zcolors.mp4'
zbin_width = 0.2

zcolors = Zcolors(filename, zbin_width)
dz = 0.2

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-1,4), ylim=(-1,4))

ax.grid()

line, = ax.plot([], [], 'o', markeredgecolor='none')
z_text = ax.text(2.2, 3.5, '', fontsize=20)

def animate(i):
    global zcolors, dz
    zcolors.step(dz)

    line.set_data(*zcolors.colors())
    z_text.set_text('z = {}-{}'.format(zcolors.zmin, zcolors.zmax))
    return line, z_text

ani = anim.FuncAnimation(fig, animate, frames=25, blit=True)

ani.save(save_file)
