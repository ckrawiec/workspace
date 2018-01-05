import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord

inputfile = '/Users/Christina/DES/data/Galaxia_SV390deg_desflux.fits'
outputfile = '/Users/Christina/DES/data/Galaxia_SV390deg_desflux_radec.fits'

tab = Table.read(inputfile)
galc = SkyCoord(frame="galactic", l=tab['glon'], b=tab['glat'], unit=u.deg)
galrd = galc.icrs
tab['RA'] = galrd.ra.degree
tab['DEC'] = galrd.dec.degree

tab.write(outputfile)
