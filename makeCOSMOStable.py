import esutil
import glob
from myutils import fitsstack
from astropy.io import fits
from astropy.table import Table, join

#cosmos15 = fits.open('/Users/Christina/COSMOS/COSMOS2015_Laigle+_v1.1.fits')[1].data

y1a1_files = glob.glob('/Users/Christina/DES/data/y1a1_gold_dfull_00000*.fits')

y1a1_hdu = fitsstack(y1a1_files)
pri_hdu = fits.Header()
hdu_list = fits.HDUList([pri_hdu, y1a1_hdu])
hdu_list.writeto('y1a1_gold_dfull.fits')

print y1a1_hdu.columns

plt.scatter(y1a1_hdu['ra'], y1a1_hdu['dec'], s=4.)
plt.show()
