import galsim
import os
import numpy as np
import math
import argparse
from astropy.io import fits

#Moffat PSF with ellipticity e2
#Galaxies have bulge and disk components
#    -exponential disk
#    -deVaucouleurs profile
#    -disk and bulge have same e
#    -uniform distribution for bulge fraction
#    -different centroids chosen from circle of radius r_e


def get_params():
    parser = argparse.ArgumentParser(description='Parameters for galaxy postage stamps.')

    #File structure
    parser.add_argument('-n_files', default=1, type=int, action='store',
                        dest='n_files',
                        help='The number of separate files to produce (default=1)' )
    parser.add_argument('-id', default='', type=str, action='store',
                        dest='id',
                        help='Output image files will be of the form stamps_<id>_01.fits (default='')')
    parser.add_argument('-nx_tiles', default=1, type=int, action='store',
                        dest='nx_tiles',
                        help='The number of postage stamps in the x-direction (default=10)')
    parser.add_argument('-ny_tiles', default=1, type=int, action='store',
                        dest='ny_tiles',
                        help='The number of postage stamps in the y-direction (default=10)')
    parser.add_argument('-stamp_size', default=48, type=int, action='store',
                        dest='stamp_size',
                        help='Size of stamp edge in pixels (default=48)')
    parser.add_argument('-pixel_scale', default=1., type=float, action='store',
                        dest='pixel_scale',
                        help='Pixel scale of the images in arcsec/pixel (default=1.)')
    parser.add_argument('-random_seed', default=0, type=int, action='store',
                        dest='random_seed',
                        help='Integer seed for the random number generator (default=0)')

    #Added noise
    parser.add_argument('-without_noise', action='store_true',
                        dest='noise_suppression',
                        help='Include for galaxies without added noise')
    parser.add_argument('-with_noise', action='store_false',
                        dest='noise_suppression',
                        help='Include for galaxies with added noise (this is the default)')
    parser.set_defaults(noise_suppression=False)
    


    #PSF
    parser.add_argument('-psf_hlr', default=1.5, type=float, action='store',
                        dest='psf_hlr',
                        help='Half-light radius of PSF in arcsec (default=1.5)')
    parser.add_argument('-psf_e2', default=0.05, type=float, action='store',
                        dest='psf_e2',
                        help='PSF e2 (default=0.05)')

    #Galaxy
    parser.add_argument('-gal_ellip_sigma', default=0.0, type=float, action='store',
                        dest='gal_ellip_sigma',
                        help='Sigma of galaxy ellipticity distribution (default=0.0)')
    parser.add_argument('-gal_flux_min', default=100., type=float, action='store',
                        dest='gal_flux_min',
                        help='Minimum galaxy flux (default=100.0)')
    parser.add_argument('-gal_flux_max', default=100., type=float, action='store',
                        dest='gal_flux_max',
                        help='Maximum galaxy flux (default=100.0)')
    parser.add_argument('-gal_snr_min', default=10., type=float, action='store',
                        dest='gal_snr_min',
                        help='Minimum S/N (default=10.0)')
    parser.add_argument('-gal_hlr_min', default=3., type=float, action='store',
                        dest='gal_hlr_min',
                        help='Minimum half-light radius (default=3.)')
    parser.add_argument('-gal_hlr_max', default=3., type=float, action='store',
                        dest='gal_hlr_max',
                        help='Maximum half-light radius (default=3.)')

    #Added shear
    parser.add_argument('-g1', default=0., type=float, action='store',
                        dest='g1',
                        help='Shear g1 (default=0.)')
    parser.add_argument('-g2', default=0., type=float, action='store',
                        dest='g2',
                        help='Shear g2 (default=0.)')

    params = parser.parse_args()
    return params

class BobsEDist:
    """
    Sets up an ellipticity distribution: exp^(-e^2/2sig^2)*(1-e^2)^2
    """
    
    def __init__(self, sigma, rng):
        self.ellip_sigma = sigma
        self.ud = rng
        self.gd = galsim.GaussianDeviate(self.ud, sigma=self.ellip_sigma)
    
    def sample(self):
        if self.ellip_sigma==0.:
            e1,e2 = 0.,0.
        else:
            while True:
                while True:
                    e1 = self.gd()
                    e2 = self.gd()
                    esq = e1*e1 + e2*e2
                    if esq < 1.:
                        break
                if self.ud() < (1-esq)*(1-esq) :
                    break
        return e1,e2
    
def main():

    p = get_params()

    rng = galsim.UniformDeviate(p.random_seed)
 
    #where to save images
    image_dir = os.path.join('/home','ckrawiec','sims','images')
    
    target_file_name = "stamps_{0}_.fits".format(p.id)
    template_file_name = "template_stamps_{0}_.fits".format(p.id)

    if p.noise_suppression==False:
        gal_name = os.path.splitext(target_file_name)
    else:
        gal_name = os.path.splitext(template_file_name)
        
    #initial noise
    noise = galsim.GaussianNoise(rng)

    #create test image to determine noise level
    if p.noise_suppression==False:
        test_psf = galsim.Moffat(3.5,half_light_radius=p.psf_hlr*p.pixel_scale)
        test_psf = test_psf.shear(e2=p.psf_e2)

        test_bulge = galsim.DeVaucouleurs(half_light_radius=((p.gal_hlr_max-p.gal_hlr_min)/2.
                                          +p.gal_hlr_min)*p.pixel_scale)
        test_disk = galsim.Exponential(half_light_radius=((p.gal_hlr_max-p.gal_hlr_min)/2.
                                       +p.gal_hlr_min)*p.pixel_scale)
        test_gal = 0.5*test_bulge + 0.5*test_disk
        test_gal = test_gal.withFlux(p.gal_flux_min)

        test_image = galsim.ImageF(p.stamp_size,
                                   p.stamp_size,
                                   scale=p.pixel_scale)

        test_final = galsim.Convolve([test_psf,test_gal])
        test_final.drawImage(test_image)

        noise_var = test_image.addNoiseSNR(noise,
                                           p.gal_snr_min,
                                           preserve_flux=True)
        print noise_var

        #final noise
        noise = galsim.GaussianNoise(rng, sigma=np.sqrt(noise_var))
            
    #make galaxy image
    gal_image = galsim.ImageF(p.nx_tiles*p.stamp_size,
                              p.ny_tiles*p.stamp_size,
                              scale=p.pixel_scale)

    #make ellipticity distribution from sigma
    ellip_dist = BobsEDist(p.gal_ellip_sigma, rng)

    for ifile in range(p.n_files):
        #all files have same psf
        psf = galsim.Moffat(3.5,half_light_radius=p.psf_hlr*p.pixel_scale)
        psf = psf.shear(e2=p.psf_e2)
        psf_image = galsim.ImageF(p.stamp_size, p.stamp_size, scale=p.pixel_scale)
        psf.drawImage(psf_image)
                
        for iy in range(p.ny_tiles):
            for ix in range(p.nx_tiles):
                
                b = galsim.BoundsI(ix*p.stamp_size+1, (ix+1)*p.stamp_size,
                                   iy*p.stamp_size+1, (iy+1)*p.stamp_size)
                sub_gal_image = gal_image[b]
                    
                hlr =  (rng() * (p.gal_hlr_max-p.gal_hlr_min) + p.gal_hlr_min) * p.pixel_scale
                flux = rng() * (p.gal_flux_max-p.gal_flux_min) + p.gal_flux_min
                bulge_frac = rng()

                e1,e2 = ellip_dist.sample()
    
                #galaxy = ExponentialDisk + DeVaucouleurs
                this_bulge = galsim.DeVaucouleurs(half_light_radius = hlr)
                this_disk = galsim.Exponential(half_light_radius = hlr)

                #same ellipticity
                this_bulge = this_bulge.shear(e1=e1,e2=e2)
                this_disk = this_disk.shear(e1=e1,e2=e2)

                #Apply a shift in disk center within one pixel
                shift_r = p.pixel_scale/2.
                dx = shift_r * (rng()*2.-1.)
                dy = shift_r * (rng()*2.-1.)
                this_disk = this_disk.shift(dx,dy)
                
                #Apply a shift in bulge center within circle of radius hlr
                theta = rng() * 2. * np.pi
                r = rng() * hlr
                dx = r * np.cos(theta)
                dy = r * np.sin(theta)
                this_bulge = this_bulge.shift(dx,dy)

                this_gal = bulge_frac*this_bulge + (1-bulge_frac)*this_disk
                this_gal = this_gal.withFlux(flux)
                
                #All galaxies have same applied shear
                this_gal = this_gal.shear(g1=p.g1,g2=p.g2)
                
                final = galsim.Convolve([psf,this_gal])

                final.drawImage(sub_gal_image)            

                if p.noise_suppression==False:
                    sub_gal_image.addNoise(noise)

                
        gal_file = os.path.join(image_dir,''.join([gal_name[0],
                                                str(ifile+1).zfill(2),
                                                gal_name[1]]))
        psf_file = os.path.join(image_dir,''.join([gal_name[0],
                                                str(ifile+1).zfill(2),
                                                '.psf']))
        gal_image.write(gal_file)
        psf_image.write(psf_file)


if __name__ == "__main__":
    main()


