#!/usr/bin/env python
import os
import subprocess
import glob
import numpy
import socket

#hostname
print socket.gethostname()

#qsub jobs
task_id = str(os.environ['SGE_TASK_ID'])

#galsim script parameters
stamp_size = 48
pixel_scale = 1.
gal_ellip_sigma = 0.2
psf_hlr = 1.5
psf_e2 = 0.05
gal_snr_min = 5
gal_flux_min = 1000.
gal_flux_max = 5000.
gal_hlr_min = 1.*psf_hlr
gal_hlr_max = 2.*psf_hlr

#target specific
n_files_targ = 50
nx_tiles_targ = 100
ny_tiles_targ = 100
g1_targ = 0.06
g2_targ = 0.0

#template specific
n_files_temp = 5
nx_tiles_temp = 50
ny_tiles_temp = 100
g1_temp = 0.0
g2_temp = 0.0

#absolute paths to directories
home_dir = os.path.join('/home','ckrawiec')
sim_dir = os.path.join(home_dir,'sims')
image_dir = os.path.join('/data3','garyb','BFD','sims','images','g06')

#names based on job task id
id = task_id.zfill(4)

#absolute paths to files
image_script = os.path.join(sim_dir, 'stamps_final.py')
image_files_temp = os.path.join(image_dir, 'template_stamps_{0}_*.fits'.format(id))
psf_file_temp = os.path.join(image_dir, 'template_stamps_{0}_01.psf'.format(id))
image_files_targ = os.path.join(image_dir, 'stamps_{0}_*.fits'.format(id))
psf_file_targ = os.path.join(image_dir, 'stamps_{0}_01.psf'.format(id))
psf_files_targ = os.path.join(image_dir, 'stamps_{0}_*.psf'.format(id))

#make template images
subprocess.call(['python',image_script,
                 '-file_dir',image_dir,
                 '-id',id,
                 '-n_files', str(n_files_temp),
                 '-nx_tiles', str(nx_tiles_temp),
                 '-ny_tiles', str(ny_tiles_temp),
                 '-without_noise',
                 '-stamp_size', str(stamp_size),
                 '-pixel_scale', str(pixel_scale),
                 '-psf_hlr', str(psf_hlr),
                 '-psf_e2', str(psf_e2),
                 '-g1', str(g1_temp),
                 '-g2', str(g2_temp),
                 '-gal_ellip_sigma', str(gal_ellip_sigma),
                 '-gal_flux_min', str(gal_flux_min),
                 '-gal_flux_max', str(gal_flux_max),
                 '-gal_hlr_min', str(gal_hlr_min),
                 '-gal_hlr_max', str(gal_hlr_max)])

#make target images and save variance of added noise (in text file)
subprocess.call(['python',image_script,
                 '-file_dir',image_dir,
                 '-id',id,
                 '-n_files', str(n_files_targ),
                 '-nx_tiles', str(nx_tiles_targ),
                 '-ny_tiles', str(ny_tiles_targ),
                 '-stamp_size', str(stamp_size),
                 '-pixel_scale', str(pixel_scale),
                 '-with_noise',
                 '-psf_hlr', str(psf_hlr),
                 '-psf_e2', str(psf_e2),
                 '-g1', str(g1_targ),
                 '-g2', str(g2_targ),
                 '-gal_ellip_sigma', str(gal_ellip_sigma),
                 '-gal_snr_min', str(gal_snr_min),
                 '-gal_flux_min', str(gal_flux_min),
                 '-gal_flux_max', str(gal_flux_max),
                 '-gal_hlr_min', str(gal_hlr_min),
                 '-gal_hlr_max', str(gal_hlr_max)])
    
