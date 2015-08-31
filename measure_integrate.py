tory #!/usr/bin/env python
import os
import subprocess
import glob
import numpy

#qsub -pe omp 8 -l h_vmem=300M -t 1-5:1 -o pixtestoutputs/ -e pixtestoutputs/ ./measure_integrate.py

task_id = str(os.environ['SGE_TASK_ID']).zfill(3)
#task_id = '_'
os.environ['OMP_NUM_THREADS']="12"

n_batch = 1

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
g1_targ = 0.02
g2_targ = 0.0

#template specific
n_files_temp = 1
nx_tiles_temp = 100
ny_tiles_temp = 100
g1_temp = 0.0
g2_temp = 0.0

#greatMeasure parameters
weight_sigma = 3.5

#greatIntegrate parameters
n_sample = 20000
select_sn_min = 8
select_sn_max = 20
prior_sigma_cutoff = 6.

#absolute paths to directories
home_dir = os.path.join('/home','ckrawiec')
sim_dir = os.path.join(home_dir,'sims')
image_dir = os.path.join(sim_dir,'images')
measure_dir = os.path.join(home_dir,'git','bfd')
integrate_dir = os.path.join(home_dir,'git','bfd')
moment_dir = os.path.join(home_dir,'moments')

#absolute paths to files
image_script = os.path.join(sim_dir, 'stamps2.py')
measure_prog = os.path.join(measure_dir, 'greatMeasure')
integrate_prog = os.path.join(integrate_dir, 'greatIntegrate')
constant_prog = os.path.join(integrate_dir, 'greatConstant')
image_files_temp = os.path.join(image_dir, 'template_stamps_{0}_*.fits'.format(task_id))
psf_file_temp = os.path.join(image_dir, 'template_stamps_{0}_01.psf'.format(task_id))
moment_files = os.path.join(moment_dir, 'stamps_moments_{0}_*.fits'.format(task_id))
nominal_moment_file = os.path.join(moment_dir, 'stamps_moments_{0}_01.fits'.format(task_id))



#make template images
subprocess.call(['python',image_script,
                 '-id',task_id,
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




for i in range(n_batch):
    #names based on job task id
    id = task_id+'_'+str(i+1).zfill(2)

    image_files_targ = os.path.join(image_dir, 'stamps_{0}_*.fits'.format(id))
    psf_file_targ = os.path.join(image_dir, 'stamps_{0}_01.psf'.format(id))
    psf_files_targ = os.path.join(image_dir, 'stamps_{0}_*.psf'.format(id))

    moment_file = os.path.join(moment_dir,
                               'stamps_moments_{0}.fits'.format(id))


    #make target images and save variance of added noise
    noise_var = subprocess.check_output(['python',image_script,
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

    #parameters for greatMeasure
    measure_params = ['-targetFiles={0}'.format(image_files_targ),
                      '-targetPsf={0}'.format(psf_file_targ),
                      '-outfile={0}'.format(moment_file),
                      '-targetNoise={0}'.format(noise_var),
                      '-targetStampSize={0}'.format(stamp_size),
                      '-weightKSigma={0}'.format(weight_sigma)]
    

    #measure moments of images using greatMeasure
    subprocess.call([measure_prog]+measure_params)

    #delete images
    for image in glob.glob(image_files_targ):
        os.remove(image)
    for psf in glob.glob(psf_files_targ):
        os.remove(psf)

#parameters for greatIntegrate
integrate_params = ['-templateFiles={0}'.format(image_files_temp),
                    '-templatePsfs={0}'.format(psf_file_temp),
                    '-momentFiles={0}'.format(moment_files),
                    '-templateStampSize={0}'.format(stamp_size),
                    '-nominalMomentFile={0}'.format(nominal_moment_file),
                    '-nSample={0}'.format(n_sample),
                    '-selectSnMin={0}'.format(select_sn_min),
                    '-selectSnMax={0}'.format(select_sn_max),
                    '-weightKSigma={0}'.format(weight_sigma),
                    '-priorSigmaCutoff={0}'.format(prior_sigma_cutoff)]

#call greatIntegrate
subprocess.call([integrate_prog]+integrate_params)

#call greatConstant
bfile = open(os.path.join('/home','ckrawiec','moments','copies',
                          'constant_before_deselect_{0}.txt'.format(task_id)),'w')

for mb in glob.glob(moment_files):
    out = subprocess.check_output([constant_prog, mb])
    print>>bfile, out

bfile.close()


#parameters for greatIntegrate with Deselection
integrate_params_deselect = ['-templateFiles={0}'.format(image_files_temp),
                             '-templatePsfs={0}'.format(psf_file_temp),
                             '-momentFiles={0}'.format(moment_files),
                             '-templateStampSize={0}'.format(stamp_size),
                             '-nominalMomentFile={0}'.format(nominal_moment_file),
                             '-nSample={0}'.format(n_sample),
                             '-selectSnMin={0}'.format(select_sn_min),
                             '-selectSnMax={0}'.format(select_sn_max),
                             '-weightKSigma={0}'.format(weight_sigma),
                             '-priorSigmaCutoff={0}'.format(prior_sigma_cutoff),
                             '-doDeselection=T']

#call greatIntegrate with Deselection
subprocess.call([integrate_prog]+integrate_params_deselect)

#call greatConstant again
afile = open(os.path.join('/home','ckrawiec','moments','copies',
                          'constant_after_deselect_{0}.txt'.format(task_id)),'w')

for ma in glob.glob(moment_files):
    out = subprocess.check_output([constant_prog, ma])
    print>>afile, out

afile.close()

#delete template images
#for image in glob.glob(image_files_temp):
#    os.remove(image)
#os.remove(psf_file_temp)
    
