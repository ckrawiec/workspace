#!/usr/bin/env python
import os
import subprocess
import glob
import numpy

#qsub -pe omp 12 -l h_vmem=2G -t 1-1:1 -o mptest1.out -e mptest1.err ./measure_integrate_constant.py

#qsub jobs
task_id = str(os.environ['SGE_TASK_ID'])
os.environ['OMP_NUM_THREADS']="12"

#file names based on job task id
id = task_id.zfill(4)

#galsim script parameters
stamp_size = 48
pixel_scale = 1.

#greatMeasure parameters
weight_sigma = 3.5

#greatIntegrate parameters
n_sample = 50000
select_sn_min = 8
select_sn_max = 20
prior_sigma_step = 1.1
prior_sigma_cutoff = 6.

#absolute paths to directories
home_dir = os.path.join('/home','ckrawiec')
data3_dir = os.path.join('/data3','garyb','BFD')
image_dir = os.path.join(data3_dir,'sims','images','g06')
measure_dir = os.path.join(home_dir,'git','bfd')
integrate_dir = os.path.join(home_dir,'git','bfd')
moment_dir = os.path.join(data3_dir,'sims','moments','g06')
constant_dir = os.path.join(home_dir,'sims','constants','g06')

#absolute paths to files
measure_prog = os.path.join(measure_dir, 'greatMeasure')
integrate_prog = os.path.join(integrate_dir, 'greatIntegrate')
constant_prog = os.path.join(integrate_dir, 'greatConstant')

image_files_temp = os.path.join(image_dir, 'template_stamps_{0}_*.fits'.format(id))
psf_file_temp = os.path.join(image_dir, 'template_stamps_{0}_01.psf'.format(id))
psf_files_temp = os.path.join(image_dir, 'template_stamps_{0}_*.psf'.format(id))
nominal_moment_file = os.path.join(moment_dir, 'stamps_moments_{0}.fits'.format(id))

image_files_targ = os.path.join(image_dir, 'stamps_{0}_*.fits'.format(id))
psf_file_targ = os.path.join(image_dir, 'stamps_{0}_01.psf'.format(id))
psf_files_targ = os.path.join(image_dir, 'stamps_{0}_*.psf'.format(id))

moment_file = os.path.join(moment_dir,
                           'stamps_moments_{0}.fits'.format(id))

f = open(os.path.join(image_dir,'noise_var_{}.txt'.format(id)))
noise_var = f.read()
f.close()

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
                    '-momentFiles={0}'.format(moment_file),
                    '-templateStampSize={0}'.format(stamp_size),
                    '-nominalMomentFile={0}'.format(nominal_moment_file),
                    '-nSample={0}'.format(n_sample),
                    '-selectSnMin={0}'.format(select_sn_min),
                    '-selectSnMax={0}'.format(select_sn_max),
                    '-weightKSigma={0}'.format(weight_sigma),
                    '-priorSigmaStep={0}'.format(prior_sigma_step),
                    '-priorSigmaCutoff={0}'.format(prior_sigma_cutoff)]

#call greatIntegrate
subprocess.call([integrate_prog]+integrate_params)

#call greatConstant
bfile = open(os.path.join(constant_dir,
                          'constant_before_deselect_{0}.txt'.format(id)),'w')

for mb in glob.glob(moment_file):
    out = subprocess.check_output([constant_prog, mb])
    print>>bfile, out

bfile.close()


#parameters for greatIntegrate with Deselection
integrate_params_deselect = ['-templateFiles={0}'.format(image_files_temp),
                             '-templatePsfs={0}'.format(psf_file_temp),
                             '-momentFiles={0}'.format(moment_file),
                             '-templateStampSize={0}'.format(stamp_size),
                             '-nominalMomentFile={0}'.format(nominal_moment_file),
                             '-nSample={0}'.format(n_sample),
                             '-selectSnMin={0}'.format(select_sn_min),
                             '-selectSnMax={0}'.format(select_sn_max),
                             '-weightKSigma={0}'.format(weight_sigma),
                             '-priorSigmaStep={0}'.format(prior_sigma_step),
                             '-priorSigmaCutoff={0}'.format(prior_sigma_cutoff),
                             '-doDeselection=T']

#call greatIntegrate with Deselection
subprocess.call([integrate_prog]+integrate_params_deselect)

#call greatConstant again
afile = open(os.path.join(constant_dir,
                          'constant_after_deselect_{0}.txt'.format(id)),'w')

for ma in glob.glob(moment_file):
    out = subprocess.check_output([constant_prog, ma])
    print>>afile, out

afile.close()

#delete template images
for image in glob.glob(image_files_temp):
    os.remove(image)
for psf in glob.glob(psf_files_temp):
    os.remove(psf)
    
