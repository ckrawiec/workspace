import meds
import sys
import matplotlib.pyplot as plt
import numpy as np

if __name__=="__main__":
    filename = sys.argv[1]

    # create a MEDS object for the given MEDS file
    m=meds.MEDS(filename)

    # open log file
    log = open('testDESMEDS.log','w')
    
    log.write('MEDS file name: {}\n'.format(filename))

    # print the number of objects
    nobj = m.size
    log.write('number of objects: {}\n'.format(nobj))

    # read a cutout for object 35, cutout index 2 (0-indexed)
    index = 35
    cutout_index = 2
    image = m.get_cutout(index, cutout_index)
    
    # get other image types
    seg  = m.get_cutout(index, cutout_index, type='seg')
    wt   = m.get_cutout(index, cutout_index, type='weight')
    mask = m.get_cutout(index, cutout_index, type='bmask')

    # get a python list of all cutouts for this object
    imlist  = m.get_cutout_list(index)
    seglist = m.get_cutout_list(index,type='seg')
    
    # save images to files
    cutout_file = '/Users/Christina/DES/testDESMEDS_cutoutpython'
    seg_file = '/Users/Christina/DES/testDESMEDS_segpython'
    weight_file = '/Users/Christina/DES/testDESMEDS_weightpython'
    mask_file = '/Users/Christina/DES/testDESMEDS_maskpython'
    np.save(cutout_file, image)
    np.save(seg_file, seg)
    np.save(weight_file, wt)
    np.save(mask_file, mask)

    # report file names of saved cutouts
    log.write('\n')
    log.write('opening cutout index {} from object index {}\n'.format(cutout_index, index))
    log.write('cutout array written to: {}\n'.format(cutout_file))
    log.write('seg array written to: {}\n'.format(seg_file))
    log.write('weight array written to: {}\n'.format(weight_file))
    log.write('mask array written to: {}\n'.format(mask_file))

    # get the jacobian of the WCS transformation
    # as a dict
    j = m.get_jacobian(index, cutout_index)
        
    # as a numpy matrix
    jmatrix = m.get_jacobian_matrix(index, cutout_index)
        
    # get the "ubserseg" weight map
    wt = m.get_cweight_cutout_nearest(index, cutout_index)

    # report affine and noise for test cutout
    log.write('\njacobian of the WCS transformation: {}\n\n'.format(j))

    # print number of cutouts per object
    log.write('object #, number of cutouts\n')
    log.write('---------------------------\n')
    for iobj in range(nobj):
        ncutout = m['ncutout'][iobj]
        log.write('{}, {}\n'.format(iobj, ncutout))


    log.close()
