from astropy.table import Table, join

mu = 1.04

mag_file = '/home/ckrawiec/DES/magnification/balrog_y1a1_truth_noiseless_flux_size_objtype1_80M_mag_04.fits'
truth_file = '/home/ckrawiec/DES/data/balrog_y1a1_truth_noiseless_flux_size_objtype1.fits'
sim_file = '/home/ckrawiec/DES/data/balrog_y1a1_truth_sim_flux_detmodel.fits'
zprob_file = '/home/ckrawiec/DES/magnification/lbgselect/zprob_balrog_y1a1_z25_3bins_sigma_tree_griz.fits'

def truth_plus_mag_by_orig():
    """
    Return the truth table inner joined to the
    magnified table by original balrog index
    """
    bmag = Table.read(mag_file)
    bal = Table.read(truth_file)
    print "joining truth table to magnified table by original index..."
    return join(bal,bmag)

def truth_plus_mag_by_mag():
    """
    Return the truth table inner joined to the
    magnified table by balrog index of
    magnified match
    """
    bal = Table.read(truth_file)
    bal.rename_column('BALROG_INDEX', 'BALROG_INDEX_MAG'+str(mu))
    bmag = Table.read(mag_file)
    print "joining truth table to magnified table by magnified index..."
    return join(bal,bmag)

def sim_plus_mag_by_orig():
    """
    Return the sim table outer joined
    to the magnified table by balrog index
    of the original object
    """
    sim = Table.read(sim_file)
    bmag = Table.read(mag_file)
    print "joining sim table to magnified table by original index..."
    return join(bmag, sim, join_type='outer')

def sim_plus_mag_by_mag():
    """
    Return the sim table outer joined
    to the magnified table by balrog index
    of the magnified match
    """
    sim = Table.read(sim_file)
    sim.rename_column('BALROG_INDEX', 'BALROG_INDEX_MAG'+str(mu))
    bmag = Table.read(mag_file)
    print "joining sim table to magnified table by magnified index..."
    return join(bmag, sim, join_type='outer')

def unlensed_detection():
    det = sim_plus_mag_by_orig()
    m = det['FLUX_DETMODEL_G'].mask
    n0 = len(det[~m])
    print "{}/{} detected".format(n0, len(det))
    return n0

def zprob_plus_sim_plus_mag_by_orig():
    det = sim_plus_mag_by_orig()
    zbal = Table.read(zprob_file)
    print "joining sim, mag, and zprob tables by original index..."
    return join(det, zbal)

def magnified_detection():
    zdet = zprob_plus_sim_plus_mag_by_orig()
    src = zdet['P[2.5,9.9]']>0.6
    m = zdet['FLUX_DETMODEL_G'].mask & src
    n0 = len(zdet[~m])

    print "For galaxies with P[2.5,9.9]>0.6: "
    print "   Number detected with original flux: ", n0

def zprob_plus_sim_plus_mag_by_mag():
    det = sim_plus_mag_by_mag()
    zbal = Table.read(zprob_file)
    print "joining sim, mag, and zprob tables by magnified index..."
    return join(det, zbal)

#joining on original balrog_index to get unmagnified zsrc
magzdet = join(magdet, zbal)
del magdet, zbal

print "   #Number detected with magnified flux: ", nmag
print "   #alpha = "#, (nmag-n0)/(mu-1)


#did zprob use all or just detected

del det, zdet

nmag = len(magdet[~magdet['FLUX_DETMODEL_G'].mask])
magsrc = magzdet['P[2.5,9.9]']>0.6
nmag = len(magzdet[~magzdet['FLUX_DETMODEL_G'].mask & magsrc])

print "#   Number detected with magnified flux: ", nmag
print "   alpha = "#, (nmag-n0)/(mu-1)

#for a given truth object... was it detected
#was its magnified version detected

