import zprobability as zprob
import ConfigParser
import time
import os

def parseconfig(config_file):
    config = ConfigParser.SafeConfigParser(config_file)

    params = {}

    params['template_file'] = config.get('template_file')
    params['target_file'] = config.get('target_file')
    params['output_file'] = config.get('output_file')
    if os.path.exists(params['output_file']):
        sys.exit("Output file already exists. Delete before running again.\n    Exiting.")
    
    params['filters'] = config.get('filters')
    params['redshift_groups'] = config.get('redshift_groups')

    params['integration_type'] = config.get('integration_type')
    if params['integration_type'] == 'tree':
        params['k_near'] = config.getint('k_near')
        
    params['template_id_column'] = config.get('template_id_column')
    params['template_data_column'] = config.get('template_data_column')
    params['redshift_column'] = config.get('redshift_column')
    
    params['target_id_column'] = config.get('target_id_column')
    params['target_data_column'] = config.get('target_data_column')
    params['target_error_column'] = config.get('target_error_column')

    #if start/end indices not specified, use all targets
    if config.has_option('target_start_index'):
        params['target_start_index'] = config.getint('target_start_index')
    else:
        params['target_start_index'] = 0

    if config.has_option('target_end_index'):
        params['target_end_index'] = config.getint('target_end_index')
    else:
        params['target_end_index'] = None
            
    params['debug'] = config.getboolean('debug')
    if params['debug'] == True:
        params['N_debug'] = config.getint('N_debug')

def printtime():
    now = time.strftime("%Y-%m-%d %H:%M")
    print "#"+now

def main(args):
    parseconfig(args[1])

    printtime()

    print "num_threads="+str(num_threads)
    print "Calculating {} probabilities, using {} integrals".format(params['target_data_type'],
                                                                    params['ptype'])
    
    os.environ['OMP_NUM_THREADS']=str(num_threads)

    setup_start = time.time()

    templates = zprob.Templates(params['template_file'],
                                params['template_id_column'],
                                params['template_data_column'])

    targets = zprob.Targets(params['target_file'],
                            params['target_id_column'],
                            params['target_data_column'],
                            params['target_error_column'])
    

    setup_time = time.time()-setup_start
    print "Total setup time took {} s".format(setup_time)
    
    P_dict = {}
        
    if debug==True:
        targets.debug(templates,
                      params['N_debug'],
                      params['output_file'])
        



    N_try = len(data_zip)

    n_per_process = int( np.ceil(N_try/ float(num_threads)) )
    print "n_per_process = ", n_per_process
    data_chunks = [data_zip[i:i+n_per_process] for i in xrange(0, N_try, n_per_process)]
    print "Length of data_zip, N_try: {}, {}".format(len(data_zip), N_try)
    del data_zip

    err_chunks = [err_zip[i:i+n_per_process] for i in xrange(0, N_try, n_per_process)]
    del err_zip

    print "Length of data_chunks: {}".format(len(data_chunks))
    print "Working on {} galaxies (table indices {}-{}) ...".format(N_try, start_index, end_index)
    
    #multiprocessing
    pool = Pool(processes=num_threads)

    for z_group in z_groups:
        z_mask = (z >= np.min(z_group)) & (z < np.max(z_group))

        template_zip = np.array( zip( *[templates[data_type+'_'+SE_col+'_'+f][z_mask] for f in filters] ) )
        print "# of COSMOS templates in z group {}: {}".format(str(z_group), len(template_zip))
    
        start = time.time()
        
        results = pool.map(pwrapper, itertools.izip(data_chunks, err_chunks, itertools.repeat(template_zip)))
        
        P_dict[str(z_group)] = np.concatenate(results)

        work_time = time.time() - start
        print "Work completed in {} s".format(work_time)

    pool.close()
        
    #write results to fits file
    col_defs = [fits.Column(name=id_col_name, format='K', array=data_ids)]

    P_norm = np.zeros(N_try)
    for k in P_dict.keys():
        P_norm += P_dict[k]
    for z_group in z_groups:
        col_defs.append(fits.Column(name='P'+str(z_group), format='D', array=P_dict[str(z_group)]/P_norm))

    col_defs.append(fits.Column(name='PNORM', format='D', array=P_norm))

    pri_hdr = fits.Header()
    tb_hdr = fits.Header()
    tb_hdr['COMMENT'] = "Bayesian redshift probabilities for data in {} using photo-zs of templates from {}. Data vectors/errors were comprised of {}(ERR)_{}_% for % in {}. Columns reported here are \'P[zmin, zmax]\'".format(data_file, template_file, data_type, SE_col, filters)
    if ptype=='tree':
        tb_hdr['NTREE'] = str(k_near)

    pri_hdu = fits.PrimaryHDU(header=pri_hdr)
    tb_hdu = fits.BinTableHDU.from_columns(fits.ColDefs(col_defs), nrows=N_try, header=tb_hdr)
    hdu_list = fits.HDUList([pri_hdu, tb_hdu])
    hdu_list.writeto(output_file, clobber=True)
    
    now = time.strftime("%Y-%m-%d %H:%M")
    print "#"+now
        
if __name__=="__main__":
    main(sys.argv)
