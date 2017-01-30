#
#Gather and join SVA1 Balrog tables from dessci db
#
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import easyaccess as ea
from astropy.table import Table, join

#directory to store outputs
data_dir = '/Users/Christina/DES/data/balrog/sva1/'

#format of balrog tables in db
table_format = "JELENA.BALROG_SVA1_TAB{}_{}_{}"

table_nums = [str(r).zfill(2) for r in range(1,12)]
#09 is not valid
table_nums.remove('09')

bands = ['G','R','I','Z']
tables = ['SIM', 'TRUTH']#, 'NOSIM']

connection = ea.connect()

for table in tables:

    for table_num in table_nums:
        if table=='SIM':
            out_file = data_dir+'balrog_sva1_detmodel_tab{}_SIM.fits'.format(table_num)
        elif table=='TRUTH':
            out_file = data_dir+'balrog_sva1_tab{}_TRUTH.fits'.format(table_num)

        if os.path.exists(out_file):
            print "Output path {} already exists. Delete before running again. Moving to next table.".format(out_file)
            continue

        band_dfs = []

        for band in bands:            
            table_name = table_format.format(table_num, table, band)
            if table=='SIM':
                query = "select BALROG_INDEX, ALPHAMODEL_J2000, DELTAMODEL_J2000, FLUX_DETMODEL, FLUXERR_DETMODEL, MAG_DETMODEL, MAGERR_DETMODEL, FLAGS from {};".format(table_name)
                print "\nsubmitting query: "
                print "    ", query
                band_df = connection.query_to_pandas(query)
                band_df = band_df.rename(index=str, columns={'ALPHAMODEL_J2000':'ALPHAMODEL_J2000_'+band,
                                                             'DELTAMODEL_J2000':'DELTAMODEL_J2000_'+band,
                                                             'FLUX_DETMODEL':'FLUX_DETMODEL_'+band,
                                                             'FLUXERR_DETMODEL':'FLUXERR_DETMODEL_'+band,
                                                             'MAG_DETMODEL':'MAG_DETMODEL_'+band,
                                                             'MAGERR_DETMODEL':'MAGERR_DETMODEL_'+band,
                                                             'FLAGS':'FLAGS_'+band})
                band_dfs.append(band_df)

            elif table=='TRUTH':
                query = "select BALROG_INDEX, ID, OBJTYPE, RA, DEC, Z, ZEROPOINT, MAG, FLUX_0, HALFLIGHTRADIUS_0, FLUX_NOISELESS, FLUX_NOISED from {};".format(table_name)
                band_df = connection.query_to_pandas(query)
                band_df = band_df.rename(index=str, columns={'ZEROPOINT':'ZEROPOINT_'+band,
                                                             'MAG':'MAG_'+band,
                                                             'FLUX_0':'FLUX_0_'+band,
                                                             'HALFLIGHTRADIUS_0':'HALFLIGHTRADIUS_0_'+band,
                                                             'FLUX_NOISELESS':'FLUX_NOISELESS_'+band,
                                                             'FLUX_NOISED':'FLUX_NOISED_'+band})
                band_dfs.append(band_df)
                print "\nsubmitting query: "
                print "    ", query

        merged_df = band_dfs[0]
        for band_df in band_dfs[1:]:
            if table=='SIM':
                merged_df = pd.merge(merged_df, band_df, how='outer', on='BALROG_INDEX')
            elif table=='TRUTH':
                merged_df = pd.merge(merged_df, band_df, how='outer', on=['BALROG_INDEX',
                                                                          'ID',
                                                                          'RA',
                                                                          'DEC',
                                                                          'Z',
                                                                          'OBJTYPE'])
                                                                         
        tab = Table()
        for col in merged_df.columns:
            tab.add_column(Table.Column(merged_df[col], col))

        tab.write(out_file)
    
connection.close()
