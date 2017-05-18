#
#Gather and join SVA1 Balrog tables from dessci db
#
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import easyaccess as ea
from astropy.table import Table, join, Column

#directory to store outputs
data_dir = '/Users/Christina/DES/data/'

#format of balrog tables in db
table_format = "JELENA.BALROG_SVA1_TAB{}_{}_{}"

table_nums = [str(r).zfill(2) for r in range(1,12)]
#09 is not valid
table_nums.remove('09')

bands = ['G','R','I','Z']
tables = ['TRUTH']#, 'NOSIM']

connection = ea.connect()

for table in tables:

    for table_num in table_nums:
        if table=='SIM':
            out_file = data_dir+'balrog_sva1_aper_10_tab{}_SIM.fits'.format(table_num)
        elif table=='TRUTH':
            out_file = data_dir+'balrog_sva1_tab{}_TRUTH.fits'.format(table_num)

        if os.path.exists(out_file):
            print "Output path {} already exists. Delete before running again. Moving to next table.".format(out_file)
            continue

        band_dfs = []

        for band in bands:            
            table_name = table_format.format(table_num, table, band)
            if table=='SIM':
                query = "select BALROG_INDEX, ALPHAMODEL_J2000, DELTAMODEL_J2000, FLUX_APER_10, FLUXERR_APER_10, MAG_APER_10, MAGERR_APER_10, FLAGS from {};".format(table_name)
                print "\nsubmitting query: "
                print "    ", query
                band_df = connection.query_to_pandas(query)
                band_df = band_df.rename(index=str, columns={'ALPHAMODEL_J2000':'ALPHAMODEL_J2000_'+band,
                                                             'DELTAMODEL_J2000':'DELTAMODEL_J2000_'+band,
                                                             'FLUX_APER_10':'FLUX_APER_10_'+band,
                                                             'FLUXERR_APER_10':'FLUXERR_APER_10_'+band,
                                                             'MAG_APER_10':'MAG_APER_10_'+band,
                                                             'MAGERR_APER_10':'MAGERR_APER_10_'+band,
                                                             'FLAGS':'FLAGS_'+band})
                band_dfs.append(band_df)

            elif table=='TRUTH':
                query = "select BALROG_INDEX, ID, OBJTYPE, RA, DEC, Z, ZEROPOINT, MAG, FLUX_0, FLUX_NOISELESS, FLUX_NOISED, NOT_DRAWN, INDEXSTART, MOD, X, Y, TILENAME, SERSICINDEX_0, HALFLIGHTRADIUS_0, AXISRATIO_0, BETA_0 from {};".format(table_name)
                band_df = connection.query_to_pandas(query)
                band_df = band_df.rename(index=str, columns={'ZEROPOINT':'ZEROPOINT_'+band,
                                                             'MAG':'MAG_'+band,
                                                             'FLUX_0':'FLUX_0_'+band,
                                                             'NOT_DRAWN':'NOT_DRAWN_'+band,
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
                                                                          'INDEXSTART',
                                                                          'ID',
                                                                          'OBJTYPE',
                                                                          'MOD',
                                                                          'RA',
                                                                          'DEC',
                                                                          'X','Y',
                                                                          'Z',
                                                                          'TILENAME',
                                                                          'SERSICINDEX_0',
                                                                          'HALFLIGHTRADIUS_0',
                                                                          'AXISRATIO_0',
                                                                          'BETA_0'])


        tab = Table()
        for col in merged_df.columns:
            tab.add_column(Column(merged_df[col], col))
        
        if 'TILENAME' in tab.colnames:
            new_col = Column([str(tilename) for tilename in tab['TILENAME']], 'TILENAME')
            tab.remove_column('TILENAME')
            tab.add_column(new_col)

        tab.write(out_file)
    
connection.close()
