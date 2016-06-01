import ezgal
import numpy as np
import glob
import os

#Make different FITS table (HDUs?) for each model set?
#Make different FITS files for each IMF?

#Have column to note that models are interpolated?
#Or just have no interpolations

table_file = 'bc03_colors_01.tab'

#cosmology
h = 0.7
H0 = h*100. #km/s/Mpc
Om = 0.3
Ol = 0.7

#ranges of values used to create table
#currently uses all metallicities available with no interpolations
z_form_range = np.arange(0.1,1.,0.5)

#exponential decline (tau) and burst (burst_length) SFHs are included with EzGal
#burst has issues at this time, leave range empty
tau_range = [0.1]
burst_length_range = []

#additional SFHs (define below and specify in main() )
sfr_range = [5.0]

#dust function (define below and specify in main(), one per table)
dust_params = (1.0, 0.3)

#lines are 120 chars wide
#lines to go in table that explain the dust model and parameters being used with it
dust_lines = ['#                      Dust model                : Charlot & Fall 2000                                                  #\n',
              '#                      Parameters for dust model : (tauV, mu)                                                           #\n',
              '#                                                  tauV = V-band optical depth for ages <= 10^7 yr                      #\n',
              '#                                                  mu   = fraction of tav for ages > 10^7 yr                            #\n']


#DES filters
u = '/Users/Christina/DES/DES_tput_Ang_u.txt'
g = '/Users/Christina/DES/DES_tput_Ang_g.txt'
r = '/Users/Christina/DES/DES_tput_Ang_r.txt'
i = '/Users/Christina/DES/DES_tput_Ang_i.txt'
z = '/Users/Christina/DES/DES_tput_Ang_z.txt'
Y = '/Users/Christina/DES/DES_tput_Ang_Y.txt'

bc03_models = glob.glob('/Users/Christina/anaconda/lib/python2.7/site-packages/ezgal/data/models/bc03*')
bc03_ssp_models = glob.glob('/Users/Christina/anaconda/lib/python2.7/site-packages/ezgal/data/models/bc03_ssp_z_*_chab.model')

#EzGal utility
def setmodel(mo):
    mo.set_cosmology(Om=Om,Ol=Ol,h=h)
    mo.add_filter(u,name='u')
    mo.add_filter(g,name='g')
    mo.add_filter(r,name='r')
    mo.add_filter(i,name='i')
    mo.add_filter(z,name='z')
    mo.add_filter(Y,name='Y')

#star formation laws
#exponential and burst included in ezgal
def const(t, sfr):
    return t*0 + sfr

#dust attenuation laws
def CharlotFall(age, wavelength, tauV=1.0, mu=0.3):
    """
    Returns extinction values from Charlot & Fall 2000
        in array of size (???)
    Arguments:
        age (array, in Gyr)
        wavelength (array, in Angstroms)
        tauV (float, V-band optical depth when age <= 10^7 yr)
        mu (float, fraction of tauV when age > 10^7 yr)
    """
    tau = np.array([0]*len(age))
    tau[age<=0.01]=tauV
    tau[age>0.01]=mu*tauV
    return np.array([np.exp(-tau * (wi/5500.)**-0.7) for wi in wavelength])

#functions for main
def maketableentries(model, model_id, model_name, IMF, metallicity, star_formation, dust, dust_params, composite):
    print "Working on model #", model_id
    for zf in z_form_range:
        zs = model.get_zs(zf)
        l = len(zs)

        table_dict['input_model'].extend([os.path.basename(model_name)]*l)
        table_dict['model_id'].extend([model_id]*l)
        table_dict['IMF'].extend([IMF]*l)
        table_dict['metallicity'].extend([metallicity]*l)
        table_dict['star_formation'].extend([star_formation]*l)
        table_dict['dust'].extend([dust]*l)
        table_dict['dust_params'].extend([dust_params]*l)
        table_dict['z'].extend(zs)
        table_dict['z_form'].extend([zf]*l)
        table_dict['age'].extend(model.get_age(zf, zs))
        table_dict['DES_u'].extend(model.get_apparent_mags(zf, zs=zs, filters='u'))
        table_dict['DES_g'].extend(model.get_apparent_mags(zf, zs=zs, filters='g'))
        table_dict['DES_r'].extend(model.get_apparent_mags(zf, zs=zs, filters='r'))
        table_dict['DES_i'].extend(model.get_apparent_mags(zf, zs=zs, filters='i'))
        table_dict['DES_z'].extend(model.get_apparent_mags(zf, zs=zs, filters='z'))
        table_dict['DES_Y'].extend(model.get_apparent_mags(zf, zs=zs, filters='Y'))
        table_dict['composite'].extend([composite]*l)


if __name__=="__main__":

    sf_funcs = [const]
    dust_func = CharlotFall

    #dict which will be written to file later
    table_dict = {'input_model'    : [],
                  'model_id'       : [],
                  'IMF'            : [],
                  'metallicity'    : [],
                  'star_formation' : [],
                  'dust'           : [],
                  'dust_params'    : [],
                  'z'              : [],  
                  'z_form'         : [],
                  'age'            : [],
                  'DES_u'          : [],
                  'DES_g'          : [],
                  'DES_r'          : [],
                  'DES_i'          : [],
                  'DES_z'          : [],
                  'DES_Y'          : [],
                  'composite'      : []}

    model_id = 1
    
    #get colors from included models first
    for model_name in bc03_models:
        model = ezgal.ezgal(model_name)
        setmodel(model)

        name_split = model_name.split('_')
        star_formation = "_".join(name_split[1:-3])
        metallicity = name_split[-2]
        IMF = name_split[-1].split('.')[0]

        #dust stuff untrue for some, should change
        maketableentries(model, model_id, model_name, IMF, metallicity, star_formation,
                         dust=False, dust_params=None,
                         composite=False)
        model_id +=1
        
    #cycle through star formations and their parameters
    #save models?
    #then run through those models with dust and different dust params

    #different model for each metallicity
    for ssp_model_name in bc03_ssp_models:
        model = ezgal.ezgal(ssp_model_name)
        setmodel(model)
        
        name_split = ssp_model_name.split('_')
        metallicity = name_split[-2]
        IMF = name_split[-1].split('.')[0]

        for sf_func in sf_funcs:
            if sf_func==const:
                for sfr in sfr_range:
                    star_formation = 'const_{}'.format(sfr)
                    csp = model.make_csp(const, args=(sfr,))
                    setmodel(csp)
                    maketableentries(csp, model_id, ssp_model_name, IMF, metallicity, star_formation,
                                     dust=False, dust_params=None,
                                     composite=True)
                    model_id +=1
                    
                    dusty_csp = model.make_csp(const, args=(sfr,), dust_function=dust_func, dust_args=dust_params)
                    setmodel(dusty_csp)
                    maketableentries(csp, model_id, ssp_model_name, IMF, metallicity, star_formation,
                                     dust=True, dust_params=dust_params,
                                     composite=True)
                    model_id +=1
                    
        for burst_length in burst_length_range:
            star_formation = 'burst_{}'.format(burst_length)
            csp = model.make_burst(burst_length)
            setmodel(csp)
            maketableentries(csp, model_id, ssp_model_name, IMF, metallicity, star_formation,
                             dust=False, dust_params=None,
                             composite=True)
            model_id +=1
                
            dusty_csp = model.make_burst(burst_length, dust_function=dust_func, dust_args=dust_params)
            setmodel(dusty_csp)
            maketableentries(csp, model_id, ssp_model_name, IMF, metallicity, star_formation,
                             dust=True, dust_params=dust_params,
                             composite=True)
            model_id +=1

        for tau in tau_range:
            star_formation = 'exp_{}'.format(tau)
            csp = model.make_exponential(tau)
            setmodel(csp)
            maketableentries(csp, model_id, ssp_model_name, IMF, metallicity, star_formation,
                             dust=False, dust_params=None,
                             composite=True)
            model_id +=1
            
            dusty_csp = model.make_exponential(tau, dust_function=dust_func, dust_args=dust_params)
            setmodel(dusty_csp)
            maketableentries(csp, model_id, ssp_model_name, IMF, metallicity, star_formation,
                             dust=True, dust_params=dust_params,
                             composite=True)
            model_id +=1
                    
    #write table to file
    #choose a dust model?

    
    f = open(table_file,'w')
    #lines are 120 chars wide
    f.write('########################################################################################################################\n')
    f.write('#                                                                                                                      #\n')
    f.write('#                      DES colors for spectra created with EzGal (http://www.baryons.org/ezgal/).                      #\n')
    f.write('#                                                                                                                      #\n')
    f.write('#                                 Cosmology: Omega_matter={}, Omega_lambda={}, h={}                                 #\n'.format(Om,Ol,h))
    f.write('#                                                                                                                      #\n')
    for dust_line in dust_lines:
        f.write(dust_line)
    f.write('#                                                                                                                      #\n')
    f.write('#     input_model    : Name of SPS model set.                                                                          #\n')
    f.write('#                      Example: \'BC03\' means Bruzual & Charlot 2003.                                                 #\n')
    f.write('#     model_id       : ID for this table. Each ID corresponds to a single SSP or CSP.                                  #\n')
    f.write('#     IMF            : Initial mass function.                                                                          #\n')
    f.write('#                      Example: \'chab\' means Chabrier.                                                               #\n')
    f.write('#     metallicity    : Metallicity. Z_sun=0.02.                                                                        #\n')
    f.write('#     dust           : True if dust included in modeled spectrum.                                                      #\n')
    f.write('#     dust_params    : Parameters unique to the dust attenuation law chosen for this table.                            #\n')
    f.write('#     star_formation : Star formation history.                                                                         #\n') 
    f.write('#                      Example: \'exp_1.0\' refers to an exponentially decaying star formation law with tau=1.0 Gyr.   #\n')
    f.write('#     z              : Redshift.                                                                                       #\n')
    f.write('#     z_form         : Redshift at time of formation.                                                                  #\n')
    f.write('#     age            : Age in Gyr.                                                                                     #\n')
    f.write('#     DES_u          : AB magnitude in Dark Energy Camera\'s u filter.                                                 #\n')
    f.write('#     DES_g          : AB magnitude in Dark Energy Camera\'s g filter.                                                 #\n')
    f.write('#     DES_r          : AB magnitude in Dark Energy Camera\'s r filter.                                                 #\n')
    f.write('#     DES_i          : AB magnitude in Dark Energy Camera\'s i filter.                                                 #\n')
    f.write('#     DES_z          : AB magnitude in Dark Energy Camera\'s z filter.                                                 #\n')
    f.write('#     DES_Y          : AB magnitude in Dark Energy Camera\'s Y filter.                                                 #\n')
    f.write('#     composite      : False if model is included with EzGal.                                                          #\n')
    f.write('#                                                                                                                      #\n')
    f.write('########################################################################################################################\n')
    f.write('input_model IMF metallicity dust dust_params star_formation z z_form age DES_u DES_g DES_r DES_i DES_z DES_Y composite\n')
    
    for row_num in range(len(table_dict['input_model'])):
        f.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(table_dict['input_model'][row_num],
                                                                           table_dict['model_id'][row_num],
                                                                           table_dict['IMF'][row_num],
                                                                           table_dict['metallicity'][row_num],
                                                                           table_dict['dust'][row_num],
                                                                           table_dict['dust_params'][row_num],
                                                                           table_dict['star_formation'][row_num],
                                                                           table_dict['z'][row_num],
                                                                           table_dict['z_form'][row_num],
                                                                           table_dict['age'][row_num],
                                                                           table_dict['DES_u'][row_num],
                                                                           table_dict['DES_g'][row_num],
                                                                           table_dict['DES_r'][row_num],
                                                                           table_dict['DES_i'][row_num],
                                                                           table_dict['DES_z'][row_num],
                                                                           table_dict['DES_Y'][row_num],
                                                                           table_dict['composite'][row_num]))
    f.close()
