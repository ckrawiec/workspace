import os
import sys
#from astropy.io import ascii
import numpy as np
import argparse

def main(args):
    parser = argparse.ArgumentParser(description='Run Bruzual & Charlot\'s 2003 gal\
axy evolution code cm_evolution.',
                                     fromfile_prefix_chars='@')
    parser.add_argument('-dust', action='store_true',
                        help='Include for models to use dust attenuation. If includ\
ed, --A option must be specified.')
    parser.add_argument('--IMF', action='store', type=str,
                        help='Initial mass function. Choices are \'salpeter\' or \'\
chabrier\'.')
    parser.add_argument('--H0', action='store', type=float,
                        help='Hubble constant (km/s/Mpc)')
    parser.add_argument('--Omega_m', action='store', type=float,
                        help='Omega matter today.')
    parser.add_argument('--Omega_lambda', action='store', type=float,
                        help='Omega lambda.')
    parser.add_argument('--model', action='store', type=str,
                        help='Theoretical tracks to use. Choices are \'Padova1994\'\
 or \'Padova2000\'.')
    parser.add_argument('--SFhistory', dest='starform', nargs='+',
                        help='Star formation history type. Can add any from this li\
st: SSP, EXP1, EXP2, EXP3, BURST, CONST.')
    parser.add_argument('--metal', dest='metals', nargs='+',
                        help='Metallicity in B&C notation (m42, m52..). Choices wil\
l be appended to list.')
    parser.add_argument('--age', dest='ages', type=float, nargs='+',
                        help='Ages.')
    parser.add_argument('--filtercombo', nargs='+', action='append',
                        help='Filter combinations. Format: (filter1,filter2).')
    parser.add_argument('--constSFR', nargs='+',
                        help='Star formation rates for the CONST option.')
    p = parser.parse_args()


    print "The following options have been chosen for your output:"
    print " "
    print "Starting model:"
    print "    Theoretical tracks from:  ", p.model
    print "    Initial Mass Function:    ", p.IMF
    print "    Metallicities:            ", p.metals
    print "    Star formation histories: ", p.starform
    print "Evolution parameters:"
    print "    Ages:                     ", p.ages
    print "    Filter combinations:      ", p.filtercombo

    if p.dust:
        dustdir = 'dust'
        print "    Dust?                     Yes"
    else:
        dustdir = 'nodust'
        print "    Dust?                     No"
    
    #bc03 stuff
    #filter_coms = [(121,123),(121,122),(122, 123), (123,124), (257,258),(258,259),(259,260),(260,261)]
    #path to place you want all the directories created
    outpath="./output/{}/{}/".format(p.model,dustdir) 
    cm_command= '/Users/Christina/GradCourses/Galaxies/hw2/bc03/src/cm_evolution'
    #path to whichever models you want to use 
    modelname = "./output/{}/{}/%(metal)s/%(sfdir)s/bc2003_hr_%(metal)s_{}_ssp_%(sf)s.ised".format(p.model,dustdir,p.IMF[:4]) 

    for metal in p.metals:
        if not os.path.isdir('%s/%s' %(outpath, metal)):
            os.system('mkdir %s/%s' %(outpath, metal))   
        for sf in p.starform:
            sfdir = sf
            if sf=='CONST':
                #only working for one SFR right now
                for sfr in p.constSFR:
                    sfdir = '%s/%sMsun-yr' %(sf, sfr)
                modelname = os.path.splitext(modelname)[0]+'_'+sfr+'.ised'
            if not os.path.isdir('%s/%s/%s' %(outpath, metal, sfdir)):
                os.system('mkdir %s/%s/%s' %(outpath, metal, sfdir))   
            outdirsf = '{}/{}/{}'.format(outpath,metal,sfdir)
            for age in p.ages:
                outdir = '%s/%.1f' %(outdirsf, age)
                if not os.path.isdir(outdir):
                    os.system('mkdir %s' %(outdir))
        
                for fil_com in p.filtercombo:
                    run_dict = {'metal':metal, 'filter1':int(fil_com[0]), 
                                'filter2':int(fil_com[1]), 'age':age,
                                'sf':sf,
                                'sfdir':sfdir,
                                'outpath':outpath,  
                                'H0':p.H0,'omega_m':p.Omega_m,'omega_l':p.Omega_lambda}
                    run_dict['modelname']=os.path.splitext(modelname)[0] %(run_dict)
                    run_dict['outfiles']=[modelname.replace('ised', a) %run_dict  for a in ['color_F%(filter1)03d_F%(filter2)03d', 'magnitude_F%(filter1)03d','magnitude_F%(filter2)03d']]
                    
                    runfile = 'run_%(metal)s_F%(filter1)03d_F%(filter2)03d_%(age).1f.dat' %run_dict
                    outfile = open(runfile, 'w')
                    outfile.write("""%(H0).1f, %(omega_m)0.2f, %(omega_l)0.2f
%(age)2.1f
%(modelname)s
%(filter1)03d, %(filter2)03d

""" %run_dict)
                    outfile.close()
                
                    os.system('%s < %s' %(cm_command, runfile))
                    
                    for output in run_dict['outfiles']+[runfile,]:
                        cmd = ('mv %s %s\n' %(output, outdir))
                        os.system(cmd)

if __name__=="__main__":
    main(sys.argv)
