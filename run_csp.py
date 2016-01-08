import os
import sys
import glob
import sys
import argparse
from subprocess import call

def main(args):
    parser = argparse.ArgumentParser(description='Run Bruzual & Charlot\'s 2003 galaxy evolution code csp_galaxev.',
                                     fromfile_prefix_chars='@')
    parser.add_argument('-dust', action='store_true',
                        help='Include for models to use dust attenuation. If included, --A option must be specified.')
    parser.add_argument('--IMF', action='store', type=str,
                        help='Initial mass function. Choices are \'salpeter\' or \'chabrier\'.')
    parser.add_argument('--model', action='store', type=str,
                        help='Theoretical tracks to use. Choices are \'Padova1994\' or \'Padova2000\'.')
    parser.add_argument('--SFhistory', dest='starform', nargs='+',
                        help='Star formation history type. Can add any from this list: SSP, EXP1, EXP2, EXP3, BURST, CONST.')
    parser.add_argument('--constSFR', nargs='+',
                        help='Star formation rate for use with the \'CONST\' SFhistory option.')
    parser.add_argument('--metal', dest='metals', nargs='+',
                        help='Metallicity in B&C notation (m42, m52..). Choices will be appended to list.')
    p = parser.parse_args()
        
    #bc03 stuff
    #path to whichever models you want to use
    modelname = "./models/{}/{}/bc2003_hr_%(metal)s_{}_ssp.ised".format(p.model,p.IMF,p.IMF[:4])
    outpath="./output/{}".format(p.model) #path to place you want all the directories created
    csp_command= 'src/csp_galaxev' 

    print "The following options have been chosen for your output models:"
    print "    Theoretical tracks from:  ", p.model    
    print "    Initial Mass Function:    ", p.IMF
    print "    Metallicities:            ", p.metals
    print "    Star formation histories: ", p.starform
    if p.dust:
        dust = 'Y'
        dustdir = 'dust'
        print "    Dust?                 Yes"        
    else:
        dust = 'N'
        dustdir = 'nodust'
        print "    Dust?                 No"        
    print " "
    print "Output files will be saved in ", outpath

    #make sure output directories exist
    if not os.path.isdir('%s/%s' %(outpath, dustdir)):
        os.system('mkdir %s/%s' %(outpath, dustdir))
    outdird = os.path.join(outpath, dustdir)       
    
    for metal in p.metals:
        if not os.path.isdir('%s/%s' %(outdird, metal)):
            os.system('mkdir %s/%s' %(outdird, metal))   
        outdirm = os.path.join(outdird,metal)
        for sf in p.starform:
            if not os.path.isdir('%s/%s' %(outdirm, sf)):
                os.system('mkdir %s/%s' %(outdirm, sf))   
            outdir = os.path.join(outdirm,sf)
            if sf == 'CONST':
                for sfr in p.constSFR:
                    if not os.path.isdir('%s/CONST/%sMsun-yr' %(outdirm,sfr)):
                        os.system('mkdir %s/CONST/%sMsun-yr' %(outdirm,sfr))
                    outdir = os.path.join(outdir,sfr+'Msun-yr')
                    
                    run_dict = {'metal':metal,
                                'dust':dust,
                                'sf':'CONST',
                                'sfr':sfr
                                }
                    model = os.path.splitext(modelname)[0]
                        
                    #make sure .ised file exists, if not, use bin_ised
                    if not os.path.exists(modelname %run_dict):
                        if os.path.exists(model+'.ised_ASCII'):
                            call(["src/bin_ised",modelname %run_dict])
                        else:
                            print "Model name '{}' does not exist in either .ised or .ised_ASCII form.".format(modelname%run_dict)

                    run_dict['modelname']=modelname %run_dict
                    run_dict['model']=model %run_dict
                    run_dict['outputname']='%(model)s_CONST_%(sfr)s' %run_dict
                    run_dict['outfiles']=glob.glob('%(outputname)s*'%run_dict)
                    
                    runfile = ('run_%(metal)s_'+dustdir+'_CONST_%(sfr)s.dat') %run_dict
                    outfile = open(runfile, 'w')
                    #TCUT(SFR=0) at 20Gyr                                                                                                     
                    outfile.write("""%(modelname)s
%(dust)s                                                                                                                              
                                                                                                                                      
3                                                                                                                                     
%(sfr)s                                                                                                                               
                                                                                                                                      
%(outputname)s                                                                                                                        
""" %run_dict)
                    print '%s < %s' %(csp_command, runfile)
                    os.system('src/csp_galaxev < {}'.format(runfile))
                    exit()
                    for output in run_dict['outfiles']+[runfile,]:
                        cmd = ('mv %s %s\n' %(output, outdir))
                        os.system(cmd)
            else:
                run_dict = {'metal':metal,
                            'dust':dust,
                            'sf':sf,
                            }
                model = os.path.splitext(modelname)[0]

                #make sure .ised file exists, if not, use bin_ised                                                                          
                if not os.path.exists(modelname %run_dict):
                    if os.path.exists(model+'.ised_ASCII'):
                        call(["src/bin_ised",modelname %run_dict])
                    else:
                        print "Model name '{}' does not exist in either .ised or .ised_ASCII form.".format(modelname%run_dict)

                run_dict['modelname']=modelname %run_dict
                run_dict['model']=model %run_dict
                run_dict['outputname']='%(model)s_output_%(sf)s' %run_dict
                run_dict['outfiles']=glob.glob('%(outputname)s*'%run_dict)
                
                runfile = ('run_%(metal)s_'+dustdir+'_%(sf)s.dat') %run_dict
                outfile = open(runfile, 'w')
                outfile.write("""%(modelname)s
%(dust)s

""" %run_dict)
                if sf=='SSP':
                    outfile.write("""0
%(outputname)s

""" %run_dict)
                if sf[:3]=='EXP':
                    run_dict['tau']=int(sf[-1])
                    run_dict['starnum']=1
                    
                    outfile.write("""1
%(tau)d


%(outputname)s

""" %run_dict)
                if sf=='BURST':
                    #single burst of length 2 Gyr
                    outfile.write("""2
2
%(outputname)s

""" %run_dict)
                    
                outfile.close()
            
                os.system('%s < %s' %(csp_command, runfile))
                
                for output in run_dict['outfiles']+[runfile,]:
                    cmd = ('mv %s %s\n' %(output, outdir))
                    os.system(cmd)

if __name__=="__main__":
    main(sys.argv)
