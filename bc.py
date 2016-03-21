import os
import numpy as np
import scipy.integrate as integ
from scipy import interpolate
import matplotlib.pyplot as plt
from subprocess import call
from astropy.io import ascii, fits

#constants
#flat universe
h=0.72
H0=h*100*3.24e-20 #1/s
Om = 0.27
Ol = 0.73
Or = 4.2e-5/h**2.
cH = (3.e5)/(h*100) #Mpc

cA = 3.e18 #Angstrom/s

ABflux = 3.631e-20 #ergs/s/cm^2/Hz
Lsun = 3.826e33 #erg/s
 
colorfit = True

IMF = "chabrier"
model = "Padova2000"
dust = 'nodust' #'dust' or 'nodust'
tauV="1."
mu="0.3"
metals = ["m142","m152","m162","m172"]
starform = ["BURST","CONST"]
mass = np.logspace(11.9,13.,5)
constSFR = ["10","20","30","40","50"]
burstlength = ["1","2","3"]
filtcoms = [(258,259),(259,260),(260,261)]

pltd = {'x' : 'redshift',
        'y' : 'g-r',
        'color' : 'age (Gyr)',
        'cutmag' : 'z',
        'maglimitlow' : 0.,
        'plotdir' : 'plots/highz/{}/'.format(dust),
        'redmagic' : False,
        'xlim' : [-1,6],
        'ylim' : [-1,6]}

fitdict = {'color':'g-r',
           'Nmatch':22034}#max 22034

filter_dict = {'257258':'u-g',
               '258259':'g-r',
               '259260':'r-i',
               '260261':'i-z',
               '257':'u',
               '258':'g',
               '259':'r',
               '260':'i',
               '261':'z'}

#ages you want
Etypes = ascii.read('/Users/Christina/GradCourses/Galaxies/hw2/Etypes_AstroGrad.dat')
ages_then = np.arange(0.01,6.,0.5)

#redshifts you want
#use round(z,1), step size of 0.025 to match bc outputs
zrange = np.arange(round(0.5,1),round(6.5,1),0.1)

#SVA1 RedMagic
red = fits.open('redmagic_sva1_public_v6.3_bright.fits.gz')[1].data

rmd = {'g-r':red['mag_auto_g']-red['mag_auto_r'],
       'r-i':red['mag_auto_r']-red['mag_auto_i'],
       'i-z':red['mag_auto_i']-red['mag_auto_z'],
       'g':red['mag_auto_g'],
       'r':red['mag_auto_r'],
       'i':red['mag_auto_i'],
       'z':red['mag_auto_z'],
       'redshift':red['zredmagic']}

#DES throughputs
tputtabs = {'u': ascii.read('/Users/Christina/DES/DES_tput_Ang_u.txt'),
            'g': ascii.read('/Users/Christina/DES/DES_tput_Ang_g.txt'),
            'r': ascii.read('/Users/Christina/DES/DES_tput_Ang_r.txt'),
            'i': ascii.read('/Users/Christina/DES/DES_tput_Ang_i.txt'),
            'z': ascii.read('/Users/Christina/DES/DES_tput_Ang_z.txt'),
            'Y': ascii.read('/Users/Christina/DES/DES_tput_Ang_Y.txt')}
tput = {'u' : {'wave': np.array(tputtabs['u']['wavelength']),
               'tput': np.array(tputtabs['u']['u'])},
        'g' : {'wave': np.array(tputtabs['g']['wavelength']),
               'tput': np.array(tputtabs['g']['g'])},
        'r' : {'wave': np.array(tputtabs['r']['wavelength']),
               'tput': np.array(tputtabs['r']['r'])},
        'i' : {'wave': np.array(tputtabs['i']['wavelength']),
               'tput': np.array(tputtabs['i']['i'])},
        'z' : {'wave': np.array(tputtabs['z']['wavelength']),
               'tput': np.array(tputtabs['z']['z'])},
        'Y' : {'wave': np.array(tputtabs['Y']['wavelength']),
               'tput': np.array(tputtabs['Y']['Y'])}}

def uage(z):
    #age of universe at redshift z
    return integ.quad(ageintegral,z,np.inf)[0]

def ageintegral(z):
    return 1./(H0*np.sqrt(Om*(1.+z)**5. + Or*(1.+z)**6. + Ol*(1+z)))

def I(x):
    return 1./np.sqrt((1+x)**3.*Om+(1+x)**4.*Or+Ol)

def dLum(z):
    dp = cH * integ.quad(I,0,z)[0] #Mpc
    dL = dp * (1+z) * 3.086e24 #cm
    return dL

def calcspectrum(outdir,model,age):
    runfile = outdir+'run_{}_gpl.dat'.format(model.split('.')[0])
    
    if not os.path.exists('{}{}/'.format(outdir,age)):
        os.system('mkdir {}{}'.format(outdir,age))

    specfile = os.path.splitext(model)[0]+'.sed'

    if os.path.exists('{}{}/{}'.format(outdir,age,specfile)):
        print "Spectrum {} has previously been calculated. Skipping galaxevpl.".format(specfile)

    else:
        outfile = open(runfile,'w')
        outfile.write("""%s
%s

%s
"""%(outdir+model,age,specfile))
        
        outfile.close()
        
        os.system('src/galaxevpl < {}'.format(runfile))
        os.system('mv {} {}{}'.format(specfile,outdir,age))
    return '{}{}/{}'.format(outdir,age,specfile)

def colormatch(d):
    tally = {}
    for mod in d['model']:
        tally[mod]=0

    mindict = {}
    for l in d.keys():
        mindict[l] = []
    del mindict['title']

    for i in range(fitdict['Nmatch']):
        diffs = []
        #get minimum of list of differences between model points and data point
        minarray = np.abs(np.array(d[fitdict['color']])-rmd[fitdict['color']][i])
        mindiff = np.nanmin(minarray)
        index = np.where(np.abs(np.array(d[fitdict['color']])-rmd[fitdict['color']][i])==mindiff)[0][0]
        #save model, age, redshift, etc
        for k in mindict.keys():
            mindict[k].append(d[k][index])
        #tally the model
        tally[d['model'][index]]+=1

    print "Number of color matches for each model:"
    for key in tally.keys():
        print key, ": ", tally[key]

    mindict['title']='colorfit_'+d['title']
    plotfromdict(mindict)

def getcolors(specfile,age,fdict):
    #uses global vars zrange and filts
    spectrum = ascii.read(specfile)
    model = specfile.split('/')[-1].split('.')[0]
    #units?
    lam = np.array(spectrum['col1']) #Angstrom
    flux = np.array(spectrum['col2']) * Lsun #ergs/s/Angstrom

    check=0
    fset = {}

    for combo in filtcoms:
        f1,f2 = combo
        filters = str(f1)+str(f2)

        for f in combo:
            if f not in fset.keys():
                fset[f]=[]

                for z in zrange:
                    if age > uage(z):
                        continue

                    mag = getmag(lam, flux, tput[filter_dict[str(f)]]['wave'],
                                 tput[filter_dict[str(f)]]['tput'], z)

                    for M in mass:
                        fset[f].append(mag-2.5*np.log10(M))
                        if check==0:
                            fdict['age (Gyr)'].append(float(age))
                            fdict['redshift'].append(z)
                            fdict['mass'].append(M)
                            fdict['model'].append(model)
                fdict[filter_dict[str(f)]] += fset[f]
                check=1
        color = [fset[f1][i]-fset[f2][i] for i in range(len(mass)*len(zrange))]
        fdict[filter_dict[filters]] += [c for c in color]

    return fdict


def getmag(wave, specflux, filtwave, filttput, redshift):
    #change units from Angstrom to Hz
    nu = cA/np.array(wave)
    nufilt = cA/np.array(filtwave)
    #make sure frequencies can be redshifted in interp
    znurange = np.where((nu.min()<=(nu*(1+redshift))) &
                   ((nu*(1+redshift))<=nu.max()))
    #spectrum flux in frequency units 
    nuflux = (np.array(wave)**2.) * np.array(specflux) / cA**2.

    #filter throughput vs. frequency
    fspec = interpolate.interp1d(nu,nuflux)
    specz = fspec(nu[znurange]*(1+redshift))
    #cut down frequency to filter
    filtrange = np.where((nufilt.min() <= nu[znurange]) & (nu[znurange] <= nufilt.max()))     
    nufinal = nu[znurange][filtrange]
    #interpolate and resample throughput at spectrum frequencies
    ftput = interpolate.interp1d(nufilt, filttput)

    #sort arrays for proper trapz integration
    zipp = zip(nufinal,specz[filtrange])
    zipp.sort()
    nufinal,speczfinal = zip(*zipp)

    itput = ftput(nufinal)
    filtspec = speczfinal * itput

    #calculate total observed flux
    totflux = np.trapz(filtspec * (1+redshift) / nufinal, nufinal) / (4.*np.pi*dLum(redshift)**2) / np.trapz(itput*ABflux/nufinal,nufinal)

    mag = -2.5*np.log10(totflux)#-48.60 #AB

    return mag

def makecolordict(basefile):
    fdict={'mass':[],'age (Gyr)':[],'redshift':[],'model':[]}
    for combo in filtcoms:
        f1,f2 = combo
        filters = str(f1)+str(f2)
        fdict[filter_dict[filters]]=[]
        
        for f in combo:
            if str(f) not in fdict.keys():
                fdict[filter_dict[str(f)]]=[]

    fdict['title']=basefile[7:]+' '
    return fdict

def mergedicts(dicts):
    #dict elements must be lists, not arrays
    new = dicts[0].copy()
    for d in dicts[1:]:
        for k in new.keys():
            new[k] = new[k]+ d[k]
    return new

def plotfromdict(d):                            
    colormap = plt.get_cmap('jet')
    if pltd['redmagic']:
        plt.scatter(rmd[pltd['x']],rmd[pltd['y']],c='m',edgecolor='none',s=5.)


    #magnitude cuts
#    magcut = np.where(np.array(d[pltd['cutmag']])>pltd['maglimitlow'])
    plt.scatter(np.array(d[pltd['x']]),
                np.array(d[pltd['y']]),
                c=np.array(d[pltd['color']]),cmap=colormap,edgecolor='none',s=6.)
    plt.xlim(pltd['xlim'])
    plt.ylim(pltd['ylim'])
    plt.colorbar(label=pltd['color'])
    plt.xlabel(pltd['x'])
    plt.ylabel(pltd['y'])
    plt.title(d['title'])
    savefile = pltd['plotdir']+d['title'].replace(' ','_')+pltd['y']+'_vs_'+pltd['x']
    print "Saving plot: ", savefile
    plt.savefig(savefile)
    plt.close()                                  


def main():

    cspcmd = ["python","run_csp.py","--IMF",IMF,"--model",model,
              "--SFhistory"]
    for sf in starform:
        cspcmd.append(sf)
    if 'CONST' in starform:
        cspcmd.append("--constSFR")
        for sfr in constSFR:
            cspcmd.append(sfr)
    if 'BURST' in starform:
        cspcmd.append("--burstlength")
        for l in burstlength:
            cspcmd.append(l)
    cspcmd.append("--metal")
    for metal in metals:
        cspcmd.append(metal)
    if dust=='dust':
        dustdir='dust/tauV{}mu{}'.format(tauV,mu).replace('.','')
        cspcmd.append("-dust")
        cspcmd.append("--tauV")
        cspcmd.append(tauV)
        cspcmd.append("--mu")
        cspcmd.append(mu)
    else:
        dustdir='nodust'

    call(cspcmd)

    #calculate colors at different redshifts
    alldicts = []
    for metal in metals:
        for sf in starform:
            if sf=='CONST':
                for sfr in constSFR:
                    sfdir = 'CONST/{}Msun-yr'.format(sfr)
                    cfiledir = './output/{}/{}/{}/{}/'.format(model,dustdir,metal,sfdir)
                    basefile = 'bc2003_hr_{}_{}_{}_{}_{}'.format(metal,IMF[:4],dust,sf,sfr)
                    fdict = makecolordict(basefile)
                    for a in ages_then:
                        sfile = calcspectrum(cfiledir,basefile+'.ised',a)
                        fdict = getcolors(sfile,a,fdict)

#get K-corrected spectrum!!!!!!

                    alldicts.append(fdict)
            elif sf=='BURST':
                for l in burstlength:
                    sfdir = 'BURST/{}Gyr'.format(l)
                    cfiledir = './output/{}/{}/{}/{}/'.format(model,dustdir,metal,sfdir)
                    basefile = 'bc2003_hr_{}_{}_{}_{}_{}'.format(metal,IMF[:4],dust,sf,l)
                    fdict = makecolordict(basefile)
                    for a in ages_then:
                        sfile = calcspectrum(cfiledir,basefile+'.ised',a)
                        fdict = getcolors(sfile,a,fdict)
                    alldicts.append(fdict)

            else:
                    cfiledir = './output/{}/{}/{}/{}/'.format(model,dustdir,metal,sf)
                    basefile = 'bc2003_hr_{}_{}_{}_{}'.format(metal,IMF[:4],dust,sf)
                    fdict = makecolordict(basefile)
                    for a in ages_then:
                        sfile = calcspectrum(cfiledir,basefile+'.ised',a)
                        fdict = getcolors(sfile,a,fdict)
                    alldicts.append(fdict)
 
    for m in alldicts:
        plotfromdict(m)
    d = mergedicts(alldicts)
    d['title']='combinedmodel'
    "Plotting the combined model."
    plotfromdict(d)

    if colorfit:
        colormatch(d)

    #sanity checks
    #print "age of the universe = ", age(0)*3.17e-17,  " Gigayears ?"

if __name__=="__main__":
    main()
