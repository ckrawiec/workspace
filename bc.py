import os
import numpy as np
import scipy.integrate as integ
from scipy import interpolate
import matplotlib.pyplot as plt
from subprocess import call
from astropy.io import ascii, fits

"""
For redshift range: use round(z,1), step size of 0.025 to match bc outputs
"""

title = "SSP and 1-3Gyr BURSTs with Z=0.02,0.04, Padova2000 tracks, Chabrier IMF"
ages_then = np.arange(1.,14.,0.5)
zrange = np.arange(round(0.15,2),round(0.95,2),0.005)
IMF = "chabrier"
model = "Padova2000"
dusts = ['nodust'] #'dust' or 'nodust'
tauV="1."
mu="0.3"
metals = ["m162","m172"]
starform = ["SSP","BURST"]
mass = np.logspace(11.9,13.,5)
constSFR = ["10","20","30","40","50"]
burstlengths = ["1","2","3"]
filtcoms = [(258,259),(259,260),(260,261)]

pltd = {'x' : 'redshift',
        'y' : 'i-z',
        'color' : 'age (Gyr)',
        'cutmag' : 'z',
        'maglimitlow' : 0.,
        'plotdir' : 'plots/SVA1RedMagic',
        'redmagic' : True,
        'xlim' : [0.15,0.85],
        'ylim' : [0.2,1.5]}

colorfit = True
fitdict = {'color':'g-r',
           'Nmatch':22034}#max 22034

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
#1. #Lsun (for easier Mass units) 
 
#if you want to compare with SDSS E types
#Etypes = ascii.read('/Users/Christina/GradCourses/Galaxies/hw2/Etypes_AstroGrad.dat')

#dictionary for parsing BC03 color codes
filter_dict = {'257258':'u-g',
               '258259':'g-r',
               '259260':'r-i',
               '260261':'i-z',
               '257':'u',
               '258':'g',
               '259':'r',
               '260':'i',
               '261':'z'}

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

class Template:
    """
    model = Padova1994 or Padova2000
    imf = salpeter or chabrier
    metal = m42,m52,m62,m72 (Padova1994) or m142,m152,m162,m172 (Padova2000)
    dust = dust or nodust
    sf = CONST, BURST, SSP, EXP1, EXP2, EXP3
    sfr = star formation rate for CONST (integer, Msun/yr)
    burstlength = length of starburst for BURST (integer, Gyr)
    tau = BC2003 tau for dust (1.0 in paper)
    mu = BC2003 mu for dust (0.3 in paper)
    """
    def __init__(self,model,imf,metal,dust,sf,sfr='',burstlength='',tau='',mu=''):
        self.modelname = 'bc2003_hr_{}_{}_{}_{}'.format(metal,imf[:4],dust,sf)

        self.dustdir = dust
        if tau:
            self.dustdir = os.path.join(dust,'tauV{}mu{}'.format(tauV,mu)).replace('.','')

        self.sfdir = sf
        if sfr:
            self.sfdir = os.path.join(sf,sfr+'Msun-yr')
            self.modelname+='_'+sfr
        if burstlength:
            self.sfdir = os.path.join(sf,burstlength+'Gyr')
            self.modelname+='_'+burstlength
        self.outdir = os.path.join('output',model,self.dustdir,metal,self.sfdir)            
        
    def calcspectrum(self,age):
        """
        Runs galaxevpl for the model at the given age and plots the spectrum.
        Returns wavelength and flux arrays in Angstroms and ergs/s/Angstroms.
        """
        age = str(age)
    
        if not os.path.exists(os.path.join(self.outdir,age)):
            os.system('mkdir {}/{}'.format(self.outdir,age))

        specfile = self.modelname+'.sed'

        #if os.path.exists('{}/{}/{}'.format(self.outdir,age,specfile)):
            #print "Spectrum {} at age {} Gyr has previously been calculated. Skipping galaxevpl.".format(specfile,age)
        #else:
        runfile = '{}/run_{}_gpl.dat'.format(self.outdir,self.modelname)
        outfile = open(runfile,'w')
        outfile.write("""%s
%s

%s
"""%(os.path.join(self.outdir,self.modelname),age,specfile))
        outfile.close()
        
        os.system('src/galaxevpl < {}'.format(runfile))
        os.system('mv {} {}/{}'.format(specfile,self.outdir,age))
    
        spectrum = ascii.read(os.path.join(self.outdir,age,specfile))
        lam = np.array(spectrum['col1']) #Angstrom
        flux = np.array(spectrum['col2']) * Lsun #ergs/s/Angstrom
        plt.plot(lam,flux)
        plt.xlabel('wavelength (Angstrom)')
        plt.ylabel('flux (ergs/s/Angstrom)')
        plt.xlim(2000,12000)
        plt.ylim(0.,1.5e31)
        plt.savefig('plots/spectra/'+self.modelname+age.replace('.','-'))
        plt.close()
        
        return lam, flux

def uage(z):
    #age of universe at redshift z in Gyr
    return integ.quad(ageintegral,z,np.inf)[0] / 3.2e16

def ageintegral(z):
    return 1./(H0*np.sqrt(Om*(1.+z)**5. + Or*(1.+z)**6. + Ol*(1+z)))

def I(x):
    return 1./np.sqrt((1+x)**3.*Om+(1+x)**4.*Or+Ol)

def dLum(z):
    dp = cH * integ.quad(I,0,z)[0] #Mpc
    dL = dp * (1+z) * 3.086e24 #cm
    return dL

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

def getcolors(templates,age):
    """
    For a given list of templates at the given age, returns a dictionary
    with magnitude, colors, models, ages, masses, and redshifts listed.
    
    Uses global vars zrange, mass, and filts.
    """

    fdict={'mass':[],'age (Gyr)':[],'redshift':[],'model':[]}
    for combo in filtcoms:
        f1,f2 = combo
        filters = str(f1)+str(f2)
        fdict[filter_dict[filters]]=[]
        
        for f in combo:
            if str(f) not in fdict.keys():
                fdict[filter_dict[str(f)]]=[]
    
    for t in templates:
        wave,flux = t.calcspectrum(age)
        
        fset = {}
        
        check=0
        for combo in filtcoms:
            f1,f2 = combo
            filters = str(f1)+str(f2)

            for f in combo:
                if f not in fset.keys():
                    fset[f]=[]

                    for z in zrange:
                        if age > uage(z):
                            continue

                        mag = getmag(wave, flux, f, z)

                        for M in mass:
                            fset[f].append(mag-2.5*np.log10(M))
                            if check==0:
                                fdict['age (Gyr)'].append(float(age))
                                fdict['redshift'].append(z)
                                fdict['mass'].append(M)
                                fdict['model'].append(t.modelname)
                    fdict[filter_dict[str(f)]] += fset[f]
                    check=1
            color = [fset[f1][i]-fset[f2][i] for i in range(len(fset[f1]))]
            fdict[filter_dict[filters]] += [c for c in color]

    return fdict

def getmag(wave, specflux, filter, redshift):
    filtwave = tput[filter_dict[str(filter)]]['wave']
    filttput = tput[filter_dict[str(filter)]]['tput']

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
    if redshift==0.:
        d = 10. * 3.086e18 #cm
    else:
        d = dLum(redshift) #cm
    totflux = np.trapz(filtspec * (1+redshift) / nufinal, nufinal) / (4.*np.pi*d**2) / np.trapz(itput*ABflux/nufinal,nufinal)

    mag = -2.5*np.log10(totflux)#-48.60 #AB

    return mag

def mergedicts(dicts):
    #dict elements must be lists, not arrays
    new = dicts[0].copy()
    for d in dicts[1:]:
        for k in new.keys():
            new[k] = new[k]+ d[k]
    return new

def plotfromdict(d):                            
#    for k in d.keys():
#        print k, len(d[k])
    colormap = plt.get_cmap('jet')
    if pltd['redmagic']:
        plt.scatter(rmd[pltd['x']],rmd[pltd['y']],c='m',edgecolor='none',s=5.)

    #magnitude cuts
    #magcut = np.where(np.array(d[pltd['cutmag']])>pltd['maglimitlow'])
    plt.scatter(np.array(d[pltd['x']]),
                np.array(d[pltd['y']]),
                c=np.array(d[pltd['color']]),cmap=colormap,edgecolor='none',s=6.)
    plt.xlim(pltd['xlim'])
    plt.ylim(pltd['ylim'])
    plt.colorbar(label=pltd['color'])
    plt.xlabel(pltd['x'])
    plt.ylabel(pltd['y'])
    plt.title(d['title'])
    savefile = "{}/{}_{}_vs_{}".format(pltd['plotdir'],d['title'].replace(' ','').replace(',','_').replace('.','-'),
                                       pltd['y'],pltd['x'])
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
        for l in burstlengths:
            cspcmd.append(l)

    cspcmd.append("--metal")
    for metal in metals:
        cspcmd.append(metal)

    for dust in dusts:
        if dust=='nodust':
            call(cspcmd)
        else:
            dustcmd  = ["-dust","--tauV"]
            for tau in tauV:
                dustcmd.append(tauV)
            for m in mu:
                dustcmd.append("--mu")
            cspcmd.append(mu)
            call(cspcmd+dustcmd)                        

    #calculate colors at different ages and redshifts
    templates = []
    alldicts = []
    for metal in metals:
        for du in dusts:
            for sf in starform:
                if sf=='CONST':
                    for sfr in constSFR:
                        t = Template(model,IMF,metal,du,sf,
                                     sfr=sfr)
                        templates.append(t)
                elif sf=='BURST':
                    for l in burstlengths:
                        t = Template(model,IMF,metal,du,sf,
                                     burstlength=l)
                        templates.append(t)
                else:
                    t = Template(model,IMF,metal,du,sf)
                    templates.append(t)

    for a in ages_then:
        fdict = getcolors(templates,a)
        alldicts.append(fdict) 
#    for m in alldicts:
#        plotfromdict(m)
    d = mergedicts(alldicts)
    if title:
        d['title']=title
    else:
        d['title']='combinedmodel'
    print "Plotting the combined model."
    plotfromdict(d)

    if colorfit:
        colormatch(d)

    #sanity checks
    #print "age of the universe = ", age(0)*3.17e-17,  " Gigayears ?"

if __name__=="__main__":
    main()
