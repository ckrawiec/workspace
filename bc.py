import os
import numpy as np
import math
import scipy.integrate as integ
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

colorfit = True

IMF = "chabrier"
model = "Padova2000"
dust = 'nodust' #'dust' or 'nodust'
tauV="1."
mu="0.3"
metals = ["m162","m172"]
starform = ["BURST","SSP","EXP1"]
mass = np.logspace(10.9,12.,5)
constSFR = ["10"]#,"20","30","40","50"]
burstlength = ["2"]
filtcoms = [(258,259),(260,261)]

pltd = {'x' : 'r',
        'y' : 'g-r',
        'color' : 'age (Gyr)',
        'cutmag' : 'z',
        'maglimitlow' : 0.,
        'plotdir' : 'plots/SVA1RedMagic/{}/'.format(dust),
        'redmagic' : True,
        'xlim' : [16,24],
        'ylim' : [0.9,2.2]}

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
ages_then = np.arange(1.,10.,1.)

#redshifts you want
#use round(z,1), step size of 0.025 to match bc outputs
zrange = np.arange(round(0.2,1),round(0.8,1),0.025)

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

def uage(z):
    #age of universe at redshift z
    return integ.quad(ageintegral,z,np.inf)[0]

def ageintegral(z):
    return 1./(H0*np.sqrt(Om*(1.+z)**5. + Or*(1.+z)**6. + Ol*(1+z)))

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


def getcolors(agedict,cfiledir,basefile):
    #uses global vars zrange and filts
    #output file columns: M_ne_AB:10,M_ev_AB:11,m_ev_AB:12,dm:4
    fdict={'mass':[],'age (Gyr)':[],'redshift':[],'model':[]}
    check=0    
    for combo in filtcoms:
        f1,f2 = combo
        filters = str(f1)+str(f2)
        fdict[filter_dict[filters]]=[]

        for f in combo:
            if str(f) not in fdict.keys():
                fdict[filter_dict[str(f)]]=[]

                for age in agedict.keys():
                    if len(agedict[age]['ages_now,z'])==0:
                        continue
                    zipped=agedict[age]['ages_now,z']
                    ages,zs = zip(*zipped)
                        
                    for a in ages:
                        mfile = cfiledir+'{:10.1f}/'.format(float(a)).strip()+'{}.magnitude_F{}'.format(basefile,f)
                        mtab = ascii.read(mfile)
                        zcol = mtab['col1']
                        mag = np.array([float(i) for i in mtab['col10']])
                        dm = np.array([float(i) for i in mtab['col4']])
                        for j in range(len(zcol)):
                            for M in mass:
                                fdict[filter_dict[str(f)]].append(mag[j]+dm[j]-2.5*np.log10(M))                        
                                if check==0:
                                    fdict['age (Gyr)'].append(float(age))#(uage(zcol[j])-uage(z))*3.17e-17 + float(age))
                                    fdict['redshift'].append(zcol[j])
                                    fdict['mass'].append(M)
                                    fdict['model'].append(basefile)
            check=1
        for age in agedict.keys():
            if len(agedict[age]['ages_now,z'])==0:
                continue
            zipped=agedict[age]['ages_now,z']
            ages,zs = zip(*zipped)
            for a in ages:
                cfile = cfiledir+'{:10.1f}/'.format(float(a)).strip()+'{}.color_F{}_F{}'.format(basefile,f1,f2)
                ctab = ascii.read(cfile)
                zcol = ctab['col1']
                color = np.array([float(i) for i in ctab['col10']])
                for j in range(len(zcol)):
                    for M in mass:
                        fdict[filter_dict[filters]].append(color[j])


    fdict['title']=os.path.splitext(cfile)[0].split('/')[-1][7:]+' '
    return fdict

def getspectrum(outdir,model,age):
    runfile = outdir+'run_{}_gpl.dat'.format(model.split('.')[0])
    
    if not os.path.exists('{}{}/'.format(outdir,age)):
        os.system('mkdir {}{}'.format(outdir,age))

    specfile = os.path.splitext(model)[0]+'.sed'
    
    outfile = open(runfile,'w')
    outfile.write("""%s
%s

%s
"""%(outdir+model,age,specfile))

    outfile.close()
    
    os.system('src/galaxevpl < {}'.format(runfile))
    os.system('mv {} {}{}').format(specfile,outdir,age)

def plotfromdict(d):                            
    colormap = plt.get_cmap('jet')
    if pltd['redmagic']:
        plt.scatter(rmd[pltd['x']],rmd[pltd['y']],c='m',edgecolor='none')
    #magnitude cuts
    magcut = np.where(np.array(d[pltd['cutmag']])>pltd['maglimitlow'])
    plt.scatter(np.array(d[pltd['x']])[magcut],
                np.array(d[pltd['y']])[magcut],
                c=np.array(d[pltd['color']])[magcut],cmap=colormap,edgecolor='none',s=6.)
    plt.xlim(pltd['xlim'])
    plt.ylim(pltd['ylim'])
    plt.colorbar(label=pltd['color'])
    plt.xlabel(pltd['x'])
    plt.ylabel(pltd['y'])
    plt.title(d['title'])
    print pltd['plotdir']+d['title'].replace(' ','_')+pltd['y']+'_vs_'+pltd['x']
    plt.savefig(pltd['plotdir']+d['title'].replace(' ','_')+pltd['y']+'_vs_'+pltd['x'])
    plt.close()                                  

def mergedicts(dicts):
    #dict elements must be lists, not arrays
    new = dicts[0].copy()
    for d in dicts[1:]:
        for k in new.keys():
            new[k] = new[k]+ d[k]
    return new

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
        cspcmd.append("-dust")
        cspcmd.append("--tauV")
        cspcmd.append(tauV)
        cspcmd.append("--mu")
        cspcmd.append(mu)

    call(cspcmd)

    cmevcmd = ["python","run_cm_evo.py","--IMF",IMF,"--model",model,
               "--H0",str(h*100.),"--Omega_m",str(Om),"--Omega_l",str(Ol),
               "--SFhistory"]
    for sf in starform:
        cmevcmd.append(sf)
    if 'CONST' in starform:
        cmevcmd.append("--constSFR")
        for sfr in constSFR:
            cmevcmd.append(sfr)
    if 'BURST' in starform:
        cmevcmd.append("--burstlength")
        for l in burstlength:
            cmevcmd.append(l)
    cmevcmd.append("--metal")
    for metal in metals:
        cmevcmd.append(metal)
    for combo in filtcoms:
        cmevcmd.append("--filtercombo")
        cmevcmd.append(str(combo[0]))
        cmevcmd.append(str(combo[1]))
    cmevcmd.append("--age")

    age_dict = {}




    for a in ages_then:
        age_dict[str(a)]={'ages_now,z':[]}
        for z in zrange:
        #calculate required input age in Gyr 
        #    = age at z=0 for age to be correct at z
            age_now = (uage(0)-uage(z))*3.17e-17 + a 
            if age_now<=13.3:
                age_dict[str(a)]['ages_now,z'].append([age_now,z])
                cmevcmd.append(str(age_now))
    if dust=='dust':
        cmevcmd.append("-dust")
        cmevcmd.append("--tauV")
        cmevcmd.append(tauV)
        cmevcmd.append("--mu")
        cmevcmd.append(mu)
        dustdir='{}/tauV{}mu{}'.format(dust,tauV,mu).replace('.','')
    else:
        dustdir=dust
    call(cmevcmd)

    #plot some colors
    alldicts = []
    for metal in metals:
        for sf in starform:
            if sf=='CONST':
                for sfr in constSFR:
                    sfdir = 'CONST/{}Msun-yr'.format(sfr)

                    cfiledir = './output/{}/{}/{}/{}/'.format(model,dustdir,metal,sfdir)
                    basefile = 'bc2003_hr_{}_{}_{}_{}_{}'.format(metal,IMF[:4],dust,sf,sfr)
                    for a in ages_then:
                        getspectrum(cfiledir,basefile+'.ised',age)


#                    dsf = gettcolors(age_dict, cfiledir, basefile)
#                    alldicts.append(dsf)
            if sf=='BURST':
                for l in burstlength:
                    sfdir = 'BURST/{}Gyr'.format(l)
                    cfiledir = './output/{}/{}/{}/{}/'.format(model,dustdir,metal,sfdir)
                    basefile = 'bc2003_hr_{}_{}_{}_{}_{}'.format(metal,IMF[:4],dust,sf,l)
                    dsf = getcolors(age_dict, cfiledir, basefile)
                    alldicts.append(dsf)

            else:
                    cfiledir = './output/{}/{}/{}/{}/'.format(model,dustdir,metal,sf)
                    basefile = 'bc2003_hr_{}_{}_{}_{}'.format(metal,IMF[:4],dust,sf)
 #                   dsf = getcolors(age_dict, cfiledir, basefile)
 #                   alldicts.append(dsf)
    for m in alldicts:
        plotfromdict(m)
    d = mergedicts(alldicts)
    "Plotting the combined model."
    plotfromdict(d)

    if colorfit:
        colormatch(d)

    #sanity checks
    #print "age of the universe = ", age(0)*3.17e-17,  " Gigayears ?"

if __name__=="__main__":
    main()
