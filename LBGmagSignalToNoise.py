import matplotlib.pyplot as plt
import numpy as np
import itertools
from astropy.io import fits
import matplotlib
from matplotlib import rcParams
from scipy import integrate, interpolate

#constants
c = 3.e5 #km/s
G = 4.3e-9 #km^2 * Mpc * MSun^-1 * s^-2
#cosmology
h = 0.7
H0 = h*100. #km/s/Mpc
Om = 0.3
Ol = 0.7

#data
#density of Y1A1 "High Density Catalog" (0.1<z<0.7) (reliable out to 0.65) ~1.0e-3 Mpc^-3
#same for redmagic Faint in SVA1, according to Rozo et al 2015 (1.0 1/h /Mpc^3) Table 3
Nlbg = [1000.,5000.,10000.,20000.]
red_density = 5.6e-5 # per sq. arcsec, rough estimate (rounded low)

angles_deg = np.logspace(-2,0)
zlrange = [0.2,0.5,0.5,0.2,0.8]
zsrange = [0.5,1.0,2.,4.,4.]

#Figure 8 of Clampitt et al 2016
Mlens = 2.e13
    
def main():
    plotNpairs()
    #plotmu()
    printcalc()
    plotSNR()
    
def H(z):
    return H0*np.sqrt(Om*(1+z)**3. + Ol)

def dcomI(z):
    return c/H(z)

def Dang(z1,z2):
    #Mpc
    dcom1 = integrate.quad(dcomI,0.,z1)[0]
    dcom2 = integrate.quad(dcomI,0.,z2)[0]
    return (1./(1+z2)) * (dcom2 - dcom1)

def rhocrit(z):
    # MSun/Mpc^3
    return 3 * H(z)**2. / (8*np.pi*G)

def concentration(M,z):
    #Duffy et al 2008
    A = 5.71
    B = -0.084
    C = -0.47
    Mpivot = 2.e12/h #Msun
    return A * (M/Mpivot)**B * (1+z)**C

def deltac(conc):
    num = conc**3.
    denom = np.log(1.+conc) - conc/(1.+conc)
    return (200./3.) * num / denom

def r200(M200,z):
    #Mpc
    return ((3.*M200) / (200.*4.*np.pi * rhocrit(z)))**(1./3.)

def rs(conc,r200):
    #scale radius, Mpc
    return r200/conc

def rhoNFW(R, z, dc, rs):
    return (dc * rhocrit(z)) / ((R/rs) * (1.+R/rs)**2.)

def sigmaNFW(R, z, dc, rs):
    #Msun/Mpc^2
    x = R/rs
    pc = rhocrit(z)
    coeff = 2.*rs*dc*pc / (x**2.-1.)
    if x < 1:
        return coeff * (1.-(2./np.sqrt(1.-x**2.))*np.arctanh(np.sqrt((1.-x)/(1.+x))))
    elif x == 1:
        return 2.*rs*dc*pc / 3.
    elif x > 1:
        return coeff * (1.-(2./np.sqrt(x**2.-1.))*np.arctan(np.sqrt((x-1.)/(1.+x))))
    
def kappaNFW(Ml,theta,zlens,zsrc):
    #theta in radians
    conc = concentration(Ml,zlens)
    dc = deltac(conc)
    rscale = rs(conc,r200(Ml,zlens))
    Ds = Dang(0.,zsrc)
    Dl = Dang(0.,zlens)
    Dls = Dang(zlens,zsrc)
    sigmacrit = (c**2./(4*np.pi*G)) * Ds/(Dl*Dls)
    return sigmaNFW(theta*Dl,zlens,dc,rscale)/sigmacrit

def plotNpairs():
    for N in Nlbg:
        npairs = []
        angles_rad = degtorad(angles_deg)
        for angle in angles_rad:
            area = annulus_area(angle, angle/2.)
            nred = red_density * area
            npairs.append(nred * N)
        plt.plot(angles_deg, npairs, 'o-', label='$N_{LBG}=$'+str(int(N)))
    plt.xlabel('$\\theta (\deg)$')
    plt.xscale('log')
    plt.ylabel('number of redmagic/LBG pairs')
    plt.legend(loc='best')
    plt.savefig('redmagic_lbg_Npairs')
    plt.close()

def plotDang():
    #to compare with David W. Hogg (2000) Figure 2 
    zrange = np.arange(0,5,0.01)
    y = np.array([Dang(0,zi)*H0/c for zi in zrange])
    plt.plot(zrange, y)
    plt.xlabel('z')
    plt.ylabel('dimensionless angular diameter distance ($D_{A}  H_{0} / c$)')
    plt.title('$\Omega_{M}=0.3, \Omega_{\Lambda}=0.7$')
    plt.savefig('Dang_Hogg2000')
    plt.close()

def radtoarc(theta):
    return theta * 206265.
def arctorad(theta):
    return theta / 206265.
def arctodeg(theta):
    return theta / 3600.
def degtorad(theta):
    return theta * np.pi/180
def radtodeg(theta):
    return theta * 180./np.pi
def dtheta(theta):
    return 0.5*theta
def steradtosqarcsec(omega):
    return omega * (180./np.pi)**2. * 1.296e7
def sqdegreetosqarcsec(omega):
    return omega*1.296e7

def annulus_area(theta, dtheta):
    #theta, dtheta in radians
    #returns solid angle in arcsec**2.
    Asterad = 2*np.pi * (np.cos(theta) - np.cos(theta+dtheta))
    return steradtosqarcsec(Asterad)

def muest(n,n0,alpha):
    return 1. + (n/n0 - 1.)/(alpha - 1.)

def Npairs(theta, Nsrc):
    #theta in radians
    return red_density * annulus_area(theta, dtheta(theta)) * Nsrc

def SNR(theta, Nsrc, M, alph, zl, zs):
    return 2.*kappaNFW(M, theta, zl, zs) * (alph-1.) * np.sqrt(Npairs(theta, Nsrc))

def plotSNR():
    SNRs = []
    angles_rad = degtorad(angles_deg)
    alphas = [1.1,3.]
    for zl,zs in zip(zlrange, zsrange):
        clist = rcParams['axes.color_cycle']
        cgen = itertools.cycle(clist)
        for N in Nlbg[::-1]:
            SNRs_top = [SNR(angle, N, Mlens, np.max(alphas), zl, zs)
                        for angle in angles_rad]
            SNRs_bot = [SNR(angle, N, Mlens, np.min(alphas), zl, zs)
                        for angle in angles_rad]
            plt.fill_between(angles_deg, SNRs_bot, SNRs_top, alpha=0.4, label='Nsrc='+str(int(N)), facecolor=cgen.next())
        plt.xlabel('angle (deg)')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel('S/N')
        plt.legend()
        plt.title('Mlens={} zlens={} zsrc={}, alpha={}-{}'.format(Mlens, zl, zs, np.min(alphas), np.max(alphas)))
        plt.savefig('SNR_angle_zl{}_zs{}.png'.format(zl,zs))
        plt.close()
        
def varmuest(theta,n0,alpha):
    #theta in deg
    var = 1. / (alpha-1)**2. / (n0*2.*np.pi*theta*dtheta(theta))
    return var

#From Hildebrandt et al 2009
def lumfunc(M, phistar, Mstar, alphaLF):
    return 0.4*phistar*np.log(10.) * (10.**(0.4*(Mstar-M)))**(alphaLF+1.) * np.exp(-10.**(0.4*(Mstar-M)))

def linfit(x,a,b):
    return a*x + b

def quadfit(x,a,b,c):
    return a*x**2 + b*x + c

def cubicfit(x,a,b,c,d):
    return a*x**3 + b*x**2 + c*x + d

def Hild09interp():
    #estimates from figure - g-dropouts
    n = [0.003,0.02,0.06,0.2,0.4,0.7,0.9,0.6,0.12]
    mi = [23.25,23.75,24.25,24.75,25.25,25.75,26.25,26.75,27.25]

    newmi = np.arange(23.5,24.5,0.1)
    interp = np.polyfit(mi[:3],np.log10(n[:3]),1)
    plt.semilogy(newmi,10**(linfit(newmi,*interp)))
    plt.semilogy(mi,n,'-')
    plt.xlabel('m_i_AB')
    plt.ylabel('N/arcmin^2')
    plt.show()
    
    #bin from 23.5 to 24
    Npersqdeg = np.trapz(newmi,10**(linfit(newmi,*interp))) * (60.)**2.
    lbgalpha = 2.5*interp[0]
    print Npersqdeg, lbgalpha

#integral from m1 to mlim dm/sigma^2 = (alpha-1)^2 * dN/dm
def plotmu():
    angles_rad = degtorad(angles_deg)
    for zl,zs in zip(zlrange,zsrange):
        conc = concentration(Mlens,zl)
        dconc = deltac(conc)
        rscale = rs(conc,r200(Mlens,zl))
        Dl = Dang(0.,zl)
        kappas = np.array([kappaNFW(angle,zl,zs,dconc,rscale) for angle in angles_rad])
        mus = 1 + 2.*kappas
        plt.semilogx(angles_deg, mus, label="zs={},zl={}".format(zs,zl))
    plt.legend(loc='best')
    plt.xlabel('angle (deg)')
    plt.ylabel('magnification')
    plt.savefig('magnification_vs_angle_zbins')
    plt.close()

#Change lens masses, number densities, and alphas, check alphas from Hild
#gglensing 0.1-10 Mpc/h

def printcalc():
    print "M = ", Mlens, "Msun"
    print "r200(z={}) = {} Mpc".format(0.8, r200(Mlens,0.8))
    print "r_scale = ", rs(concentration(Mlens,0.8), r200(Mlens,0.8)), "Mpc"
    print "concentration = ", concentration(Mlens,0.8)
    print "delta_concentration = ", deltac(concentration(Mlens,0.8))
    print "comoving distance at z=0.1: ", integrate.quad(dcomI, 0, 0.1)[0]
    print "comoving volume at z=0.1: ", 4*np.pi/3. * integrate.quad(dcomI, 0, 0.1)[0]
    print "comoving distance at z=0.7: ", integrate.quad(dcomI, 0, 0.7)[0]
    print "0.1-10 Mpc/h = {}-{} degrees at z=0.4".format(0.1/h/Dang(0,0.4), 10/h/Dang(0,0.4))
    
#for a given size lens
#for each magnitude bin
#get number density
#get alpha
#for each annulus
#calculate estimator using weights, etc

#need alphas ----- luminosity function for LBGs/source > z=2.5 (fit Schecter to data?)
#could use Manuel's just to see.. ours should be better than that
#number densities
#plot bands for different z range or masses, or alphas

if __name__=="__main__":
    main()
