import numpy as np
import scipy.integrate as integ
from subprocess import call

#constants
#flat universe
h=0.72
H0=h*100*3.24e-20 #1/s
Om = 0.27
Ol = 0.73
Or = 4.2e-5/h**2.

IMF = "salpeter"
model = "Padova1994"
dust = 'nodust'
metals = ["m42"]
starform = ["CONST"]
constSFR = ["10"]
filts = [(257,258),(258,259),(259,260)]

#files
#output = mag_dir + model_dir + IMF_dir + age_dir + dust_dir

#age directory = output/model/dust/metal/sf/agenow/...color_c1_c2


#ages you want
ages_then = np.arange(0.5,0.6,0.1)

#redshifts you want
zrange = np.arange(1.,1.1,0.1)

def age(z):
    return integ.quad(ageintegral,z,np.inf)[0]

def ageintegral(z):
    return 1./(H0*np.sqrt(Om*(1.+z)**5. + Or*(1.+z)**6. + Ol*(1+z)))

def main():
    cspcmd = ["python","run_csp.py","--IMF",IMF,"--model",model,
              "--SFhistory"]
    for sf in starform:
        cspcmd.append(sf)
    if sf=='CONST':
        cspcmd.append("--constSFR")
        for sfr in constSFR:
            cspcmd.append(sfr)
    cspcmd.append("--metal")
    for metal in metals:
        cspcmd.append(metal)

    call(cspcmd)

    cmevcmd = ["python","run_cm_evo.py","--IMF",IMF,"--model",model,
               "--H0",str(h*100.),"--Omega_m",str(Om),"--Omega_l",str(Ol),
               "--SFhistory"]
    for sf in starform:
        cmevcmd.append(sf)
        if sf=='CONST':
            cmevcmd.append("--constSFR")
            for sfr in constSFR:
                cmevcmd.append(sfr)
    cmevcmd.append("--metal")
    for metal in metals:
        cmevcmd.append(metal)
    cmevcmd.append("--filtercombo")
    for combo in filts:
        cmevcmd.append(str(combo[0]))
        cmevcmd.append(str(combo[1]))
    cmevcmd.append("--age")
    for a in ages_then:
        for z in zrange:
        #calculate required input age in Gyr 
        #    = age at z=0 for age to be correct at z
            age_now = (age(0)-age(z))*3.17e-17 + a 
            agez.append([age_now,z])
            cmevcmd.append(str(age_now))

    call(cmevcmd)

    #plot some colors
    for metal in metals:
        for sf in starform:
            if sf=='CONST':
                for sfr in constSFR:
                    sfdir = 'CONST/{}Msun-yr'.format(sfr)
                    for a,z in agez:
                        for f1,f2 in filts:
                            cfile = './output/{}/{}/{}/{}/{:10.1f}/bc2003_hr_{}_{}_ssp_{}_{}.color_{}_{}'.format(model,dust,metal,sfdir,a,metal,IMF[:4],sf,sfr,f1,f2)
            else:
                for a,z in agez:
                    for f1,f2 in filts:
                        cfile = './output/{}/{}/{}/{}/{:10.1f}/bc2003_hr_{}_{}_ssp_{}.color_{}_{}'.format(model,dust,metal,sf,a,metal,IMF[:4],sf,f1,f2)
    #sanity checks
    #print "age of the universe = ", age(0)*3.17e-17,  " Gigayears ?"

#select values at z row

if __name__=="__main__":
    main()
