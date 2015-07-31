import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import re
import sys
import glob

g1=0.02
g2=0.0

n_per_trial = 100e3

const_script = os.path.join('/home','ckrawiec','git','bfd','greatConstant')
moments = glob.glob(os.path.join('/home','ckrawiec','moments','stamps_moments_5*'))

def gplot(garr, errbars, title, name):
    trials = range(1,len(garr)+1)
    plt.scatter(trials, garr)
    plt.errorbar(trials, garr, yerr=errbars, linestyle="None")
    if title=='g1':
        plt.plot(trials, [g1]*len(trials), label='g1='+str(g1))
        plt.plot(trials, [np.mean(garr)]*len(trials), label='mean g1='+str(np.mean(garr)))
        plt.ylim(g1-0.2*g1,g1+0.2*g1)
    if title=='g2':
        plt.plot(trials, [g2]*len(trials), label='g2='+str(g2))
        plt.plot(trials, [np.mean(garr)]*len(trials), label='mean g2='+str(np.mean(garr)))
    plt.xlabel('trial ('+str(n_per_trial)+' targets)')
    plt.ylabel(title)
    plt.legend(loc='upper left')
    plt.savefig(os.path.join('/home','ckrawiec','sims','gplots',name+'_'+title))
    plt.close()

def gmeasure(const_script, moments):
    glist = []
    ulist = []
    used = 0
    for moment in moments:
        out = subprocess.check_output([const_script, moment])
        outre = re.search('(?<=\n).+\n.+\n.+',out)
        print outre.group(0)
        use = re.search('(?<=Use )\S+',out)
        used += int(use.group(0))
        g1re = re.search('(?<=g1: )\S+',out)
        g2re = re.search('(?<=g2: )\S+',out)
        g1 = float(g1re.group(0))
        g2 = float(g2re.group(0))
        g1u_1 = re.search('(?<=\+- ).+g',out)
        g2u_1 = re.search('g2:.+r',out)
        g2u_2 = re.search('(?<=\+-).+r',g2u_1.group(0))
        g1u = float(g1u_1.group(0)[:-1])
        g2u = float(g2u_2.group(0)[:-1])
        glist.append([g1,g2])
        ulist.append([g1u,g2u])
    return glist, ulist, used

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i+n]
        
def main(args):
    glist = []
    elist = []
    used = 0
    n_trials = len(moments)
    #c = chunks(moments, 11)
#    for chunk in c:
    glist, elist, u = gmeasure(const_script, moments)
#    glist.append(np.mean(g, axis=0))
#    elist.append(np.mean(e, axis=0))
    used += u
    print 'number of trials: ', n_trials
    print 'Used ', used
    print 'Std Dev: ', np.std(glist, axis=0)[0], np.std(glist, axis=0)[1]
    print 'average g1: ', np.mean(glist, axis=0)[0],'+-', np.mean(elist, axis=0)[0]/np.sqrt(n_trials)
    print 'average g2: ', np.mean(glist, axis=0)[1],'+-', np.mean(elist, axis=0)[1]/np.sqrt(n_trials)

#    g1arr = np.array([i[0] for i in glist])
#    g2arr = np.array([j[1] for j in glist])
#    g1errs = np.array([i[0] for i in ulist])
#    g2errs = np.array([j[1] for j in ulist])
#    gplot(g1arr, g1errs, 'g1', args[1])
#    gplot(g2arr, g2errs, 'g2', args[1])

if __name__=="__main__":
    main(sys.argv)
