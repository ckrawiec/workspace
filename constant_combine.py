import numpy as np
import re
import sys
import glob

#usage: python constant_combine.py <'text files to glob'>

def gaverage(moments):
    use = 0
    g1 = []
    g2 = []
    g1u = []
    g2u = []
    for moment in moments:
        f = open(moment, 'r')
        for line in f.readlines():
            if line[0:4]=='Use ':
                use+=int(line[4:10])
            if line[0:2]=='g1':
                g1re = re.search('(?<=g1: )\S+',line)                
                g1.append(float(g1re.group(0)))
                g1u_1 = re.search('(?<=\+- ).+g',line)
                g1u.append(float(g1u_1.group(0)[:-1]))
 #           if line[0:2]=='g2':
                g2re = re.search('(?<=g2: ).+r',line)                
                g2re_1 = re.search('.\S+',g2re.group(0))
                g2.append(float(g2re_1.group(0)))
                g2u_1 = re.search('(?<=\+- )\S+',g2re.group(0))
                g2u.append(float(g2u_1.group(0)[:-1]))
        f.close()
    print 'number of trials: ', len(g1)
    print 'Total used: ', use
    print 'Std. Dev: ', np.std(g1)
    print 'average g1: ', np.mean(g1),'+-', np.mean(g1u)/np.sqrt(len(g1))
    print 'average g2: ', np.mean(g2),'+-', np.mean(g2u)/np.sqrt(len(g2))
        
def main(args):
    moments = glob.glob(args[1])
    gaverage(moments)

if __name__=="__main__":
    main(sys.argv)
