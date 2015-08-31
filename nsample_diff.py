import numpy as np
import glob
import sys
import re
import os

"""
usage: python nsample_diff.py <start id> <one past end id>
"""

ids = range(int(sys.argv[1]),int(sys.argv[2]))

diff1,diff2 = [],[]

for id in ids:
    before_file = os.path.join('/home','ckrawiec','moments','copies',
                               'constant_after_deselect_{}.txt'.format(id))
    after_file  = os.path.join('/home','ckrawiec','moments','copies',
                               'constant_after_deselect_nsample_50000_{}.txt'.format(id))
    fb = open(before_file)
    fa = open(after_file)
    for line in fb.readlines():
        if line[0:2]=='g1':
            g1reb = re.search('(?<=g1: )\S+',line)                
            g1b = float(g1reb.group(0))
            g2reb = re.search('(?<=g2: ).+r',line)                
            g2re_1b = re.search('.\S+',g2reb.group(0))
            g2b = float(g2re_1b.group(0))
    fb.close()
    for line in fa.readlines():
        if line[0:2]=='g1':
            g1rea = re.search('(?<=g1: )\S+',line)                
            g1a = float(g1rea.group(0))
            g2rea = re.search('(?<=g2: ).+r',line)                
            g2re_1a = re.search('.\S+',g2rea.group(0))
            g2a = float(g2re_1a.group(0))
    fa.close()
    diff1.append(g1a-g1b)
    diff2.append(g2a-g2b)

print "---Differences are (nsample 50k - nsample 20k)---"
print "Number of files: ", len(ids)
print "Average difference in g1: ", np.mean(diff1)
print "Standard deviation: ", np.std(diff1)
print "Average difference in g2: ", np.mean(diff2)
print "Standard deviation: ", np.std(diff2)
