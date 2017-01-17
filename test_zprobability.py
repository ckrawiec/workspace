import time
import itertools
import numpy as np
import zprobability as zp
import matplotlib.pyplot as plt
from scipy import stats
#do test with different number of chunks?

#ptree
# try different values of knear vs full integration
#make tables of test data (short)

#break up ngmix fluxes, run on some

#test data table positions

data = np.random.rand(1000,5)
errs = data*0.1
truth = data * 1.3

def run_tests():
    p_func_test()
    gauss_test()

def p_func_test():
    print '\n Comparing full integration (p) with tree results (ptree)...'
    start_ptree = time.time()
    result_ptree = zp.ptree(data, errs, truth, knear=1000)
    end_ptree = time.time()

    start_p = time.time()
    result_p = zp.p(data, errs, truth)
    end_p = time.time()
    
    notnan_p = ~np.isnan(result_p)
    notnan_ptree = ~np.isnan(result_ptree)

    print "Length of input data: {}, length of results (p and ptree): {} and {}".format(len(data), len(result_p), len(result_ptree))
    print "NaNs in ptree and p results: ", len(result_ptree[~notnan_ptree]), len(result_p[~notnan_p])
    print "Sum of p results with NaNs: {}, without NaNs: {}".format(np.sum(result_p), np.sum(result_p[notnan_p]))
    print "Time for ptree: {}, Time for p: {}".format(end_ptree-start_ptree, end_p-start_p)
    print "Means, mins, maxs: ptree ({},{},{}) p ({},{},{})".format(np.mean(result_ptree), np.min(result_ptree), np.max(result_ptree),
                                                                    np.mean(result_p), np.min(result_p), np.max(result_p))
    print "Maximum (absolute) difference between ptree and p results: ", np.max(np.abs(result_ptree[notnan_ptree]-result_p[notnan_p]))
    

def gauss_test():
    print "\n Comparing results from scipy.stats.multivariate_normal.pdf with p function..."

    scipy_result = []

    start_scipy = time.time()
    for datum, err in itertools.izip(data, errs):
        scipy_result.append(np.sum(stats.multivariate_normal.pdf(truth, mean=datum, cov=err**2.)))
    end_scipy = time.time()

    start_p = time.time()
    p_result = zp.p(data, errs, truth)
    end_p = time.time()
    
    print "Maximum (absolute) difference between scipy.stats result and p result: ", np.max(np.abs(scipy_result-p_result))
    print "Time for scipy: {}, time for p: {}".format(end_scipy-start_scipy, end_p-start_p)

run_tests()
