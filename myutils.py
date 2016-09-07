import numpy as np
from astropy.io import fits

def match(list1, list2):
    """
    Finds items in list2 that match items in list1
        if both lists are sorted and sortable.

    Returns indices of matches in each list (index1, index2).

    Faster than using 'in' for long lists.
    For fastest performance the first argument should
        be the longer list.
    """

    index1,index2 = [],[]
    i,j = 0,0

    if not (np.all(list1 == np.sort(list1)) or np.all(list2 == np.sort(list2))):
        print "Lists not sorted. Try again with sorted lists. None returned."
        return None

    else:
        while i < len(list1):
            while j < len(list2):

                if list2[j]==list1[i]:
                    index1.append(i)
                    index2.append(j)
                
                    i+=1
                    j+=1
                    continue

                else:
                    
                    if list2[j] < list1[i]:
                        j+=1
                        continue
                
                    if list2[j] > list1[i]:
                        i+=1
                        continue
            break

        return index1, index2

def fitsstack(table_list):
    """
    Arguments:
        table_list - list of fits file names
        
    Returns:
        combined HDU using data from HDU 1
    """
    nrows=[]
    data=[]

    for table in table_list:
        t = fits.open(table)
        data.append(t[1].data)
        nrows.append(t[1].data.shape[0])
        
    hdu = fits.BinTableHDU.from_columns(t[1].columns, nrows=sum(nrows), fill=True)

    for colname in t[1].columns.names:
        for i in range(len(data)):
            hdu.data[colname][sum(nrows[:i]):sum(nrows[:i+1])] = data[i][colname]
    
    check=0
    for i in range(len(data)):
        check.append(hdu.data[colname][sum(nrows[:i])]==data[i][colname][0])
    if np.all(check):
        return hdu
    else:
        raise ValueError('Something went wrong')

