import numpy as np
from astropy.io import fits
from astropy.table import Column, Table, join

def joincorrecttype(tab1_file, tab2_file, col1_name, col2_name, col_type):
    """
    Fix tables where matching columns are in different formats
        by converting col2 to col_type.

    Returns astropy Table joined on that column.
    """
    tab1 = Table.read(tab1_file)
    tab2 = Table.read(tab2_file)
    
    keep = []
    for row in range(len(tab2[col2_name])):
        try:
            col_type(tab2[col2_name][row])
            keep.append(row)
        except ValueError:
            continue

    tab2 = tab2[keep]
    col2 = tab2[col2_name]
    tab2.remove_column(col2_name)
    tab2.add_column(Column(data=col2, name=col2_name, dtype=col_type))

    if (col2_name != col1_name):
        tab2.rename_column(col2_name, col1_name)
    joined = join(tab1, tab2, keys=col1_name)
    return joined

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
        while ((j < len(list2)) & (i < len(list1))):

            if list2[j]==list1[i]:
                index1.append(i)
                index2.append(j)
                
                i+=1
                j+=1
                continue
            
            if list2[j] < list1[i]:
                j+=1
                continue
                
            if list2[j] > list1[i]:
                i+=1
                continue

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
    
    check=[]
    for i in range(len(data)):
        new = hdu.data[colname][sum(nrows[:i])]
        orig = data[i][colname][0]
        if (np.isnan(new) or np.isnan(orig)):
            check.append(True)
        else:
            check.append(hdu.data[colname][sum(nrows[:i])]==data[i][colname][0])
    if np.all(check):
        return hdu
    else:
        raise ValueError('Something went wrong')

def writefitsstack(infiles, outfile):
    hdu = fitsstack(infiles)
    pri_hdu = fits.PrimaryHDU()
    hdu_list = fits.HDUList([pri_hdu, hdu])
    hdu_list.writeto(outfile)
