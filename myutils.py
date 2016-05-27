import numpy as np

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
