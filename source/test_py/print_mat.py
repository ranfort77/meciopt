#!/usr/bin/env python

import numpy as np

def main():
    v1 = np.array([1,2,3])
    v2 = np.array([4,5,6])
    v3 = np.array([7,8,9])
    mt = np.column_stack([v1,v2,v3])
    print mt.shape
    print mt.size
    print mt 
    for n in range(mt.shape[0]):
        print '%s %f %f %f' %tuple(['a'] + list(mt[n,:]))
    print mt
    mt = mt.flatten()
    print mt 
    print mt.reshape(len(mt)/3, 3) 

if __name__ == '__main__':
    main()
