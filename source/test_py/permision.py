#!/usr/bin/env python 

import os

tdir = '/work/tempo'
#tdir = '/work/tmp/ahs' 
#tdir = '/work/tmp/qcp' 

print os.access(tdir, os.F_OK)
print os.access(tdir, os.W_OK)


