#!/usr/bin/env python
"""
Geometry optimization of the minimum energy conical intersection

Requirements:
    1. Molpro-2015
    2. input option file 

Author: Heesun, An

Date of last edit: 2018-04-16 

"""
import os
import sys
import re
import molprolib as ml
import pprint
from optparse import OptionParser
import numpy as np
import collections
#import matplotlib as mpl
#import matplotlib.pyplot as plt
from penaltyfun import Penaltyfunction
from optimize import linesearch_halfstep, linesearch_armijo

__version__ = '1.0.0'


# ---------------------------------------------------------------------
# function hierarchy
# ---------------------------------------------------------------------
# main
#  
#   - process_cmdline
#       - print_summary
#       - get_input_params
#
#   - get_nonexistence_dirname
#   - make_input_and_run
#   - check_implement_limit
#   - setup
#
#   - runopt_penaltyfunction
#       - calc - make_input_and_run
#       - compute_G_bfgs
#       - compute_G_ms
# ---------------------------------------------------------------------

# conversion factor (molpro manual 2015 chapter 8.8.1)
AU2EV = 27.2113839
AU2KCAL = 627.5096
AU2ANG = 0.529177209
AU2CM = 219474.63067

# Atoms
ATOMS = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
         'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']

# Defined variables and values 
# format : variable:(values,...)
INPUT_PARAMS = collections.OrderedDict()
INPUT_PARAMS['PRINT'] = (1,)
INPUT_PARAMS['OPTMTD'] = ('PF',)
INPUT_PARAMS['LINESEARCH'] = ('ARMIJO', 'HALFSTEP')
INPUT_PARAMS['PFUN'] = ('jpcb2008', 'jpcl2011', 'cej2004', 'gamma')
INPUT_PARAMS['PF_ALPHA'] = (3.50,)
INPUT_PARAMS['PF_BETA'] = (0.02,)
INPUT_PARAMS['PF_MAXCYC'] = (30,)
INPUT_PARAMS['THRESH_STEP'] = (1.0e-6,)
INPUT_PARAMS['THRESH_GRAD'] = (0.005,)
INPUT_PARAMS['THRESH_EGAP'] = (0.001,)
INPUT_PARAMS['NEWG'] = ('IDENTITY', 'OLDG', 'NEWCALC')
INPUT_PARAMS['UPDATEG'] = ('BFGS', 'MS')
INPUT_PARAMS['OPTCOORD'] = ('XYZ', 'Z-matrix')
INPUT_PARAMS['DIFDR'] = (0.01,)
INPUT_PARAMS['IDXLOW'] = (1,)
INPUT_PARAMS['IDXUPP'] = (2,)
INPUT_PARAMS['INPSTR'] = ('molpro single point input, Here',)

# This templet string is printed by 'meciopt -t' command
INPUT_TEMPLET = """# MECIOPT INPUT TEMPLET
# This input file is run in meciopt script via the 'exec' command, 
# so the variables defined in this input file are equivalent  
# the variables in meciopt.
# The option values are case-sensitive.
# ---------------------------------------------------------------------
# Optimization options (str): 'PF'
OPTMTD = 'PF' 
PFUN = 'jpcb2008'
LINESEARCH = 'HALFSTEP' 
PF_ALPHA = 3.50
PF_BETA = 0.02
PF_MAXCYC = 30
THREH_STEP = 1.0e-6
THRESH_GRAD = 0.005
THRESH_EGAP = 0.001
NEWG = 'OLDG'  # IDENTITY, OLDG, NEWCALC
# update method of G matrix (str): 'BFGS', 'MS'
UPDATEG = 'BFGS' 
# Optimization coordinates (str): 'XYZ', 'Z-matrix'
OPTCOORD = 'XYZ' 
# ---------------------------------------------------------------------
# molpro single point energy input (str)
INPSTR = '''***,h2s
memory,128,m

basis=6-31g*

symmetry,x
geometry={
S
Q 1 LR
H 2 HSR 1 GM
H 2 HSR 1 GMM 3 DD}

LR = 1.9989 bohr
SR = 4.2904 bohr
HSR = SR/2
GM = 90.0 degree
GMM = 180 - GM
DD = 180.0 degree

{hf;occ,7,2;wf,18,1,0}

{multi;
    closed,4,1;
      occ,10,2;
    wf,18,2,2;state,2;weight,0.55,0.45;
    wf,18,2,0;state,2;weight,0.55,0.45;
    wf,18,1,0;state,1;weight,0.55}

---
'''
# ---------------------------------------------------------------------
# IDXUPP = index of Upper-state,  IDXLOW = index of Lower-state (int)
IDXLOW = 3
IDXUPP = 4
# numerical difference delta R in bohr (float)
DIFDR = 0.01 
# print option (int)
#   1 : print until vectors (ex. F, R ...)
#   2 : print until matries (ex. G, P ...)
#   3 : print until system massage
PRINT = 1 
"""


# ---------------------------------------------------------------------
# functions to
# process command line, get input parameters and control directory 
# ---------------------------------------------------------------------
def get_nonexistence_dirname(dirname, intsuffix=1):
    """recursive function"""
    dirname_int = '%s_%d' %(dirname, intsuffix)
    if os.access(dirname_int, os.F_OK):
        dirname_int = get_nonexistence_dirname(dirname, intsuffix + 1) 
    return dirname_int

def process_cmdline():
    """ parser of command line arguments and options"""
    parser = OptionParser(usage="%prog inputfile [options]...\n\n"
             "Example: %prog -h                 # option help\n"
             "         %prog -t                 # print inputfile templet\n"
             "         %prog -t > inputfile     # save  inputfile templet\n"
             "         %prog inputfile          # run inputfile in default workdir\n"
             "         %prog inputfile -w temp  # run inputfile in workdir is temp/", 
                version='%%prog %s' %__version__)
    parser.add_option("-t", "--templet", action='store_true', default=False,
                      dest="istemplet", 
                      help="print templet of meciopt input. " 
                "If this option is activated, all other options are ignored.")
    parser.add_option('-w', '--workdir', action='store', type='string', 
                      dest='workdir', metavar='WORKDIR', 
                      help="set work dir. "
                "If the specified WORKDIR does not exist, it is created."
                "If it already exists, the existing workdir is changed to"
                "WORKDIR[_1[_2...]] and a new WORKDIR is created.")
    parser.add_option("-m", "--molprorun", action='store', type='string',
                      dest="molpro_inputfile", metavar='MOLPRO_INPUTFILE', 
                      help="run molpro calculation") 
    parser.add_option("-s", "--summary", action='store', type='string',
                      dest="smecioptout", metavar='MECIOPTOUT', 
                      help="print geometries in g09 format and summary table") 
    if len(sys.argv) == 1:
        (options, args) = parser.parse_args(['-h'])
    else:
        (options, args) = parser.parse_args()
    if options.istemplet:
        print INPUT_TEMPLET
        sys.exit()
    elif options.molpro_inputfile:
        run_molpro(options.molpro_inputfile) 
        sys.exit()
    elif options.smecioptout:
        print_summary(options.smecioptout) 
        sys.exit()
    if len(args) < 1:
        msg_exit('There is no inputfile') 
    else:
        inputfile = args[0]
    # check the existence of inputfile 
    if not os.access(inputfile, os.F_OK):
        msg = "inputfile '%s' does not exist." %inputfile
        msg_exit(msg) 
    # check workdir
    if options.workdir is None:
        # workdir is set to the inputfile name excluding the extension. 
        workdir = os.path.splitext(os.path.basename(inputfile))[0]
    else:
        workdir = options.workdir
    options = get_input_params(inputfile)
    return options, workdir

def get_input_params(inputfile):
    """get input parameters from input file"""
    fp = open(inputfile, 'r')
    inpstr = fp.read()
    fp.close()
    del inputfile
    del fp
    exec inpstr
    del inpstr
    params = locals()
    # check input params
    for var in params:
        # check whether defined option
        if var not in INPUT_PARAMS:
            msg = 'Unknown input option: "%s"' %var
            msg_exit(msg)
        # check data type of option value
        curval = params[var]
        regvals = INPUT_PARAMS[var]
        if type(curval) != type(regvals[0]):
            msg = 'value type of "%s" must be "%s"' %(var, type(regvals[0]))
            msg_exit(msg)
        # check whether defined value
        if len(regvals) > 1:
            if curval not in regvals:
                msg = 'Unknown value of "%s" option: "%s"' %(var, curval)
                msg_exit(msg)
    #endfor
    # ordering params
    oparams = collections.OrderedDict()
    for var in INPUT_PARAMS:
        if var in params:
            oparams[var] = params[var]
    return oparams

def run_molpro(inpfilename):
    inpobj = ml.Input(inpfilename)
    runmsg = inpobj.run()
    println('stdout of molpro ...\n    %s' %runmsg['stdout'])
    println('stderr of molpro ...\n    %s' %runmsg['stderr'])
    
def print_summary(mecioptout):
    # parsing data
    atom = []
    xyzs = []
    normsd = []  # norm of search direction (size = k - 1)
    normns = []  # norm of next step (size = k - 1)
    alpha = []
    beta = []
    sigma = []
    para = []; ispara = []
    perp = []; isperp = []
    step = []; isstep = []
    gamma = []; isgamma = []
    nls = []  # number of line search (size = k - 1) 
    steplength = [] #(size = k - 1)
    L = []
    P = []
    # start parsing
    fp = open(mecioptout, 'r')
    for line in fp:
        if line[3:33] == 'Summary of Initial Calculation':
            [fp.next() for n in range(6)] # skip 6 lines
            xyz = []
            for line in fp:
                foo = line.split()
                if len(foo) != 5:
                    break 
                else:
                    atom.append(foo[1]) 
                    xyz.append([float(e) for e in foo[2:5]])
            xyzs.append(xyz) 
        if line[0:25] == ' Calculates the next step':
            fp.next(); fp.next(); fp.next()
            xyz = []
            for n in range(len(atom)):
                x = float(fp.next().split()[4])
                y = float(fp.next().split()[4])
                z = float(fp.next().split()[4])
                xyz.append([x, y, z])
            xyzs.append(xyz)
            fp.next()
            normsd.append(float(fp.next().split('=')[1]))
            normns.append(float(fp.next().split('=')[1]))
        if line[6:18] == 'MECIOPT (k =': 
            alpha.append(float(fp.next().split()[3])) 
            beta.append(float(fp.next().split()[3])) 
            fp.next()
            fp.next()
            sigma.append(float(fp.next().split()[2])) 
            fp.next() 
            fp.next()
            L.append(float(fp.next().split()[3])) 
            P.append(float(fp.next().split()[3])) 
            fp.next()
            fp.next()
            foo = fp.next().split()
            para.append(float(foo[2][:-1]))
            ispara.append(foo[4])
            foo = fp.next().split()
            perp.append(float(foo[1]))
            isperp.append(foo[3])
            foo = fp.next().split()
            step.append(float(foo[2][:-1]))
            isstep.append(foo[4])
            foo = fp.next().split()
            gamma.append(float(foo[2]))
            isgamma.append(foo[4])
            ils = None
            for line in fp:
                if ',i=' in line: 
                    ils = int(line.split(',')[1].split('=')[1].split(')')[0])
                    nls.append(ils) 
                    fp.next(); fp.next()
                    steplength.append(float(fp.next().split('=')[1]))
                    break
    fp.close()
    # print geometry in g09 format
    def g09_stdori(atom, xyzbohr):
        res = []
        header = []
        header.append('                         Standard orientation:')
        header.append(' ---------------------------------------------------------------------')
        header.append(' Center     Atomic     Atomic              Coordinates (Angstroms)')
        header.append(' Number     Number      Type              X           Y           Z')
        header.append(' ---------------------------------------------------------------------')
        res.append('%s' %'\n'.join(header))
        for n in range(len(atom)):
            lb = n + 1
            atomNum = ATOMS.index(atom[n]) + 1
            atomTyp = 0
            x = xyzbohr[n][0] * AU2ANG
            y = xyzbohr[n][1] * AU2ANG
            z = xyzbohr[n][2] * AU2ANG
            res.append(' %4d  %9d  %12d     %11.6f %11.6f %11.6f' %(lb, 
                             atomNum, atomTyp, x, y, z))
        res.append(' ---------------------------------------------------------------------')
        return '\n'.join(res)
    sgrad = ' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad' 
    g09str = []
    g09str.append('     1 basis functions,     1 primitive gaussians,     1 cartesian basis functions')
    g09str.append('     1 alpha electrons        1 beta electrons')
    g09str.append(sgrad)
    g09str.append(sgrad)
    for n, xyz in enumerate(xyzs):
        k = n + 1
        g09str.append(g09_stdori(atom, xyz))
        g09str.append(sgrad)
        g09str.append(' Step number %3d' %k)
        g09str.append(sgrad)
    g09str.append('') 
    g09str.append(' Normal termination of Gaussian 09')
    g09str.append('') 
    print '\n'.join(g09str)
    print '\n'
    # print summary table
    nint = '%4d'      # '   1'
    nene = '%16.8f'   # '  -1234.12345678'
    nden = '%13.8f'   #    '  -0.12345678'
    npar = '%10.4f'   #  ' +1234.1234'
    sint = '%4s'
    sene = '%16s'
    sden = '%13s'
    spar = '%10s'
    # table 1 
    tabone = []
    tabone.append(' Summary Table 1') 
    sfmt1 = ''.join([sint, sene, sden, sene, sden, spar, spar])
    head1 = sfmt1 %('k','Sigma','Gamma', 'L','P','alpha','beta')
    tabone.append(head1)
    nfmt1 = ''.join([nint, nene, nden, nene, nden, npar, npar])
    # table 2
    tabtwo = []
    tabtwo.append(' Summary Table 2') 
    sfmt2 = ''.join([sint, sden, sden, sden, sint, 
                     sden, sden, sene, sden, 
                     spar, spar, spar, spar])
    head2 = sfmt2 %('k', '|p(k)|', '|dR(k)|', 'steplen', '#LS',
                 'para', 'perp', 'step', 'egap', 
                 '?para','?perp','?step','?egap')
    tabtwo.append(head2)
    nfmt2 = ''.join([nint, nden, nden, nden, nint, 
                     nden, nden, nene, nden, 
                     spar, spar, spar, spar])
    # record data 
    siz = len(sigma)
    for n in range(siz):
        k = n + 1
        sg = sigma[n]
        gm = gamma[n]
        el = L[n]
        pi = P[n]
        ap = alpha[n]
        bt = beta[n]
        pa = para[n]
        pe = perp[n]
        st = step[n]
        ispa = ispara[n]
        ispe = isperp[n]
        isst = isstep[n]
        isgm = isgamma[n]
        if k == siz:
            sl = 0. 
            nl = 0. 
            nsd = 0. 
            nns = 0. 
        else:
            sl = steplength[n]
            nl = nls[n]
            nsd = normsd[n]
            nns = normns[n]
        tabone.append(nfmt1 %(k, sg, gm, el, pi, ap, bt))
        tabtwo.append(nfmt2 %(k, nsd, nns, sl, nl, 
                              pa, pe, st, gm,
                              ispa, ispe, isst, isgm))
    print '\n'.join(tabone) 
    print '\n'
    print '\n'.join(tabtwo) 
    print '\n'

# ---------------------------------------------------------------------
# functions to print output
# ---------------------------------------------------------------------

def msg_exit(msg):
    print msg 
    print '**** exit meciopt ****'
    sys.exit()

def println(s):
    print ' %s' %s

def print_section_title(title, symbol, maxcol=None):
    siz = len(title)
    if maxcol is None:
        maxcol = siz + 5
    print 
    print ' %s' %((maxcol-1)*symbol)
    print ' %s %s %s' %(symbol, title.center(maxcol-5), symbol)
    print ' %s' %((maxcol-1)*symbol)
    print 

def print_underline_title(title):
    println(title)
    println('%s\n' %(len(title)*'-'))

def print_params(params):
    for var in params:
        if 'INPSTR' in var:
            print ' %20s =\n%s' %(var, params[var]) 
        else:
            print ' %20s = %-51s' %(var, params[var]) 

def print_geometry(atom, xyz):
    print ' %8s%8s%14s%14s%14s' %('LABEL', 'ATOM', 'X', 'Y', 'Z') 
    for n in range(len(atom)):
        print ' %8d%8s%14.8f%14.8f%14.8f' %(n+1, atom[n], 
                    xyz[n][0], xyz[n][1], xyz[n][2])
    print 

def print_energy(molprolib_final_energy, idxupp, idxlow):
    order = ['SCF', 'MULTI', 'CI']
    finprog = None
    for substr in order:
        for program in molprolib_final_energy:
            if substr in program:
                finprog = program
                print ' %10s%10s%10s%10s%18s' %('PROGRAM', 'INDEX', 'STATE',
                                            'SPIN', 'ENERGY(a.u.)')
                d = molprolib_final_energy[program]
                for n in range(len(d['state'])):
                    state = d['state'][n]
                    try:
                        spin = d['spin'][n]
                    except IndexError:  # spin 
                        # For HF program, spin is output only when 
                        # explicitly specfied in input.
                        # 'spin' key exists but element is empty list
                        spin = ''
                    energy = d['energy'][n]
                    print ' %10s%10d%10s%10s%18.8f' %(program, n + 1, 
                                                  state, spin, energy)
                #endfor
                print 
        #endfor
    #endfor
    d = molprolib_final_energy[finprog]
    up = d['energy'][idxupp-1] # python is zero-base
    lo = d['energy'][idxlow-1]
    sg = (up + lo) / 2.0
    gm = (up - lo)
    print ' %10s of Upper = %d, Lower = %d' %('INDEX', idxupp, idxlow)
    print ' %10s%18s%18s%18s' %('STATE', 'ENERGY(a.u.)', '(eV)', '(kcal/mol)') 
    print ' %10s%18.8f' %('Upper', up)
    print ' %10s%18.8f' %('Lower', lo)
    print ' %10s%18.8f' %('SIGMA', sg)
    print ' %10s%18.8f%18.8f%18.8f' %('GAMMA', gm, gm*AU2EV, gm*AU2KCAL)
    print 

def print_sp_summary(params, outobj):
    # get xyz coordinates from output 
    pg = outobj.get_final_pointgroup() 
    symel = outobj.get_final_symmetry_elements() 
    atom, xyz = outobj.get_xyzcoords('final')
    println('Point group      : %s' %pg)
    println('Symmetry Elements: %s\n' %','.join(symel))
    # print atomic coordinates and energies
    print_geometry(atom, xyz)
    print_energy(outobj.get_final_energies(), 
                 params['IDXUPP'], params['IDXLOW'])

def print_vectors(array2d, colnames=None, rownames='3N'):
    """
    >>> mt = np.array([[1,2,3,4],[5,6,7,8],[9,10,11,12]])
    >>> print_vectors(mt, ['abc','def','ghi','jkl'])
    """
    nrow = array2d.shape[0]
    ncol = array2d.shape[1]
    if rownames == '3N':
        natom = nrow / 3
        rownames = []
        for n in range(natom):
            for axis in ['x', 'y', 'z']:
                rownames.append('%s%d' %(axis, n+1))
    if colnames is None:
        colnames = ['c%d' %(n+1) for n in range(ncol)]
    colnames.insert(0, '')
    print ('%5s' + ncol*'%18s') %tuple(colnames)
    for n in range(nrow):
        print ('%5s' + ncol*'%18.8f') %tuple([rownames[n]] 
                                         + list(array2d[n,:]))
    print  

def print_3n3n_matrix(mat):
    """
    >>> natom = 2 
    >>> mat = np.arange(3*natom*3*natom).reshape(3*natom, 3*natom)
    >>> print mat
    >>> print_3n3n_matrix(mat)
    """
    natom = mat.shape[0] / 3 
    rownames = []
    for n in range(natom):
        for axis in ['x', 'y', 'z']:
            rownames.append('%s%d' %(axis, n+1))
    for n in range(natom):
        colnames = ['%s%d' %(axis, n+1) for axis in ['x', 'y', 'z']]
        colnames.insert(0, '')
        print '%5s%18s%18s%18s' %tuple(colnames) 
        for nrow in range(natom*3):
            print ('%5s' + 3*'%18.8f') %tuple([rownames[nrow]] 
                                + list(mat[nrow,(3*n+0):(3*n+3)]))
        print

# ---------------------------------------------------------------------
# interface for ab-initio programs 
# ---------------------------------------------------------------------
def make_input_and_run(inpname, inpstr, copywfuname=None, iprint=1):
    if copywfuname is None:
        igobj = ml.InputGenerator(inpname, inpstr)
    else:
        igobj = ml.InputGenerator(inpname, inpstr, copywfuname)
    igobj.make_input()
    inpobj = ml.Input(inpname)
    runmsg = inpobj.run()
    outobj = ml.Output(inpobj.out)
    normalterm = outobj.check_normalterm()
    if not normalterm:
        return normalterm, inpobj, outobj
    println('make & run input ...\n    %s' %inpname)
    if iprint > 2:
        println('%s\n' %repr(inpobj))
        println('stdout of molpro ...\n    %s' %runmsg['stdout'])
        println('stderr of molpro ...\n    %s' %runmsg['stderr'])
    println('total wall time: %s sec\n' %outobj.get_realtime())
    return normalterm, inpobj, outobj

def check_implement_limit(params, outobj):
    pg = outobj.get_final_pointgroup() 
    symel = outobj.get_final_symmetry_elements() 
    if params['OPTCOORD'] == 'XYZ':
        if pg not in ['C1', 'Cs']:
            msg = 'Only C1 or Cs can be calculated'
            msg = msg + '\ncurrent point group = %s' %pg
            msg_exit(msg) 
        elif pg == 'Cs': 
            if symel[0] != 'X':
                msg = 'If Cs point group, symmetry element must be X' 
                msg = msg + '\ncurrent symmetry element = %s' %','.join(symel)
                msg_exit(msg)
    elif params['OPTCOORD'] == 'Z-matrix':
        msg_exit('not implement, yet')

# ---------------------------------------------------------------------
# main 
# ---------------------------------------------------------------------

def main():
    
    # get input parameters and workdir from command line
    params, workdir = process_cmdline()
    print_section_title('START MECIOPT', '*', 72)
    print_section_title('Input Parameters', '-')
    print_params(params)

    # make work directory.
    # If the workdir already exists, rename old-workdir to workdir_[int]
    if os.access(workdir, os.F_OK):
        workdir_int = get_nonexistence_dirname(workdir)
        os.rename(workdir, workdir_int)
    print_section_title('Working Directory', '-')
    os.mkdir(workdir)
    println('make work directory ...\n    mkdir %s' %workdir)
    os.chdir(workdir)
    println('move to work directory ...\n    cd %s' %workdir)

    # make initial input and run
    print_section_title('Run Initial Calculation', '-')
    normalterm, inpobj, outobj = make_input_and_run('step_001_a.inp', 
                                 params['INPSTR'], None, params['PRINT'])
    if not normalterm:
        msg = 'The qc program not conversed.'
        msg_exit(msg)    

    # print summary of single point calculation
    print_section_title('Summary of Initial Calculation', '-')
    print_sp_summary(params, outobj)

    # limitation of current implementation
    check_implement_limit(params, outobj)

    # setup for meciopt
    optparams = setup(params, inpobj, outobj)
    print_section_title('Setup Parameters for Meciopt', '-')
    print_params(optparams)

    # start meciopt 
    if params['OPTMTD'] == 'PF': 
        print_section_title('Start Meciopt with Penalty-Function', '-')
        runopt_penaltyfunction(optparams)

    print_section_title('END MECIOPT', '*', 72)


# ---------------------------------------------------------------------
# functions for MECIOPT 
# ---------------------------------------------------------------------
def setup(params, inpobj, outobj):
    pfparams = collections.OrderedDict() 
    pfparams['PFUN'] = params['PFUN']
    pfparams['LINESEARCH'] = params['LINESEARCH']
    pfparams['PRINT'] = params['PRINT']
    pfparams['UPDATEG'] = params['UPDATEG']
    pfparams['PF_ALPHA'] = params['PF_ALPHA']
    pfparams['PF_BETA'] = params['PF_BETA']
    pfparams['PF_MAXCYC'] = params['PF_MAXCYC']
    pfparams['THRESH_STEP'] = params['THRESH_STEP']
    pfparams['THRESH_GRAD'] = params['THRESH_GRAD']
    pfparams['THRESH_EGAP'] = params['THRESH_EGAP']
    pfparams['NEWG'] = params['NEWG']
    pfparams['DIFDR'] = params['DIFDR']
    pfparams['IDXLOW'] = params['IDXLOW']
    pfparams['IDXUPP'] = params['IDXUPP']
    atom, xyz = outobj.get_xyzcoords('final')
    pfparams['ATOM'] = atom 
    pfparams['NATOM'] = len(xyz)
    pfparams['INIXYZ'] = xyz
    symel = outobj.get_final_symmetry_elements() 
    pfparams['SYMEL'] = symel 
    blocks = inpobj.get_blocks_inpstr()
    # get input string of energy
    inpstr_multi = '' 
    inpstr_ci = '' 
    for block in blocks:
        if block['keyword'] in 'multi':
            inpstr_multi = block['inpstr'] 
        if block['keyword'] in 'ci':
            inpstr_ci = block['inpstr'] 
    # get nstates 
    if inpstr_ci == '':
        nstates = re.findall('state,\w*(\d)+', inpstr_multi)
    else:
        nstates = re.findall('state,\w*(\d)+', inpstr_ci)
    pfparams['NSTATES'] = sum([int(e) for e in nstates])
    pfparams['INPSTR_MULTI'] = inpstr_multi  
    pfparams['INPSTR_CI'] = inpstr_ci 
    if nstates < 2:
        msg = 'The last energy input must be included more than one state'
        msg_exit(msg)
    return pfparams


def runopt_penaltyfunction(optparams):

    pf = Penaltyfunction(optparams['PFUN'])
    println('----------------------------')
    println(' Penalty function %s' %optparams['PFUN'])
    println('----------------------------')

    thresh_step = optparams['THRESH_STEP'] 
    thresh_grad = optparams['THRESH_GRAD'] 
    thresh_egap = optparams['THRESH_EGAP'] 

    kmax = optparams['PF_MAXCYC'] 
    println('----------------------------')
    println(' thresholds for convergence ')
    println('----------------------------')
    println(' %8s = %14.8f' %('step', thresh_step))
    println(' %8s = %14.8f' %('grad', thresh_grad))
    println(' %8s = %14.8f' %('EGap', thresh_egap))
    println(' %8s = %4d\n' %('kmax', kmax))

    pf.alpha = optparams['PF_ALPHA'] 
    pf.beta  = optparams['PF_BETA'] 

    R = np.array(optparams['INIXYZ']).flatten()
    k = 1
    Lprev = 0.0
    I = np.eye(3*optparams['NATOM'])
    kscr = Kscratch(kmax)
    # main loop
    while True:

        # stop at maximum iteration
        if k > kmax:
            println('!BREAK: kmax = %d was reached' %kmax)
            break

        # calculate numerical gradient 
        try:
            if k == 1:
                copywfuname = 'step_001_a.wfu'
            else:
                copywfuname = 'step_%03d_nd.wfu'  %(k-1)
            inpname = 'step_%03d_nd.inp' %k
            sg, gm, sgp, gmp = calc('nd', optparams, R, inpname, copywfuname)
        except TypeError:
            println('!BREAK: qc program not converged.')
            break

        # compute data of penalty function  
        pf.sigma = sg
        pf.gamma = gm
        pf.sigmaprime = sgp
        pf.gammaprime = gmp 
        pf.update() 

        # check convergence
        step = pf.L - Lprev
        pass_step = abs(step) <= thresh_step
        pass_gradpara = abs(pf.gradpara) <= thresh_grad 
        pass_gradperp = pf.gradperp <= thresh_grad 
        pass_egap = pf.gamma <= thresh_egap
                
        # print convergence info.
        lw = sg - gm/2.0
        up = sg + gm/2.0
        fmt = '%20s = %18.8f'
        println('---- MECIOPT (k = %2d) ----' %k)
        println(fmt %('PF, alpha', pf.alpha))
        println(fmt %('PF, beta', pf.beta))
        println((fmt + ' a.u.') %('Lower state', lw))
        println((fmt + ' a.u.') %('Upper state', up))
        println((fmt + ' a.u.') %('Sigma', sg))
        println((fmt + ' a.u. (%15.5f eV)\n') %('Energy gap', gm, gm*AU2EV))
        println((fmt + ' a.u.') %('objective L', pf.L))
        println((fmt + ' a.u.\n') %('penalty P', pf.P))
        println('%20s%18s%14s%8s' %('check convergence', 'value', 
                                    'thresh', 'convsed'))
        println('%20s%2s%16.8f%-2s%12.8f%8s' %('parallel', 
                                        '|', pf.gradpara, '|', 
                                        thresh_grad, pass_gradpara))
        println('%20s%18.8f%14.8f%8s' %('perpendicular', pf.gradperp, 
                                        thresh_grad, pass_gradperp))
        println('%20s%2s%16.8f%-2s%12.8f%8s' %('step', '|', step, '|', 
                                        thresh_step, pass_step))
        println('%20s%18.8f%14.8f%8s\n' %('Energy gap', pf.gamma, 
                                        thresh_egap, pass_egap))

        # braching 
        if pass_gradpara and pass_gradperp and pass_step: 
            if pass_egap:
                println('!BREAK: convergence reached!')
                break
            else:
                println('!PASS: step and grad!')
                println('!NOT PASS: egap\n')
                # update alpha
                pf.alpha = 2.0*pf.alpha
                println('!UPDATE: alpha = 2*alpha\n')
                # compute G(k-1) for three categories
                if optparams['NEWG'] == 'IDENTITY':
                    println('!SET: G(k-1) = Identity matrix\n')
                    G = I
                elif optparams['NEWG'] == 'NEWCALC':
                    println('!SET: G(k-1) = NEWCALC \n')
                    for ikm1 in range(k-1): # ikm1 is 0, 1, ..., k-2 
                        pf.sigma = kscr.sigma[ikm1]
                        pf.gamma = kscr.gamma[ikm1]
                        pf.sigmaprime = kscr.sigmaprime[ikm1]
                        pf.gammaprime = kscr.gammaprime[ikm1]
                        pf.update()
                        iF = pf.Lprime 
                        if ikm1 == 0:
                            iG = I
                        else:
                            iG = compute_G_bfgs(idRp, iF-iFp, iG) 
                        iFp = iF
                        idRp = kscr.dR[ikm1]
                    G = iG 
                elif optparams['NEWG'] == 'OLDG':
                    # variable G is already G(k-1) of old alpha 
                    println('!SET: G(k-1) = OLDG \n')
                # compute "Fprev" for new alpha
                #     notice: k start 1 but list index start 0.
                #             Therefore, index of (k-1)th geometry is k-2. 
                println('!SET: F(k-1) for new alpha \n')
                ix = k - 2 
                pf.sigma = kscr.sigma[ix] 
                pf.gamma = kscr.gamma[ix] 
                pf.sigmaprime = kscr.sigmaprime[ix] 
                pf.gammaprime = kscr.gammaprime[ix] 
                pf.update()
                Fp = pf.Lprime
                # re-set to compute "F" for new alpha 
                pf.sigma = sg 
                pf.gamma = gm 
                pf.sigmaprime = sgp 
                pf.gammaprime = gmp 
                pf.update()
        else:
            println('!NOT PASS: step and grad!\n')

         
        F = pf.Lprime # F(k)
        # update inverse Hessian G 
        if k == 1:
            println('!Initial G is Identity matrix')
            G = I
        else:     
            println('!UPDATE: G using %s\n' %optparams['UPDATEG'])
            G = compute_G_bfgs(dR, F-Fp, G, optparams['PRINT'])

        pk = -np.matmul(G, F)  # search direction

        # line search
        print_underline_title('Line search')
        Lk = pf.L
        iname = 'linesearch.inp'
        cwfuname = 'step_%03d_nd.wfu'  %k
        def phi(tsl): # tsl = trial step length
            tsg, tgm = calc('sp', optparams, R + tsl*pk, iname, cwfuname)
            pf.sigma = tsg
            pf.gamma = tgm
            pf.update()
            return pf.L
        if optparams['LINESEARCH'] == 'ARMIJO':
            println('linesearch by armijo\n')
            steplength, Lki, i = linesearch_armijo(phi, Lk, np.dot(F, pk))
        elif optparams['LINESEARCH'] == 'HALFSTEP':
            println('linesearch by half step\n')
            steplength, Lki, i = linesearch_halfstep(phi, Lk)
        println('%20s = %18.8f' %('L(k=%d)' %k, Lk))
        println('%20s = %18.8f' %('L(k=%d,i=%d)' %(k,i), Lki))
        println('%20s = %18.8f' %('DIFF', Lki-Lk))
        if steplength is not None:
            println('!FIND step length')
            println('!step length = %14.8f\n' %steplength) 
        else:
            msg_exit('line-search failed') 

        print_underline_title('Calculates the next step.')
        dR = steplength*pk
        print_vectors(np.column_stack([pk, R, dR, R+dR]), 
                   ['p(k)', 'R(k)', 'dR(k)', 'R(k+1)'])
        println('%40s = %18.8f' %('norm of search direction, |p(k)|', 
                                   np.linalg.norm(pk)))
        println('%40s = %18.8f' %('norm of next step,       |dR(k)|', 
                                   np.linalg.norm(dR)))
        println('')
        
        # save data to scratch for BFGS-G of new alpha and beta
        ix = k-1
        kscr.sigma[ix] = sg 
        kscr.gamma[ix] = gm
        kscr.sigmaprime[ix] = sgp 
        kscr.gammaprime[ix] = gmp
        kscr.dR[ix] = dR

        # update for next iteration
        Lprev = Lk
        R = R + dR
        Fp = F 
        k += 1
    #end-while

    print_section_title('End of optimization with Penalty-Function method', '-')
#end-funcion

def compute_G_bfgs(dR, D, Gold, iprint=1):
    # out = outer(a, b) -> out[i,j] = a[i]*b[j]
    QD = np.outer(dR, D)   # 3n x 3n 
    QQ = np.outer(dR, dR)  # 3n x 3n
    nom = np.matmul(dR, D) # scalar
    P = np.eye(dR.size) - (QD / nom)  # 3n x 3n 
    G = P.dot(Gold).dot(P.T) + QQ/nom # 3n x 3n 
    if iprint > 1:
        print ' BFGS, QD matrix\n' 
        print_3n3n_matrix(QD) 
        print ' BFGS, QQ matrix\n' 
        print_3n3n_matrix(QQ) 
        print ' BFGS,  P matrix\n' 
        print_3n3n_matrix(P) 
        print ' BFGS,  G matrix\n' 
        print_3n3n_matrix(G) 
    return G

def compute_G_ms(dR, D, Gold, iprint=1):
    U = dR - np.matmul(Gold, D)  # 3n vector 
    nom = np.matmul(U, D)        # scalar
    P = np.outer(U, U) / nom     # 3n x 3n
    G = Gold + P
    if iprint > 1:
        print '   MS,  P matrix\n' 
        print_3n3n_matrix(P) 
        print '   MS,  G matrix\n' 
        print_3n3n_matrix(G) 
    return G

def calc(mode, optparams, R, inpname, copywfuname):
    # unchanged params
    iprint = optparams['PRINT']
    idxlow = optparams['IDXLOW']
    idxupp = optparams['IDXUPP']
    atom = optparams['ATOM']
    symel = optparams['SYMEL']
    ns = optparams['NSTATES']
    difdr = optparams['DIFDR']
    multi = optparams['INPSTR_MULTI'].replace('}', 
                      '\nstart,2140.2;orbital,2140.2}')
    ci = optparams['INPSTR_CI']
    # changed param 
    xyz = R.reshape(optparams['NATOM'], 3)
    inpstr_e = '\n'.join([multi, ci])
    
    # get input string
    if mode == 'sp':
        inpstr = ml.InputGenerator.get_xyz_single(xyz, atom, symel, 
                                                  inpstr_e, True) 
    elif mode == 'nd':
        println('calculate gradients using central difference ...\n')
        inpstr_sp = ml.InputGenerator.get_xyz_single(xyz, atom, symel, 
                                                     inpstr_e, False) 
        multi = multi.replace('orbital,2140.2', 'orbital,2141.2')
        inpstr_e = '\n'.join([multi, ci])
        inpstr_do = ml.InputGenerator.get_xyz_numdif(xyz, atom, symel, 
                                                     inpstr_e, ns, difdr) 
        inpstr = '\n'.join([inpstr_sp, inpstr_do]) 
    # make input and run
    normalterm, inpobj, outobj = make_input_and_run(inpname, inpstr, copywfuname, iprint) 
    if not normalterm:
        raise TypeError

    # return output data
    if mode == 'sp':
        energy = outobj.get_formtable(1, 20, 
                   'single point energy', 2, 'f')
        lw = energy[0][idxlow-1]
        up = energy[0][idxupp-1]
        sg0 = (up + lw) / 2.0 
        gm0 = (up - lw) 
        print_sp_summary(optparams, outobj)
        return sg0, gm0
    elif mode == 'nd':
        # parsing data
        tablow = outobj.get_formtable(1, 15, 
                   'state %d energy' %idxlow, 2, 9*'f')
        tabupp = outobj.get_formtable(1, 15, 
                   'state %d energy' %idxupp, 2, 9*'f')
        # matsize: natom x 9 [x-,x0,x+,y-,y0,y+,z-,z0,z+] of each atom 
        # idx of zero-base   [ 0, 1, 2, 3, 4, 5, 6, 7, 8]
        lw = np.array(tablow[0])
        up = np.array(tabupp[0])
        # change to vector format: [x1,y1,z1,...,xn,yn,zn]
        LWM = lw[:,(0,3,6)].flatten()
        LWC = lw[:,(1,4,7)].flatten()
        LWP = lw[:,(2,5,8)].flatten()
        UPM = up[:,(0,3,6)].flatten()
        UPC = up[:,(1,4,7)].flatten()
        UPP = up[:,(2,5,8)].flatten()
        SGM = (UPM + LWM) / 2.0
        SGC = (UPC + LWC) / 2.0
        SGP = (UPP + LWP) / 2.0
        GMM = (UPM - LWM)
        GMC = (UPC - LWC)
        GMP = (UPP - LWP)
        # compute the values of SG and GM at k-th point
        lw0 = LWC[0] # constant
        up0 = UPC[0] # constant
        sg0 = SGC[0] # constant  
        gm0 = GMC[0] # constant 
        # compute the first-order deriv., SGIMA and GAMMA
        SG1 = (SGP - SGM) / (2.0 * difdr) # vector
        GM1 = (GMP - GMM) / (2.0 * difdr) # vector 
        print '%5s%18s%18s%18s%18s' %('', 'UPPER', 'LOWER', 'SIGMA', 'GAMMA')
        print '%5s%18.8f%18.8f%18.8f%18.8f' %('', up0, lw0, sg0, gm0)
        print 
        print_vectors(np.column_stack([UPP, LWP, UPM, LWM]), 
                       ['UPPER+', 'LOWER+', 'UPPER-', 'LOWER-'])
        print_vectors(np.column_stack([SGP, SGM, SG1]), 
                       ['SIGMA+', 'SIGMA-', 'DIFF/(2*dr)'])
        print_vectors(np.column_stack([GMP, GMM, GM1]), 
                       ['GAMMA+', 'GAMMA-', 'DIFF/(2*dr)'])
        return sg0, gm0, SG1, GM1


class Kscratch(object):
    def __init__(self, kmax):
        self.sigma = [None for n in range(kmax)]
        self.gamma = [None for n in range(kmax)]
        self.sigmaprime = [None for n in range(kmax)]
        self.gammaprime = [None for n in range(kmax)]
        self.dR = [None for n in range(kmax)] 
    


if __name__ == '__main__':
    main()




