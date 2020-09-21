#!/usr/bin/env python

import numpy as np


class PFException(Exception):
    pass 


class Desc(object):
    def __init__(self, name):
        # The name and the instance-name must be different. 
        # i.e. good: alpha = Desc('_alpha')
        #      good: alpha = Desc('ap')
        #     wrong: alpha = Desc('alpha')
        self.name = name

    def __get__(self, ins, cls):
        if not hasattr(ins, self.name):
            raise PFException('no value of attribute')
        return getattr(ins, self.name)


class Scalar(Desc):
    def __set__(self, ins, value):
        if not isinstance(value, float):
            raise PFException('value must be %s' %float) 
        setattr(ins, self.name, value)


class Vector(Desc):
    def __set__(self, ins, value):
        if not isinstance(value, np.ndarray):
            raise PFException('value must be %s' %np.ndarray) 
        setattr(ins, self.name, value)


class Penaltyfunction(object):
    defkinds = ('jpcb2008', 'jpcl2011', 'cej2004', 'gamma')
    __slots__ = ('_kind', '_ap', '_bt', '_sg', '_gm', '_sgp', '_gmp',
                 '_P', '_L', '_Pp', '_Lp',
                 '_gradpara', '_gradperp', '_U')
    # input properties
    alpha = Scalar('_ap')
    beta  = Scalar('_bt')
    sigma = Scalar('_sg')
    gamma = Scalar('_gm')
    sigmaprime = Vector('_sgp') 
    gammaprime = Vector('_gmp') 
    # output properties: set or recalc by update method 
    P = Scalar('_P')
    L = Scalar('_L')
    Pprime = Vector('_Pp')
    Lprime = Vector('_Lp')
    gradpara = Scalar('_gradpara')
    gradperp = Scalar('_gradperp')
    U = Vector('_U')

    def __init__(self, kind='jpcb2008', alpha=None, beta=None):
        # We defined that alpha is coef. in front of the penalty function,
        #             and  beta is value in the penalty function. 
        # If jpcl2011, initial alpha is 1/(0.04) and beta is None. 
        self.kind = kind 
        if alpha is not None:
            self.alpha = alpha
        else:
            if kind == 'jpcb2008':
                self.alpha = 3.5
            elif kind == 'jpcl2011':
                self.alpha = 1./0.04 
            elif kind == 'cej2004':
                self.alpha = 0.008
            elif kind == 'gamma':
                self.alpha = 3.5 
        if beta is not None:
            self.beta = beta
        else: 
            if kind == 'jpcb2008':
                self.beta = 0.02
            elif kind == 'jpcl2011':
                self.beta = 0.0  # dummy number 
            elif kind == 'cej2004':
                self.beta = 0.008
            elif kind == 'gamma':
                self.beta = 0.02

    def update(self): 
        ap = self.alpha
        bt = self.beta
        sg = self.sigma
        gm = self.gamma
        sgp = self.sigmaprime
        gmp = self.gammaprime
        if self.kind == 'jpcb2008':
            P = (gm**2)/(gm + bt) 
            Pp = (gmp*(gm**2 + 2.0*bt*gm)) / ((gm + bt)**2)
        elif self.kind == 'jpcl2011':
            P = gm**2
            Pp = 2.0*gm*gmp
        elif self.kind == 'cej2004':
            P = (bt**2)*np.log(1.0 + (gm/bt)**2)
            Pp = (2.0*gm*gmp) / (1.0 + (gm/bt)**2) 
        elif self.kind == 'gamma':
            P = gm
            Pp = gmp
        L = sg + ap*P 
        Lp = sgp + ap*Pp
        U = Pp / np.linalg.norm(Pp)
        self.P = P
        self.Pprime = Pp
        self.L = L
        self.Lprime = Lp
        self.U = U 
        self.gradpara = Lp.dot(U) / ap
        self.gradperp = np.linalg.norm(Lp - Lp.dot(U)*U)
        #self.gradperp = np.linalg.norm(Lp - (Lp.dot(U)/ap)*U)
        
    @property
    def kind(self): 
        if not hasattr(self, '_kind'):
            raise PFException('no value of attribute')
        return getattr(self, '_kind')
    @kind.setter
    def kind(self, kind):
        if kind not in self.defkinds:
            msg = '"%s" not defined.\n' %kind
            msg = msg + 'The defined kinds are %s' %', '.join(self.defkinds)
            raise PFException(msg) 
        setattr(self, '_kind', kind)

    @property
    def info(self):
        s = []
        s.append(' %20s = %14s' %('kind', self.kind))
        s.append(' %20s = %14.8f' %('alpha', self.alpha))
        s.append(' %20s = %14.8f' %('beta', self.beta))
        s.append(' %20s = %14.8f' %('sigma', self.sigma))
        s.append(' %20s = %14.8f' %('gamma', self.gamma))
        s.append(' %20s = %14.8f' %('L', self.L))
        s.append(' %20s = %14.8f' %('P', self.P))
        s.append(' %20s = %14.8f' %('gradpara', self.gradpara))
        s.append(' %20s = %14.8f' %('gradperp', self.gradperp))
        return '\n'.join(s)

# ---------------------------------------------------------------------
def tests():
    
    ins = Penaltyfunction('cej2004')
    print ins.info
    print ins.alpha
    print ins.beta 
    ins.sigma = 1.
    print ins.sigma
    ins.sigmaprime = np.array([1,2,3])
    print ins.sigmaprime
    print ins.thresh_step
    print ins.thresh_grad
    print ins.thresh_egap

def numeric_test():
    # k=2 of h2s_s0eq_bfgs_avdz.out
    class Sol:
        L = -398.65499724
        P = 0.00067023/3.5 
        gradpara =  -0.00337586
        gradperp = 0.02666273 

    Lprev = -398.64015016
    step = -0.01484708
    pf = Penaltyfunction('jpcb2008')
    pf.sigma =  -398.65566747 
    pf.gamma = 0.00205510 
    pf.sigmaprime = np.array([0.00000000,
                    0.00000000,
                    0.04553100,
                    0.00000000,
                   -0.01708175,
                   -0.02276600,
                    0.00000000,
                    0.01708175,
                   -0.02276600])
    pf.gammaprime = np.array([0.00000000,
                    0.00000000,
                   -0.05725700,
                    0.00000000,
                   -0.00216550,
                    0.02690700,
                    0.00000000,
                    0.00216550,
                    0.02690700])
    pf.update()
    fmt = '%10s: calc=%14.8f, sol=%14.8f' 
    print fmt %('L', pf.L, Sol.L)
    print fmt %('P', pf.P, Sol.P)
    print fmt %('step', pf.L-Lprev, step)
    print fmt %('gradpara', pf.gradpara, Sol.gradpara)
    print fmt %('gradperp', pf.gradperp, Sol.gradperp)
    print ''
    print pf.info

def test03():
    pf = Penaltyfunction('jpcl2011')
    print pf.alpha
    print pf.beta

if __name__ == '__main__':
    #numeric_test()

    test03()



