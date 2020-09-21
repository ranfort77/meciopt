#!/usr/bin/env python

class Des(object):
    def __init__(self, name):
        self.name = '_' + name 

    def __get__(self, instance, owner):
        print 'called __get__ %s' %self.name 
        if not hasattr(instance, self.name):
            return 'no attribute %s' %self.name
        #return instance.__dict__[self.name] 
        return getattr(instance, self.name) 

    def __set__(self, instance, value):
        print 'called __set__ %s = %s' %(self.name, value)
        #instance.__dict__[self.name] = value
        setattr(instance, self.name, value) 

    def __delete__(self, instance):
        print 'called __delete__ %s' %self.name
        #del instance.__dict__[self.name]
        delattr(instance, self.name)
        

class C(object):
    __slots__ = ('_geom', '_basis')  # encapsulation, but can't use __dict__
    geom = Des('geom')
    basis = Des('basis')

    def __init__(self, geom):
        self.geom = geom
        self.basis = None


if __name__ == '__main__':
    ins = C('O\nH 1 1.0\nH 1 1.0 2 104.5\n')
    print ins.geom
    #print ins.__dict__
    print 
    ins.geom = 'O 0.0 0.0 0.0\nH 0.0 0.8 0.8\nH 0.8 0.0 0.8\n'
    print ins.geom
    #print ins.__dict__
    print 
    del ins.geom
    print ins.geom
    #print ins.__dict__
    print 
    ins.basis = 'aug-cc-pvdz'
    
    ins2 = C('z-matrix')
    print ins2.geom
    ins2.basis = '6-311+G(d)'

    print ins.basis

    ins.mol = 'abc'
