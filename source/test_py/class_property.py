#!/usr/bin/env python

class C(object):  # must inherit object for property
    __slots__ = ('_geom', '_basis')

    def __init__(self, geom=None):
        #self._geom = geom 
        self.geom = geom  # feel a little better...  # but the existence of 
                          # _geom is not explicitly known 
                          # if __slots__ is used, instance workspace can 
                          # be protected and attributes are explicitly known.
        self.basis = None

    @property
    def geom(self):
        print 'get called'
        return self._geom
    @geom.setter
    def geom(self, value):
        print 'set called'
        self._geom = value 
    @geom.deleter
    def geom(self):
        print 'del called'
        del self._geom

    @property
    def basis(self):
        return self._basis
    @basis.setter
    def basis(self, value):
        self._basis = value 
        
if __name__ == '__main__':
    ins = C('O\nH 1 1.0\nH 1 1.0 2 104.5\n')
    print ins.geom 
    ins.geom = 'O 0.0 0.0 0.0\nH 0.0 0.8 0.8\nH 0.8 0.0 0.8\n'
    print ins.geom
    del ins.geom
    #print ins.geom  # raise AttributeError

    # refer interal attribute
    ins.geom = 'geometry'
    print ins._geom

    # encapsulate by __slots__
    #ins.foo = 1  # raise AttributeError 
 
    ins.basis = '6-31g(d)'
    print ins.basis




