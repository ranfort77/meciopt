#!/usr/bin/env python

class C:

    def __init__(self, geom=None):
        self.geom = geom
        self.basis = None
    
    def __getattr__(self, name):
        print '__getattr__(%s) called' %name
        raise AttributeError

    def __setattr__(self, name, value):
        print '__setattr__(%s)=%s called' %(name, value)
        self.__dict__[name] = value

#class New(object):
#    def __getattribute__(self, x):
#        print '__getattribute__ called', x 
#        return object.__getattribute__(self, x)


if __name__ == '__main__':
    ins = C()
    ins.geom = 'a b c d'
    ins.basis = 'aug-cc-pvdz'
    print ins.geom
    print ins.basis
    
    print ins.method
