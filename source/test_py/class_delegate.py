#!/usr/bin/env python

class C:

    def __init__(self):
        self.stack = [] 
    
    def __getattr__(self, name):
        print '%s' %name
        print self.stack
        print 
        return getattr(self.stack, name)


if __name__ == '__main__':
    ins = C()
    ins.append  (1)
    ins.append  (2)
    ins.append  (3)

