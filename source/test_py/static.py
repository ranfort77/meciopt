#!/usr/bin/env python

class Foo:
    def __init__(self, inp):
        self.inp = inp

    @staticmethod
    def smtd1(s):
        print 'smtd1 called: %s' %s

    @staticmethod
    def smtd2(s):
        Foo.smtd1(s)    
        print 'smtd2 called: %s' %s



if __name__ == '__main__':
    #obj = Foo('abc') 
    Foo.smtd2('abcd')
