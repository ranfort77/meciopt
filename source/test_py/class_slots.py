#!/usr/bin/env python

class Person(object):
    __slots__ = ['name', 'tel']

    def method(self):
        print 'method called'


if __name__ == '__main__':
    a = Person()
    a.name = 'abc'
    a.tel = 'def' 
    a.method()
    a.value = 1
