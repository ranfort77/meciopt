#!/usr/bin/env python


class Foo:
    def mtd(self, x):
        print 'Foo:mtd called, self=%s, x=%s' %(self, x)


class Bar:
    def mtd(x):
        print 'Bar:mtd called, x=%s' %x
    mtd = staticmethod(mtd)


class Baz:
    @staticmethod
    def mtd(x):
        print 'Baz:mtd called, x=%s' %x

if __name__ == '__main__':
    ins = Foo() 
    ins.mtd('abc')

    Bar.mtd('abc')

    Baz.mtd('abc')
