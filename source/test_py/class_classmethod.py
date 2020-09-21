#!/usr/bin/env python


class Foo:
    def mtd(self, x):
        print 'Foo:mtd called, self=%s, x=%s' %(self, x)


class Bar:
    def mtd(cls, x):
        print 'Bar:mtd called, cls=%s, x=%s' %(cls, x)
    mtd = classmethod(mtd)


class Baz:
    @classmethod
    def mtd(cls, x):
        print 'Baz:mtd called, cls=%s, x=%s' %(cls, x)

if __name__ == '__main__':
    ins = Foo() 
    ins.mtd('abc')

    Bar.mtd('abc')

    Baz.mtd('abc')
