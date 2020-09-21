#!/usr/bin/env python

def func(kind='abc', *args):
    print kind
    print args 


def main():
    func()
    func('abc')
    func('def', 1, 2, 3, 4)


if __name__ == '__main__':
    main() 
