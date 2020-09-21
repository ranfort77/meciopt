#!/usr/bin/env python

class Accumulator:
    def __init__(self):
        self.sum = 0

    def __call__(self, *args):
        self.sum += sum(args)
        return self.sum

if __name__ == '__main__':
    acc = Accumulator()
    print acc(1,2,3,4,5) 
