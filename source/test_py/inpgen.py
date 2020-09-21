#!/usr/bin/env python


class Procedure:
    def __init__(self, name):
        self.name = name
        self.commands = []
        self.directives = []

class Data:
    def __init__(self):
        self.numbers = []
        self.expressions = []
        self.strings = []  
        self.sep = [',', ' ']

class Option:
    def __init__(self, name):
        self.name = name
        self.value = None


class Directive:
    def __init__(self, name):
        self.name = name
        self.data = []
        self.options = []


class Command:
    def __init__(self, name):
        self.name = name
        self.options = []
        self.directives = []
        self.data = []

main():
    pass



if __name__ == '__main__':
    main()
