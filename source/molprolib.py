#!/usr/bin/env python
"""
Molpro library
test : Python, 2.7.14 

version : 1.0

author : Heesun, An

TODO : check normal termination
"""

import os
import sys
import re
import fileinput
import pprint
import subprocess
import shutil


class MolproLibException(Exception):
    pass 

class OutputException(MolproLibException):
    pass

class InputException(MolproLibException):
    pass

class InputGeneratorException(MolproLibException):
    pass

USER = os.environ['USER']
DEFAULT_SCRDIR = '/work/tmp/%s' %USER 

class Input:
    def __init__(self, inp, out=None, scrdir=None):
        """
        inp and out include abs-path or rel-path 
        ex) inp = '/home/user/ahs/jobs/h2s.inp'
            out = 'tests/h2s_first.out'
            scrdir = '/work/tmp'
        """
        self.inp = os.path.abspath(inp) 
        fp = open(self.inp, 'r')
        self.inpstr = fp.read()
        fp.close()
        self.curdir = os.getcwd()
        if out is None: 
            out = os.path.splitext(os.path.basename(self.inp))[0] + '.out'
            self.out = '/'.join([self.curdir, out])
            self.wfudir = self.curdir
	else: 
            self.out = out
	    self.wfudir = os.path.dirname(out)
        self.curpid = os.getpid()
        if scrdir is None:
            self.scrdir = DEFAULT_SCRDIR
        else:
            self.scrdir = scrdir
        if os.access(self.scrdir, os.F_OK): # is exist
            self.scrdir = self.scrdir + '/molpro_%d' %self.curpid 
        else:
            errmsg = 'scratch directory "%s" does not exist' %self.scrdir
            raise InputException(errmsg)      
        self.cmd = 'molpro -s --no-xml-output -d %s -W %s -o %s %s' \
            %(self.scrdir, self.wfudir, self.out, self.inp)

    def __repr__(self):
        s = []
        s.append('%15s: %-60s' %('input', self.inp))
        s.append('%15s: %-60s' %('output', self.out))
        s.append('%15s: %-60s' %('scratch dir', self.scrdir))
        s.append('%15s: %-s' %('run command', self.cmd))
        return '\n'.join(s)

    def run(self):
        p = subprocess.Popen(self.cmd, shell=True, 
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        p.wait()
        res = {'stdout':p.stdout.read(), 'stderr':p.stderr.read()}
        p = subprocess.Popen('rm -rf %s' %self.scrdir, shell=True)
        p.wait()
        return res 

    def analize_lines_of_inpstr(self):
        """This function analyzes each line of the molpro input string
        to distinguish the first keyword and parenthesis block. 
        It should only be used on simple molpro inputs. 
        Not only do-enddo, if-else, and proc-endproc block but also
        parenthesis({}) of two or more levels are not considered. 
        """
        lines = self.inpstr.strip().split('\n')
        delim = re.compile('(,|=|;)')
        clines = [None for n in range(len(lines))]
        comments = [None for n in range(len(lines))]
        keywords = [None for n in range(len(lines))]
        inblocks = [None for n in range(len(lines))]
        ifcases = [None for n in range(len(lines))]
        for n, line in enumerate(lines):
            # split input and comment
            if '!' in line:
                foo = line.strip().split('!')
                cline = foo[0].strip() 
                comment = foo[1] 
            else:
                cline = line.strip()
                comment = ''
            clines[n] = cline
            comments[n] = comment
            clineList = [e.strip() for e in delim.split(cline)] 
            # analize first keyword and {} block
            if cline == '':
                ifcases[n] = 'empty line'
                keywords[n] = None
                inblocks[n] = False
            elif cline[0] == '{' and cline[-1] == '}':
                ifcases[n] = '{hf} or {hf,...}'
                keywords[n] = clineList[0][1:].replace('}','').strip()
                inblocks[n] = False
            elif cline[0] == '{' and cline[-1] != '}':
                ifcases[n] = '{hf,...'
                keywords[n] = clineList[0][1:].strip()
                inblocks[n] = True
            elif '}' in line and inblocks[n-1] is True:
                ifcases[n] = '...}'
                keywords[n] = None
                inblocks[n] = False
            elif '{' in cline and '}' in cline and clineList[1] == '=':
                ifcases[n] = 'geometry = {...}'
                keywords[n] = clineList[0]
                inblocks[n] = False
            elif '{' in cline and '}' not in cline and clineList[1] == '=':
                ifcases[n] = 'geometry = {...'
                keywords[n] = clineList[0]
                inblocks[n] = True
            elif cline[-1] == ';' and len(clineList) == 3 and inblocks[n-1] is False:
                ifcases[n] = 'hf;' # this case clineList = ['hf',';','']
                keywords[n] = clineList[0]
                inblocks[n] = False 
            elif len(clineList) == 1 and inblocks[n-1] is False: 
                ifcases[n] = '---'
                keywords[n] = clineList[0]
                inblocks[n] = False
            elif inblocks[n-1] is True:
                ifcases[n] = 'directive between { and }'
                keywords[n] = None
                inblocks[n] = True
            elif clineList[1] == '=':
                ifcases[n] = 'R1 = 1.0 bohr'
                keywords[n] = clineList[0]
                inblocks[n] = False
            elif clineList[1] == ',':
                ifcases[n] = 'memory,128,m   or   ***,h2s'
                keywords[n] = clineList[0]
                inblocks[n] = False
            else:
                errmsg = 'Unknown case line input "%s"' %cline
                raise InputException(errmsg)      
            #endif    
        return clines, comments, keywords, inblocks, ifcases

    def get_blocks_inpstr(self):
        clines, comments, keywords, inblocks, ifcases = \
                self.analize_lines_of_inpstr()
        data = zip(inblocks, keywords, clines)
        it = data.__iter__()
        blocks = []
        for a, b, c in it:
            if c == '':  # empty line
                continue
            if a is True:
                blocks.append({'keyword':b, 'inpstr':[c]})
                for a, b, c in it:
                    blocks[-1]['inpstr'].append(c)  
                    if a is False:
                        break
            else:
                blocks.append({'keyword':b, 'inpstr':[c]})
        for n in range(len(blocks)):
            blocks[n]['inpstr'] = '\n'.join(blocks[n]['inpstr'])
        return blocks
    

# ---------------------------------------------------------------------
REGEXP_ATOM = '[A-Z][a-z]?'
REGEXP_FLOAT = '[-+]?[0-9]*\.?[0-9]+'
REGEXP_FLOAT_SCINOTATION = '%s[ED]?[-+]?[0-9]*' %REGEXP_FLOAT
REGEXP_STATESYM = '[0-9]+.[0-8]'
REGEXP_LENGTH = '(\d+)-\s?(\d+)\s+(%s)' %REGEXP_FLOAT
REGEXP_ANGLE = '(\d+)-\s?(\d+)-\s?(\d+)\s+(%s)' %REGEXP_FLOAT
REGEXP_ENERGY = \
    lambda method: '!%s STATE (%s) Energy\s+(%s)' \
        %(method, REGEXP_STATESYM, REGEXP_FLOAT)
REGEXP_DIPOLE = \
    lambda method: '!%s STATE (%s) Dipole moment\s+(%s)\s+(%s)\s+(%s)' \
        %((method, REGEXP_STATESYM) + 3*(REGEXP_FLOAT,))
REGEXP_DENMAT = lambda method: '!%s.+\<(%s)\|[DM]*([XYZ]?)\|(%s)\>\s+(%s)' \
        %(method, REGEXP_STATESYM, REGEXP_STATESYM, REGEXP_FLOAT) 
REGEXP_BASIS = '[A-Z0-9+,\-()]+'
REGEXP_VARIABLE = 'SETTING ([A-Z][A-Z0-9_]*)\(?(\d*)\)?\s+=\s+(%s)?\s+(%s)?' \
        %(REGEXP_FLOAT_SCINOTATION, REGEXP_BASIS)

class Output:
    def __init__(self, out):
        self.out = out
        self.version = self._get_version() 
        if self.version is None:
            errmsg = '"%s" is not a molpro output' %self.out
            raise OutputException(errmsg)     

    def __repr__(self):
        """
        pretty print data of all programs
        """
        alldata = self.get('all')
        s = []
        width = 60
        s.append('PARSING MOLPRO OUTPUT'.center(width))
        s.append('%s' %(width*'*')) 
        s.append('%10s: %-48s' %('file name', self.out))
        s.append('%10s: %-48s' %('version', self.version))
        for pdat in alldata:
           s.append('%s' %(width*'*')) 
           s.append('Line Number : %d' %pdat['_linenumber'])
           for key in pdat['_orderedkey']:
               s.append('%s' %key)
               s.append(pprint.pformat(pdat[key], width=width))
               s.append('') 
        return '\n'.join(s)

    def check_normalterm(self):
        normalterm = False
        fp = open(self.out, 'r')
        for line in fp:
            if line[1:25] == 'Variable memory released': 
                normalterm = True
                break
        fp.close()
        return normalterm

    def get(self, token):
        if token == 'all':
            if not hasattr(self, '_alldata'):
                self._alldata = self._parse_allprograms()
            return self._alldata
        elif token == 'var':
            if not hasattr(self, '_variables'):
                self._variables = self._parse_variables()
            return self._variables
        else:
            errmsg = '"%s" Unknown' %(token)
            raise OutputException(errmsg)     

    def get_realtime(self, token=None):
        alldata = self.get('all') 
        program = []
        cumtime = []
        eachtime = []
        for n, progdata in enumerate(alldata):
            program.append(progdata['program'])
            cumtime.append(progdata['_realtime'])
            if n == 0: 
                eachtime.append(progdata['_realtime'])
            else:
                eachtime.append(cumtime[-1] - cumtime[-2])
        if token is None: 
            return cumtime[-1] # total real time
        else:
            if token == 'each':
                return program, eachtime
            else:
                errmsg = '"%s" Unknown' %(token)
                raise OutputException(errmsg)     

    def get_xyzcoords(self, token):
        alldata = self.get('all') 
        indice = []
        for n, progdata in enumerate(alldata):
            if 'xyz' in progdata:
                 indice.append(n) 
        if token == 'initial':
            progindex = indice[0]
        elif token == 'final':
            progindex = indice[-1]
        elif type(token) is int:
            if token <= len(indice):
                progindex = indice[token-1]
            else:
                errmsg = 'xyz data exist %d. You entered %d' %(len(indice), token)
                raise OutputException(errmsg)     
        atom = alldata[progindex]['atom']
        xyz = alldata[progindex]['xyz']
        return atom, xyz

    def get_final_energies(self):
        alldata = self.get('all')
        energy = {}
        for n, progdata in enumerate(alldata):
            if 'energy' in progdata:
                energy[progdata['program']] = {'state':progdata['state'], 
                                               'spin':progdata['spin'],
                                               'energy':progdata['energy']} 
        return energy

    def get_final_symmetry_elements(self):
        alldata = self.get('all')
        symel = None
        for progdata in alldata:
            if 'symel' in progdata: 
                symel = progdata['symel']
        return symel

    def get_final_pointgroup(self):
        alldata = self.get('all')
        pg = None
        for progdata in alldata:
            if 'pg' in progdata: 
                pg = progdata['pg']
        return pg 

    def get_formtable(self, start, end, header, nskip, formats):
        """
        compares each line[start:end] of molpro output with header,
        if it is the same, then skip n lines, 
        parses table data according to format.

        For example, to get ATOMIC COORDINATES data in SEWARD,
        >> obj = Output(molpro_output) 
        >> header = 'ATOMIC COORDINATES'
        >> formats = 'd s f fff' 
        >> table = obj.get_formattedtable(1, 19, header, 3, formats)
        """ 
        table = None
        fp = open(self.out, 'r')
        for line in fp:
            table = self._parse_table(fp, line[start:end], 
                                      header, nskip, formats)
            if table is not None:
                break
        fp.close() 
        return table

    def _get_version(self):
        version = None
        fp = open(self.out, 'r')
        for line in fp:
            if line[41:72] == '***  PROGRAM SYSTEM MOLPRO  ***': 
                robj = re.compile('Version ([0-9]{4}\.[0-9]{1}) linked')
                for line in fp:
                    fobj = robj.search(line) 
                    if fobj is not None:
                        version = fobj.group(1)
                        break
                #endfor
                break
        fp.close() 
        return version

    def _parse_variables(self):
        re_var = re.compile(REGEXP_VARIABLE)
        groups = []
        fp = open(self.out, 'r')
        for line in fp:
            if line[1:8] == 'SETTING':
                # 'SETTING E1(1)          =      -306.63455068  AU' 
                # >> [('E1', '1', '-306.63455068', 'AU')] 
                # 'SETTING BASIS          =    6-311++G(D,P)'
                # >> [('BASIS', '', '', '6-311++G(D,P)')]
                groups.append(re_var.findall(line)[0]) 
        fp.close()
        # check variable size
        varsize = {} 
        for group in groups:
            varname = group[0]
            index = group[1]
            if varname not in varsize:
                varsize[varname] = 1
            else:
                if index:
                    varsize[varname] = int(index) 
        # pre-allocation
        variables = {} 
        for varname in varsize:
            variables[varname] = [None for n in range(varsize[varname])] 
        # assign value
        for group in groups:
            varname = group[0]
            index = group[1] # '' or integer
            value = group[2].replace('D', 'E')
            basis = group[3] # or unit 
            if index:
                variables[varname][int(index)-1] = float(value) 
            else:
                if value:
                    variables[varname] = float(value) 
                else:
                    variables[varname] = basis 
        return variables

    def _parse_allprograms(self):
        """
        parses the info defined for each program in the molpro output
        """
        seqdata = []
        fip = fileinput.input(self.out, mode='r')
        for line in fip:
            program = None
            if line[1:10] == 'PROGRAM *':
                program = line.split()[2] 
            if line[1:24] == 'Construct non-adiabatic':
                program = 'DDR'
            if program:
                progdata = self._program_parser(fip, program)
                seqdata.append(progdata)
        fip.close()
        return seqdata

    def _program_parser(self, fip, program):
        data = {'_linenumber':fip.lineno(), '_orderedkey':['program'], 
                'program':program, '_realtime':None} 
        # -------------------------------------------------------------
        if program == 'SEWARD':
            keys = ['symel', 'pg', 'atom', 'charge', 'xyz', 'length', 'angle']
            vals = [[], None, None, None, None, dict(), dict()]
            data['_orderedkey'].extend(keys)
            data.update(zip(keys, vals)) 
            re_len = re.compile(REGEXP_LENGTH)
            re_ang = re.compile(REGEXP_ANGLE)
            for line in fip:
                if line[1:13] == 'REAL TIME  *':
                    data['_realtime'] = float(line.split()[3]) 
                    break
                if line[1:18] == 'Symmetry elements':
                    data['symel'] = [e.strip() 
                                 for e in line.split(':')[1].split(',')]
                if line[1:12] == 'Point group':
                    data['pg'] = line.split()[2] 
                table = self._parse_table(fip, line[1:19], 
                        'ATOMIC COORDINATES', 3, 'd s f fff') 
                if table is not None:
                    data['atom'] = table[1] 
                    data['charge'] = table[2] 
                    data['xyz'] = table[3]
                if line[1:13] == 'Bond lengths':
                    while 1:
                        line = fip.next()
                        line = fip.next()
                        groups = re_len.findall(line)
                        if not groups:
                            break
                        for group in groups:
                            key = (int(group[0]), int(group[1])) 
                            data['length'][key] = float(group[2])
                        line = fip.next()
                if line[1:12] == 'Bond angles':
                    while 1:
                        line = fip.next()
                        line = fip.next()
                        groups = re_ang.findall(line)
                        if not groups:
                            break
                        for group in groups:
                            key = (int(group[0]), int(group[1]), int(group[2])) 
                            data['angle'][key] = float(group[3])
        # -------------------------------------------------------------
        elif program == 'RHF-SCF':
            keys = ['occ', 'spin', 'state', 'energy', 'dipole']
            vals = [None, [], [], [], []]
            data['_orderedkey'].extend(keys)
            data.update(zip(keys, vals)) 
            re_ene = re.compile(REGEXP_ENERGY('RHF'))
            re_dip = re.compile(REGEXP_DIPOLE('RHF'))
            for line in fip:
                if line[1:13] == 'REAL TIME  *':
                    data['_realtime'] = float(line.split()[3]) 
                    break
                if line[60:73] == 'SPIN SYMMETRY':
                    # molpro v2008 is "SPIN SYMMETRY=Singlet"
                    # v2012 & v2015 are "SPIN SYMMETRY: Singlet" 
                    groups = re.findall('SPIN SYMMETRY[:=]\s*(\w+)', line)
                    data['spin'].append(groups[0]) 
                if line[1:17] == 'Final occupancy:':
                    data['occ'] = [int(e) for e in line.split()[2:]] 
                groups = re_ene.findall(line)
                if groups:
                    data['state'].append(groups[0][0])
                    data['energy'].append(float(groups[0][1]))
                groups = re_dip.findall(line)
                if groups:
                    data['dipole'].append([float(e) for e in groups[0][1:]])
        # -------------------------------------------------------------
        elif program == 'MULTI':
            keys = ['closed', 'active', 'occ', 'spin', 'state', 
                    'energy', 'dipole']
            vals = [None, None, None, [], [], [], []]
            data['_orderedkey'].extend(keys)
            data.update(zip(keys, vals)) 
            re_ene = re.compile(REGEXP_ENERGY('MCSCF'))
            re_dip = re.compile(REGEXP_DIPOLE('MCSCF'))
            for line in fip:
                if line[1:13] == 'REAL TIME  *':
                    data['_realtime'] = float(line.split()[3]) 
                    break
                if line[1:17] == 'Number of closed':
                    groups = re.findall('\d+', line)
                    data['closed'] = [int(e) for e in groups[1:]] 
                if line[1:17] == 'Number of active':
                    groups = re.findall('\d+', line)
                    data['active'] = [int(e) for e in groups[1:]] 
                if line[31:45] == 'Spin symmetry=':
                    nstates = int(fip.next().split(':')[1])
                    groups = re.findall('Spin symmetry=(\w+)', line)
                    for ns in range(nstates): 
                        data['spin'].append(groups[0])
                groups = re_ene.findall(line)
                if groups:
                    data['state'].append(groups[0][0])
                    data['energy'].append(float(groups[0][1]))
                groups = re_dip.findall(line)
                if groups:
                    data['dipole'].append([float(e) for e in groups[0][1:]])
            #endfor
            data['occ'] = map(lambda x, y: x + y, 
                              data['closed'], data['active'])
        # -------------------------------------------------------------
        elif program == 'CI':
            fip.next()
            fip.next()
            line = fip.next()
            if line[3:20] == 'Transition moment':
                istr = True 
                keys = ['state', 'overlap', 'dmx', 'dmy', 'dmz'] 
                vals = [[], [], [], [], []]
                re_dm = re.compile(REGEXP_DENMAT('MRCI'))
                dmtemp = []
            else: 
                istr = False
                keys = ['core', 'active', 'spin', 'state', 
                        'energy', 'dipole']
                vals = [None, None, [], [], [], []]
            data['_orderedkey'].extend(keys)
            data.update(zip(keys, vals)) 
            re_ene = re.compile(REGEXP_ENERGY('MRCI'))
            re_dip = re.compile(REGEXP_DIPOLE('MRCI'))
            for line in fip:
                if line[1:13] == 'REAL TIME  *':
                    data['_realtime'] = float(line.split()[3]) 
                    break
                if line[1:20] == 'Reference symmetry:':
                    spin = line.split()[3] 
                if line[1:15] == 'Number of core':
                    groups = re.findall('\d+', line)
                    data['core'] = [int(e) for e in groups[1:]] 
                if line[1:17] == 'Number of active':
                    groups = re.findall('\d+', line)
                    data['active'] = [int(e) for e in groups[1:]] 
                groups = re_ene.findall(line)
                if groups:
                    data['state'].append(groups[0][0])
                    data['energy'].append(float(groups[0][1]))
                groups = re_dip.findall(line)
                if groups:
                    data['dipole'].append([float(e) for e in groups[0][1:]])
                if istr:
                    #'!MRCI overlap           <2.2||1.2>     0.006184435466 '
                    # >> ('2.2', '', '1.2', '0.006184435466') 
                    # '!MRCI trans          <2.2|DMZ|1.2>     1.607403924322'
                    # >> ('2.2', 'Z', '1.2', '1.607403924322')
                    groups = re_dm.findall(line)
                    if groups:
                        dmtemp.append(groups[0])
            #endfor
            if istr:
                states = [] 
                for dminfo in dmtemp:
                    if dminfo[0] not in states:
                        states.append(dminfo[0])
                    if dminfo[2] not in states:
                        states.append(dminfo[2])
                ns = len(states)
                ind = dict(zip(states, range(ns)))
                ove = [[None for m in range(ns)] for n in range(ns)] 
                dmx = [[None for m in range(ns)] for n in range(ns)] 
                dmy = [[None for m in range(ns)] for n in range(ns)] 
                dmz = [[None for m in range(ns)] for n in range(ns)] 
                for dminfo in dmtemp:
                    ii = dminfo[0]
                    jj = dminfo[2]
                    if dminfo[1] == '':
                        ove[ind[ii]][ind[jj]] = float(dminfo[3])
                    elif dminfo[1] == 'X':
                        dmx[ind[ii]][ind[jj]] = float(dminfo[3])
                    elif dminfo[1] == 'Y':
                        dmy[ind[ii]][ind[jj]] = float(dminfo[3])
                    elif dminfo[1] == 'Z':
                        dmz[ind[ii]][ind[jj]] = float(dminfo[3])
                data['state'] = states
                data['overlap'] = ove
                data['dmx'] = dmx 
                data['dmy'] = dmy 
                data['dmz'] = dmz 
            else:
                nstates = len(data['state'])
                for ns in range(nstates): 
                    data['spin'].append(spin)
        # -------------------------------------------------------------
        elif program == 'DDR':
            keys = ['nacme']
            vals = [None]
            data['_orderedkey'].extend(keys)
            data.update(zip(keys, vals)) 
            for line in fip:
                if line[1:13] == 'REAL TIME  *':
                    data['_realtime'] = float(line.split()[3]) 
                    break
                if line[1:13] == '!Total NACME':
                    data['nacme'] = float(line.split()[2])
        # -------------------------------------------------------------
        elif program == 'CCSD':
            keys = ['spin', 'state', 'energy']
            vals = [[], [], []]
            data['_orderedkey'].extend(keys)
            data.update(zip(keys, vals)) 
            for line in fip:
                if line[1:13] == 'REAL TIME  *':
                    data['_realtime'] = float(line.split()[3]) 
                    break
                table = self._parse_table(fip, line[2:15], 
                        'Final Results', 3, 's f f') 
                if table is not None:
                    data['state'] = table[0] 
                    data['energy'] = table[2] 
        # -------------------------------------------------------------
        elif program == 'RESTART':
            keys = []
            vals = []
            data['_orderedkey'].extend(keys)
            data.update(zip(keys, vals)) 
            for line in fip:
                if line[1:13] == 'REAL TIME  *':
                    data['_realtime'] = float(line.split()[3]) 
                    break
        # -------------------------------------------------------------
        elif program == 'TEMPLET':
            keys = []
            vals = []
            data['_orderedkey'].extend(keys)
            data.update(zip(keys, vals)) 
            for line in fip:
                if line[1:13] == 'REAL TIME  *':
                    data['_realtime'] = float(line.split()[3]) 
                    break
        # -------------------------------------------------------------
        else:
            errmsg = 'program "%s" parser is not defined' %(program)
            raise OutputException(errmsg)     
        # -------------------------------------------------------------
        return data

    def _parse_table(self, fp, subline, header, nskip, recordformats):
        """
        Example: 
          ATOMIC COORDINATES
          NR  ATOM    CHARGE       X              Y              Z
          1  S      16.00    0.000000000    0.000000000   -0.118252046
          2  H       1.00    0.000000000    2.145200000    1.880647954
          3  H       1.00    0.000000000   -2.145200000    1.880647954

        >> self._parse_table(fp, line[1:13], 'ATOMIC COORD', 1, 'd s f fff') 
        return, 
          ([1, 2, 3], ['S', 'H', 'H'], [16.0, 1.0, 1.0], 
           [[0.0, 0.0, -0.118252046], 
           [0.0, 2.1452, 1.880647954], 
           [0.0, -2.1452, 1.880647954]])
        """ 
        table = None
        if subline == header: 
            conv = {'f':lambda x: float(x), 
                    'd':lambda x: int(x), 
                    's':lambda x: x}
            for n in range(nskip):
                fp.next() 
            rfs = recordformats.split()
            table = [[] for n in range(len(rfs))]
            for line in fp:
                line = line.strip()
                if not line:
                    break
                line = line.split()
                sp = 0
                for k, rf in enumerate(rfs):
                    if len(rf) == 1:
                        table[k].append(conv[rf](line[sp]))
                        sp += 1
                    else:
                        foo = []
                        for cf in rf:
                            foo.append(conv[cf](line[sp])) 
                            sp += 1
                        table[k].append(foo) 
            #endfor
        return table 
   

class InputGenerator:
    def __init__(self, inpname, inpstr, prevwfuname=None):
        """If prevwfuname is None, then set file,2,...,new """
        self.inpname = inpname
        prefix = os.path.splitext(os.path.basename(self.inpname))[0]
        self.wfuname = prefix + '.wfu'
        self.prevwfuname = prevwfuname
        self.inpstr = inpstr

    def make_input(self):
        """make input file and copy wfu file"""
        self.inpstr = self._sync_wfuname()
        fp = open(self.inpname, 'w')
        fp.write(self.inpstr)
        fp.close()
        if self.prevwfuname is not None:
            shutil.copyfile(self.prevwfuname, self.wfuname)

    def _sync_wfuname(self):
        if self.prevwfuname is None:
            wfuline = 'file,2,%s,new' %self.wfuname
        else:
            wfuline = 'file,2,%s,old' %self.wfuname
        # check existence of ***, memory, and file2 lines
        # and insert wfu-line
        lines = self.inpstr.split('\n') 
        idx_sss = None
        idx_mem = None
        idx_file2 = None
        for idx, line in enumerate(lines): 
            line = line.strip().split('!')[0] # delete comment string
            line = ''.join(line.split())  # delete white space
            if '***,' in line:
                idx_sss = idx 
            elif 'memory,' in line:
                idx_mem = idx 
            elif 'file,2' in line:
                idx_file2 = idx 
            else:
                continue
        if type(idx_file2) == int:
            lines[idx_file2] = wfuline
        else:
            if type(idx_mem) == int:
                lines.insert(idx_mem + 1, wfuline)
            else:
                if type(idx_sss) == int:  # because idx_sss is possible 0 
                    lines.insert(idx_sss + 1, wfuline)
                else:
                    lines.insert(0, wfuline)
        inpstr = '\n'.join(lines) 
        # make input file
        return inpstr

    @staticmethod
    def get_xyz_single(xyz, atom, symel, energyinpstr, table=False):
        """>>> CLS.gen_xyz_single([[1, 2, 3],[3, 4, 5]], ['S', 'H'], ['X'], 
               'hf\nmulti', True) """
        natom = len(xyz)
        s = []
        s.append('***,FromInputGenerator')
        s.append('memory,128,m')
        s.append('')
        s.append(InputGenerator._var_ctxyz(xyz))
        s.append('')
        s.append(InputGenerator._symmetry(symel))
        s.append(InputGenerator._geometry_cxyz(atom))
        s.append('natom = %d' %natom)
        s.append('')
        s.append(InputGenerator._proc_setct())
        s.append('')
        s.append('setct')
        s.append(energyinpstr)
        if table:
            s.append('{table,energy')
            s.append('  title,single point energy')
            s.append("  format,'f16.8'}")
        return '\n'.join(s)

    @staticmethod
    def get_xyz_numdif(xyz, atom, symel, energyinpstr, nstates, difdr):
        s = []
        s.append('do n = 1,natom')
        s.append(InputGenerator._var_state_energy(nstates, 'st%dxc(n)', 'energy(%d)'))
        s.append(InputGenerator._var_state_energy(nstates, 'st%dyc(n)', 'energy(%d)'))
        s.append(InputGenerator._var_state_energy(nstates, 'st%dzc(n)', 'energy(%d)'))
        s.append('enddo')
        s.append('') 
        s.append('dr = %8.4f bohr' %difdr)
        s.append('do n = 1,natom')
        for axis in ['X', 'Y', 'Z']:
            for sign in ['+', '-']:
                s.append('setct')
                if axis in symel:
                    varformat = 'st%%d%s%s(n)' %(axis.lower(), 
                                  {'+':'p', '-':'m'}[sign])
                    s.append(InputGenerator._var_state_energy(nstates, varformat))
                    s.append('')
                else:
                    s.append('C%s(n) = C%s(n) %s dr' %(axis, axis, sign))    
                    s.append(energyinpstr)
                    varformat = 'st%%d%s%s(n)' %(axis.lower(), 
                                  {'+':'p', '-':'m'}[sign])
                    valueformat = 'energy(%d)'
                    s.append(InputGenerator._var_state_energy(nstates, varformat, 
                                                               valueformat))
                    s.append('')
        s.append('enddo')
        s.append('')
        tables = []
        for ns in range(nstates):
            varlist = []
            for axis in ['x', 'y', 'z']:
                for pos in ['m', 'c', 'p']:
                    varlist.append('st%d%s%s' %(ns+1, axis, pos))
            tables.append('{table,%s' %','.join(varlist))
            tables.append('  title,state %d energy' %(ns+1))
            tables.append("  format,'%df16.8'}" %len(varlist))
            tables.append('')
        s.append('\n'.join(tables))
        return '\n'.join(s)

    @staticmethod
    def _symmetry(symel):
        if len(symel) == 0:
            s = 'symmetry,nosym'
        else:
            s = 'symmetry,%s' %','.join(symel)
        return s.lower()

    @staticmethod
    def _var_ctxyz(xyz):
        """>>> CLS._var_ctxyz([[1,2,3],[4,5,6],[7,8,9]]) """ 
        natom = len(xyz)
        ctxyz = []
        for n in range(natom): 
            s = []
            lb = n + 1
            s.append('CTX(%2d)= %14.8f bohr' %(lb, xyz[n][0])) 
            s.append('CTY(%2d)= %14.8f bohr' %(lb, xyz[n][1])) 
            s.append('CTZ(%2d)= %14.8f bohr' %(lb, xyz[n][2])) 
            '; '.join(s)
            ctxyz.append('; '.join(s))
        return '\n'.join(ctxyz)

    @staticmethod
    def _geometry_cxyz(atom):
        """>>> CLS._geometry_cxyz(['S', 'H', 'H']) """ 
        natom = len(atom)
        s = []
        s.append('geometry={')
        for n in range(natom):
            lb = n + 1
            s.append('%s,, CX(%2d), CY(%2d), CZ(%2d)' %(atom[n], lb, lb, lb))
        s.append('}')
        return '\n'.join(s)

    @staticmethod
    def _proc_setct():
        """ >>> CLS._proc_setct() """
        s = []
        s.append('proc setct')
        s.append('  do k=1,natom')
        s.append('    CX(k) = CTX(k); CY(k) = CTY(k); CZ(k) = CTZ(k);')
        s.append('  enddo')
        s.append('endproc')
        return '\n'.join(s)

    @staticmethod
    def _var_state_energy(nstates, varformat, valueformat=None):
        """ >>> CLS._var_state_energy(5, 'st%dyp(n)') 
            >>> CLS._var_state_energy(5, 'st%dyp(n)', 'energy(%d)') """
        s= []
        for n in range(nstates):
            lb = n + 1
            if valueformat is None:
                s.append((varformat %lb) + ' = 0.0') 
            else:
                s.append((varformat %lb) + ' = ' + (valueformat %lb)) 
        return '\n'.join(s)


# ---------------------------------------------------------------------
def usage_example_Output():
    fnames = ['test_outs/t001_phenol_v08.out', 
              'test_outs/t002_nh3cl_c2v_v12.out',
              'test_outs/t003_h2s_cs_v15.out',
              'test_outs/t004_ethylene_nac_v15.out',
              'test_outs/t005_h2s_program_restart.out']
    fname = sys.argv[1]

    obj = Output(fname)

    print 'get main infomations in each program'
    allprog = obj.get('all')
    print allprog
    for prog in allprog:
        print '%15s: %-14.2f sec' %(prog['program'], prog['_realtime'])
    print 

    print 'pretty print main infomations by __str__ method'
    print obj
    print

    print 'get variables in molpro output'
    variables = obj.get('var')
    pprint.pprint(variables, width=60)
    print 

    print 'get initial atomic coordinates'
    atom, xyz = obj.get_xyzcoords('initial')
    pprint.pprint([atom, xyz], width=60)
    print 'get final atomic coordinates'
    atom, xyz = obj.get_xyzcoords('final')
    pprint.pprint([atom, xyz], width=60)
    print 'get n-th atomic coordinates'
    atom, xyz = obj.get_xyzcoords(1)
    pprint.pprint([atom, xyz], width=60)
    print 

    print 'get any table data by get_formtable method' 
    index = fnames.index(fname)
    if index == 0: 
        format = 'f f f f f f f f'
        table = obj.get_formtable(0, 8, '   R(OH)', 0, format)
    elif index == 1:
        format = 'ff ffff ff ff'
        table = obj.get_formtable(0, 14, ' Non-adiabatic', 2, format)
    elif index == 2:
        format = 'f f f' 
        table = obj.get_formtable(0, 10, ' h2s table', 2, format)
    elif index == 3:
        format = 'fff fff fff' 
        table = obj.get_formtable(0, 12, ' s0-s1 nacme', 2, format)
    elif index == 4:
        format = 'fff fff fff' 
        table = obj.get_formtable(0, 15, ' state 5 energy', 2, format)
    pprint.pprint(table, width=60) 

    print obj.get_final_energies()

def t01_Output():
    out = sys.argv[1]
    outobj = Output(out)
    print outobj.get_final_pointgroup()
    print outobj.get_final_symmetry_elements()

def t02_Output():
    out = sys.argv[1]
    outobj = Output(out)
    a = outobj.get_final_energies()
    print a

def usage_example_Input():
    # out and wfu are created in the current path 
    #obj = Input('tests/t005_h2s.inp')
    #obj.run()

    # out and scrdir can modify 
    outdir = 'temp' 
    os.mkdir(outdir)
    out = '%s/t005_h2s_hf.out' %outdir
    scrdir = '/work/tmp'
    obj = Input('tests/t005_h2s.inp', out, scrdir)
    print obj.curpid
    obj.run()

def t01_analize_inpstr():
    obj = Input(sys.argv[1])
    clines, comments, keywords, inblocks, ifcases = obj.analize_lines_of_inpstr()
    for n in range(len(clines)):
        print 'comments: %s' %comments[n]
        print 'clines  : %s' %clines[n]
        print 'keywords: %s' %keywords[n]
        print 'inblocks: %s' %inblocks[n]
        print 'ifcases : %s' %ifcases[n]
        print '-------------' 
    blocks = obj.get_blocks_inpstr()
    for block in blocks:
        print 'keyword = %s' %block['keyword']
        print '%s' %block['inpstr']
        print 


if __name__ == '__main__':
    #usage_example_Output()
    #t01_Output()
    #usage_example_Input()
    #t01_analize_inpstr()
    #t02_Output()

    #print InputGenerator._var_ctxyz([[1,2,3],[4,5,6],[7,8,9]]) 
    #print InputGenerator._geometry_cxyz(['S', 'H', 'H']) 
    #print InputGenerator._proc_setct()
    #print InputGenerator._var_state_energy(5, 'st%dyp(n)') 
    #print InputGenerator._var_state_energy(5, 'st%dyp(n)', 'energy(%d)')
    #print InputGenerator.get_xyz_single([[1, 2, 3],[3, 4, 5]], ['S', 'H'], ['X'], 'hf\nmulti')
    #print InputGenerator.get_xyz_numdif([[1, 2, 3],[3, 4, 5]], ['S', 'H'], ['X'], 'hf\nmulti', 2, 0.01)

   
    #obj = InputGenerator('abcd.inp', 'hf')
    #obj.make_input()

    outobj = Output(sys.argv[1])
    print outobj.check_normalterm()


