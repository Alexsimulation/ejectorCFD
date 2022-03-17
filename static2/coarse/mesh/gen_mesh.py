import os
import numpy as np



def pad(s, l):
    iter = 0
    while (len(s)<l)&(iter < 1e6):
        s += ' '
        iter += 1
    return s


def print_dict(d):
    keys = d.keys()
    values = d.values()
    l = 0
    for k in keys:
        if len(k) > l:
            l = len(k)
    for k, v in zip(keys, values):
        print(pad(k,l), ':', v)


def editVars(variables, values, filename):
    with open(filename,'rt') as file:
        lines = file.readlines()
    
    for k in range(len(lines)):
        l = lines[k].lstrip()
        if l[:len(variable)] in variables:
            variable = l[:len(variable)]
            lines[k] = variable + '\t\t\t' + str(values[variable]) + '\n'
    
    with open(filename, 'w') as file:
        file.writelines(lines)


def readVars(filename):
    with open(filename,'rt') as file:
        lines = file.readlines()
    d = {}
    for k in range(len(lines)):
        l = lines[k].lstrip()
        l = l.replace('\t',' ')
        l = l.replace(';','')
        l = ' '.join(l.split())
        l = l.split(' ')
        if l[0] not in  ('//', ''):
            d[l[0]] = float(l[1])
    return d



def geomCalcs(d, filename='../system/geomCalcs'):

    d['ej_thickB'] = ( np.sqrt(d['ej_area_ratio']) - 1 ) * d['ej_radi'] / (1 + d['ej_thickTB'])
    d['ej_thickT'] = d['ej_thickB'] * d['ej_thickTB']

    v = {}
    ints = {}

    v['x3']	= (d['inlet_radi']*3 + d['hub_radi'])/(3+1)
    v['y3']	= v['x3'] * d['wedge_factor']
    v['y3m']	= -1.0 * v['y3']

    v['x5']	= d['hub_radi']
    v['y5']	= d['hub_radi'] * d['wedge_factor']
    v['y5m']	= -1.0 * v['y5']

    v['z0']	= d['hub_len']
    v['z1']	= d['hub_len'] + d['ej_back_len']
    v['z2']	= v['z0'] * (1 + d['backL_totalL'])
    v['z5']	= d['hub_len'] - d['hub_end']

    v['x7']	= (d['hub_radi'] + d['noz_radi']*d['hub_end_fact'])/(d['hub_end_fact'] + 1)
    v['z7']	= (v['z5'] + v['z0'])/(1 + 1)
    v['y7']	=  v['y3'] / v['x3'] * v['x7']
    v['y7m']	= -1.0*v['y7']

    v['x9']	= d['noz_radi']
    v['y9']	=  v['y3'] / v['x3'] * v['x9']
    v['y9m']	= -1.0*v['y9']

    v['x21']	= d['ej_radi']
    v['z21']	= v['z0'] + d['ej_z_delta']
    v['y21']	=  v['y3'] / v['x3'] * v['x21']
    v['y21m']	= -1.0*v['y21']

    v['x19']	= d['ej_radi'] + d['ej_thickB']
    v['y19']	=  v['y3'] / v['x3'] * v['x19']
    v['y19m']	= -1.0*v['y19']
    v['z19']	= v['z21'] - d['ej_front_len']

    v['x23']	= d['ej_radi'] + d['ej_thickB'] + d['ej_thickT']
    v['y23']	=  v['y3'] /v['x3']*v['x23']
    v['y23m']	= -1.0*v['y23']

    v['x29']	= d['hub_radi'] * (1 + d['outR_ej_R']) 
    v['y29']	=  v['y3'] / v['x3'] * v['x29']
    v['y29m']	= -1.0 * v['y29']

    v['x15']	= v['x5'] + d['ej_radi'] - d['noz_radi']
    v['y15']	=  v['y3'] / v['x3'] * v['x15']
    v['y15m']	= -1.0 * v['y15']


    v['z41']	= d['inlet_len']
    v['z43']	= -1.0*d['hub_len']*d['frontL_totalL']
    v['x44']	= d['inlet_radi']
    v['y44']	=  v['y3'] / v['x3'] * v['x44']
    v['y44m']	= -1.0 * v['y44']

    v['z15']	= -1.0 * d['inlet_len'] / 2



    # Spline calculations

    v['z_3_5_bs']	=  v['z5'] / 3
    v['z_3_5_bs2']	=  v['z5'] / 2
    v['z_5_7_arc']	= v['z0'] + (v['z0'] - v['z7'])*(v['x9'] - v['x5'])/(v['x7'] - v['x9'])

    v['z21arc']	= (v['z19'] + v['z21'])/(1+1)
    v['x21arc']	= (v['x19'] + v['x21'])/(1+1)
    v['y21arc']	=  v['y3'] / v['x3'] * v['x21arc']
    v['y21arcm']	= -1.0*v['y21arc']

    v['z27arc']	= (v['z19'] + v['z21'])/(1+1)
    v['x27arc']	= (v['x19'] + v['x23'])/(1+1)
    v['y27arc']	=  v['y3'] / v['x3'] * v['x27arc']
    v['y27arcm']	= -1.0*v['y27arc']

    v['z23arc'] = (v['z21']*4 + v['z1'])/(4+1)

    v['x3arc']		= (v['x3']*2 + v['x5'])/(2+1)
    v['y3arc']		=  v['x3arc'] * d['wedge_factor']
    v['y3arcm']	= -1.0*v['y3arc']


    # Number of elements
    ints['base_x_3_15']	= 65
    ints['nelm_x_3_15']	= round( (v['x21'] - v['x9'])/0.5 * ints['base_x_3_15'] )
    

    ints['base_z_0_1']		= 95
    ints['nelm_z_0_1']		= round( (v['z1'] - v['z0'])/2.0 * ints['base_z_0_1'] )


    with open(filename, 'w') as f:
        f.write('\n#include "../../nozzleDict"\n\n')
        for var in v:
            s = var + '\t' + '{:.16e}'.format(v[var]) + ';\n'
            f.write( s )
        for var in ints:
            s = var + '\t' + '{:d}'.format(ints[var]) + ';\n'
            f.write( s )
        f.write('\n')



def main():
    d = readVars('../../nozzleDict')
    geomCalcs(d)



if __name__=='__main__':
    main()



