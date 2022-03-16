import os
import numpy as np


def pad(s, l):
    while len(s)<l:
        s += ' '
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


def split_data(d):
    x = []
    for di in d:
        v = [float(t) for t in di.replace('\t','').replace('\n','').split(' ')]
        if x == []:
            for t in range(len(v)):
                x.append([])
        for t in range(len(x)):
            x[t].append(v[t])
    return x

    


def get_end_time_str():
    dirs = [name for name in os.listdir("./coarse/") if os.path.isdir("./coarse/" + name)]
    times_str = []
    times = []
    for dir in dirs:
        try:
            time = float(dir)
        except:
            time = -1
        if time != -1:
            times_str.append(dir)
            times.append(time)

    times, times_str = (list(t) for t in zip(*sorted(zip(times, times_str))))
    return times_str[-1]


# Import data
def import_data(end_time, line_tag, line, filepath='./coarse/postProcessing/'):

    # Read string data
    with open(filepath + line_tag + '/' + end_time + '/' + line + '_p_rho.xy') as f:
        p_rho = f.readlines()

    with open(filepath + line_tag + '/' + end_time + '/' + line + '_U.xy') as f:
        u = f.readlines()
    
    # Format to float
    p_rho = split_data(p_rho)
    u = split_data(u)
    x = u[0]
    y = u[1]
    z = u[2]
    u_x = u[3]
    u_z = u[5]
    p = p_rho[3]
    rho = p_rho[4]
    data = {'x':x, 'y':y, 'z':z, 'p':p, 'u_x':u_x, 'u_z':u_z, 'rho':rho}
    return data

def areaAverage(r, v):
    A = 0
    vA = 0
    for k in range(len(r)-1):
        vk = (v[k+1]+v[k])/2
        dA = np.pi*(r[k+1]**2 + r[k]**2)
        A += dA
        vA += vk * dA
    return vA/A



# Integrate variables to find ISP
def int_var(d, kind="pressure-momentum", normal="z", norm_fact=1.0, pressure_offset=0):
    r = d['x']
    z = d['z']
    p = d['p']
    u_x = d['u_x']
    u_z = d['u_z']
    rho = d['rho']
    F = 0
    m_dot = 0
    for t in range(len(r)-1):
        rhot = (rho[t]+rho[t+1])/2
        pt = (p[t]+p[t+1])/2
        uxt = (u_x[t]+u_x[t+1])/2
        uzt = (u_z[t]+u_z[t+1])/2
        if normal == "z":
            dA = np.abs(r[t+1]**2 - r[t]**2)*np.pi
            dm_dot = rhot * uzt * dA
        elif normal == "x":
            dA = np.abs(z[t+1] - z[t])*2*np.pi*(r[t]+r[t+1])/2
            dm_dot = rhot * uxt * dA
        dF = 0
        if "momentum" in kind:
            dF += uzt*dm_dot * norm_fact
        if "pressure" in kind and normal=="z":
            dF += (pt-pressure_offset)*dA * norm_fact
        F += dF
        m_dot += dm_dot
    return F, m_dot


def get_perf(verb=False):
    # Find largest time step folder name
    Patm = 101325
    last_time = get_end_time_str()

    # Find thrust and mass flow rates
    F_i0, m_dot_i0 = int_var( import_data(last_time, 'getMomentum', 'line0') , norm_fact=-1.0)
    F_i1, m_dot_i1 = int_var( import_data(last_time, 'getMomentum', 'line1') , normal="x" )
    F_i2, m_dot_i2 = int_var( import_data(last_time, 'getMomentum', 'line2') )

    d_core_in = import_data(last_time, 'getCoreMomentum', 'line0')
    d_core_out = import_data(last_time, 'getCoreMomentum', 'line1')
    F_core_in, m_dot_core_in =   int_var( d_core_in , norm_fact=-1.0, pressure_offset=Patm )
    F_core_out, m_dot_core_out = int_var( d_core_out, pressure_offset=Patm )

    F = F_i0 + F_i1 + F_i2

    F_core = F_core_out - F_core_in
    
    res = {
        'bypass_ratio':(m_dot_i2-m_dot_core_out)/m_dot_core_out,
        'thrust_ratio':F/F_core,
        'total_thrust':F,
        'core_thrust':F_core,
        'mass_flow_in':m_dot_core_in,
        'mass_flow_out':m_dot_core_out,
        'mass_flow_total':m_dot_i2,
        'mass_flow_bypass':m_dot_i2-m_dot_core_out,
        'mass_flow_error':np.abs(m_dot_core_in - m_dot_core_out)/m_dot_core_out,
    }
    return res




def editVar(variable, value, filename='./nozzleDict'):
    with open(filename,'rt') as file:
        lines = file.readlines()
    
    for k in range(len(lines)):
        l = lines[k].lstrip()
        if l[:len(variable)] == variable:
            lines[k] = variable + ' ' + str(value) + ';\n'
    
    with open(filename, 'w') as file:
        file.writelines(lines)


def run(var, values):


    step = 0
    for v in values:
        editVar(var, v)
        # Run simulation
        print(var + ' ' + str(v))
        os.system('./coarseRun > log.coarse')
        # Get results
        res = get_perf()

        this_strs = [var] + list(res.keys())
        this_values = ['{:.6e}'.format(k) for k in [v] + list(res.values())]
        with open('./results/batchResults.dat','a') as f:
            if step == 0:
                f.write(' '.join(this_strs) + '\n')
            f.write(' '.join(this_values) + '\n')
        step += 1


if __name__=='__main__':
    var = 'ej_radi'
    values = [0.4, 0.5, 0.6, 0.7, 0.8]
    run(var, values)

