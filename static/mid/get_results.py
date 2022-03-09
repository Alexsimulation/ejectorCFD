import os
import numpy as np


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
    dirs = [name for name in os.listdir(".") if os.path.isdir(name)]
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
def import_data(end_time, line_tag, line):

    # Read string data
    with open('./postProcessing/' + line_tag + '/' + end_time + '/' + line + '_p_rho.xy') as f:
        p_rho = f.readlines()

    with open('./postProcessing/' + line_tag + '/' + end_time + '/' + line + '_U.xy') as f:
        u = f.readlines()
    
    # Format to float
    p_rho = split_data(p_rho)
    u = split_data(u)
    r = u[0]
    u_z = u[3]
    p = p_rho[1]
    rho = p_rho[2]
    data = {'r':r, 'p':p, 'u':u_z, 'rho':rho}
    return data


# Integrate variables to find ISP
def int_var(d, kind="pressure-momentum"):
    r = d['r']
    p = d['p']
    u = d['u']
    rho = d['rho']
    F = 0
    m_dot = 0
    for t in range(len(r)-1):
        dA = np.abs(r[t+1]**2 - r[t]**2)*np.pi
        rhot = (rho[t]+rho[t+1])/2
        pt = (p[t]+p[t+1])/2
        ut = (u[t]+u[t+1])/2
        dm_dot = rhot * ut * dA
        dF = 0
        if "momentum" in kind:
            dF += ut*dm_dot
        if "pressure" in kind:
            dF += pt*dA
        F += dF
        m_dot += dm_dot
    return F, m_dot


def get_perf(verb=False):
    # Find largest time step folder name
    last_time = get_end_time_str()

    # Find thrust and mass flow rates
    F_e, m_dot_e = int_var( import_data(last_time, 'getExhaustMomentum', 'line') )
    F_j, m_dot_j = int_var( import_data(last_time, 'getJetMomentum', 'line') )
    F_i0, m_dot_i0 = int_var( import_data(last_time, 'getInletMomentum', 'line0') , kind="pressure")
    F_i1, m_dot_i1 = int_var( import_data(last_time, 'getInletMomentum', 'line1') )
    F_i2, m_dot_i2 = int_var( import_data(last_time, 'getInletMomentum', 'line2') )
    F_i = F_i0 + F_i1 + F_i2

    entrainement_ratio = (m_dot_e - m_dot_j)/m_dot_j
    thrust_ratio = (F_e-F_i)/(F_j)

    if verb:
        print('entrainement_ratio:','{:.8e}'.format(entrainement_ratio),' | thrust_ratio:','{:.8e}'.format(thrust_ratio))
    
    res = {
        'm_dot':entrainement_ratio,
        'thrust':thrust_ratio,
    }
    return res




if __name__=='__main__':
    get_perf(verb=True)
