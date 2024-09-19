

import os
from math import *
import numpy as np
import datetime as dt
from scipy.integrate import odeint
from load_gmat import *    # GMAT API, need to download


mu = 0.00490e6 # km^3/s^2


def dcm_mat_gen(ang, ax):
    c = cos(ang)
    s = sin(ang)
    
    if ax == 1:
        dcm_mat = np.array([[1, 0, 0],
                            [0, c, s],
                            [0, -s, c]])
    elif ax == 2:
        dcm_mat = np.array([[c, 0, -s],
                            [0, 1, 0],
                            [s, 0, c]])
    elif ax == 3:
        dcm_mat = np.array([[c, s, 0],
                            [-s, c, 0],
                            [0, 0, 1]])
    
    return dcm_mat


def orb2st(a, ecc, inc, raan, arg_p, f):
    
    raan = radians(raan)
    arg_p = radians(arg_p)
    f = radians(f)

    p = a * (1 - ecc**2)

    rmag = p / (1 + (ecc * cos(f)))
    r_fdot = sqrt(mu/p) * (1 + (ecc * cos(f)))
    rmag_dot = sqrt(mu/p) * ecc * sin(f)

    r_pqw = np.array([[rmag*cos(f)], [rmag*sin(f)], [0]])
    v_pqw = np.array([[(rmag_dot*cos(f)) - (r_fdot*sin(f))], [((rmag_dot*sin(f)) + (r_fdot*cos(f)))], [0]])

    dcm_mat = dcm_mat_gen(-raan, 3) @ dcm_mat_gen(-inc, 1) @ dcm_mat_gen(-arg_p, 3)
    r = dcm_mat * r_pqw
    v = dcm_mat * v_pqw

    return r, v


def twoBodyProb(y, t):

    const = -mu / (np.norm(y[:3])**3)

    rx = y[0]
    ry = y[1]
    rz = y[2]
    vx = y[3]
    vy = y[4]
    vz = y[5]

    dydt = np.zeros(6)

    dydt[0] = vx
    dydt[1] = vy
    dydt[2] = vz

    dydt[3] = const * rx
    dydt[4] = const * ry
    dydt[5] = const * rz

    return dydt


def propagate2Bod(r0, v0, ts_2b):
    
    y0 = np.concatenate(r0, v0)

    sv_2b = odeint(twoBodyProb, y0, ts_2b)

    return sv_2b


def getDM(time_to_return, period):
    frac = time_to_return / period
    part = frac - np.floor(frac)

    if part < 0.5:
        dm = 1
    else:
        dm = -1
    
    return dm


def t2z(t, r1_vec, r2_vec, DM, num_revs, mu, tol):
    
    r1 = np.linalg.norm(r1_vec)
    r2 = np.linalg.norm(r2_vec)
    
    df = acos(np.dot(r1_vec, r2_vec) / (r1 * r2))
    if DM == -1:
        df = (2.*pi) - df    
    
    A = DM * sqrt(r1 * r2 * (1 + cos(df)))
    
    if DM == 1:
        z_n = (2 * num_revs * pi)**2 + 1
    elif DM == -1:
        if num_revs == 0:
            z_n = 1
        else:
            z_n = (2 * (num_revs+1) * pi)**2 - 1
    
    err = 1e10
    while err > tol:
        alph = 0.01
        
        t_n = z2t(z_n, r1_vec, r2_vec, DM, mu)
        
        c, s = z2sc(z_n)
        cp, sp = z2spcp(z_n)

        y = r1 + r2 - (A * (1 - (z_n * s))) / sqrt(c)
        x = sqrt(y / c)
        
        dtdz = (1 / sqrt(mu)) * ((x**3) * (sp - (3*s*cp) / (2*c)) + (A/8) * (((3*s*sqrt(y))/c) + (A/x)))
        
        z_n = z_n + (alph * ((t - t_n) / dtdz))
        
        err = abs(t - t_n)    
    
    z = z_n
    return z


def z2t(z, r1_vec, r2_vec, DM, mu):
    
    r1 = np.linalg.norm(r1_vec)
    r2 = np.linalg.norm(r2_vec)
    
    df = acos(np.dot(r1_vec, r2_vec) / (r1*r2))
    if DM == -1:
        df = (2*pi) - df
    
    c, s = z2sc(z)
    
    A = DM * sqrt(r1 * r2 * (1 + cos(df)))
    
    y = r1 + r2 - (A * (1 - (z*s)))/sqrt(c)
    
    x = sqrt(y/c)
    
    t = (1/sqrt(mu)) * ((x**3)*s + (A*sqrt(y)))
    return t


def z2sc(z):
    if z > 1e-6:
        c = (1 - cos(sqrt(z))) / z
        s = (sqrt(z) - sin(sqrt(z)))/sqrt(z**3)
    elif z < -1e-6:
        c = (1 - cosh(sqrt(-z))) / z
        s = (sinh(sqrt(-z)) - sqrt(-z)) / sqrt((-z)**3)
    else:
        c = (1/2) - (z/factorial(4)) + ((z**2)/factorial(6)) - ((z**3)/factorial(8))
        s = (1/6) - (z/factorial(5)) + ((z**2)/factorial(7)) - ((z**3)/factorial(9))

    return c, s


def z2spcp(z):
    
    c, s = z2sc(z)
    
    if abs(z) > 1e-6:
        cp = (1/(2*z)) * (1 - (z*s) - (2*c))
        sp = (1/(2*z)) * (c - (3*s))
    else:
        cp = (-1/factorial(4)) + ((2*z)/factorial(6)) - ((3*(z**2))/factorial(8)) + ((4*(z**3))/factorial(10))
        sp = (-1/factorial(5)) + ((2*z)/factorial(7)) - ((3*(z**2))/factorial(9)) + ((4*(z**3))/factorial(11))

    return cp, sp


def z2y(z, r1_vec, r2_vec, DM):
    
    r1 = np.linalg.norm(r1_vec)
    r2 = np.linalg.norm(r2_vec)
    
    df = acos(np.dot(r1_vec, r2_vec) / (r1*r2))
    if DM == -1:
        df = (2*pi) - df
    
    c, s = z2sc(z)
    
    A = DM * sqrt(r1 * r2 * (1 + cos(df)))
    
    y = r1 + r2 - (A * (1 - (z*s))) / sqrt(c)

    return y


def FandG_z(z, r1_vec, r2_vec, DM, mu):
    
    r1 = np.linalg.norm(r1_vec)
    r2 = np.linalg.norm(r2_vec)
    
    y = z2y(z, r1_vec, r2_vec, DM)
    
    df = acos(np.dot(r1_vec, r2_vec) / (r1*r2))
    if DM == -1:
        df = (2*pi) - df
    
    A = DM * sqrt(r1 * r2 * (1 + cos(df)))
    
    F = 1 - (y/r1)
    G = A*sqrt(y/mu)
    Gdot = 1 - (y/r2)
    Fdot = ((F*Gdot)-1)/G

    return F, G, Fdot, Gdot


def lambert_z(r1_vec, r2_vec, t, dm, num_revs, mu):

    tol = 1e-6
    z = t2z(t, r1_vec, r2_vec, dm, num_revs, mu, tol)
    
    F, G, Fdot, Gdot = FandG_z(z, r1_vec, r2_vec, dm, mu)
    
    v1_vec = (r2_vec - (F * r1_vec)) / G
    v2_vec = ((Gdot * r2_vec) - r1_vec) / G

    return v1_vec, v2_vec



def generateScript(a, ecc, inc, raan, arg_p, f, epoch, step_size, t_end, fname, ind):
    # Name of GMAT variable corresponding to each orbital element
    lines = ['GMAT Orbiter.SMA', 'GMAT Orbiter.ECC', 'GMAT Orbiter.INC',
            'GMAT Orbiter.RAAN', 'GMAT Orbiter.AOP', 'GMAT Orbiter.TA']
    
    # function to generate new file path + name for GMAT report file
    report_path = "GMAT ReportFile.Filename = '" + os.getcwd() + "\gmat_reports\\" + fname + str(ind) + "_MoonOrb.txt'\n"
    
    # Store script in array
    temp = []
    # open template file and store all the lines in program memory
    with open('station_keeping/gmat_scripts/temp_MoonOrb.script', 'r') as template:
        for line in template:
            temp.append(line)

    # path + name of new script to make
    name = 'station_keeping/gmat_scripts/' + fname + str(ind) + '_MoonOrb.script'

    with open(name, 'w') as file:
        # replace each line that starts with one of our variables with new value
        for line in temp:
            if line.startswith(lines[0]):
                file.write(lines[0] + ' = ' + str(a) + ';\n')
            elif line.startswith(lines[1]):
                file.write(lines[1] + ' = ' + str(ecc) + ';\n')
            elif line.startswith(lines[2]):
                file.write(lines[2] + ' = ' + str(inc) + ';\n')
            elif line.startswith(lines[3]):
                file.write(lines[3] + ' = ' + str(raan) + ';\n')
            elif line.startswith(lines[4]):
                file.write(lines[4] + ' = ' + str(arg_p) + ';\n')
            elif line.startswith(lines[5]):
                file.write(lines[5] + ' = ' + str(f) + ';\n')
            elif line.startswith('GMAT ReportFile.Filename'):
                file.write(report_path)
            elif line.startswith('GMAT DefaultProp.InitialStepSize'):
                file.write('GMAT DefaultProp.InitialStepSize = ' + str(step_size) + ';\n')
            elif line.startswith('GMAT Orbiter.Epoch'):
                file.write('GMAT Orbiter.Epoch = ' + epoch + ';\n')
            elif line.startswith('Propagate DefaultProp(Orbiter)'):
                file.write('Propagate DefaultProp(Orbiter) {Orbiter.ElapsedDays = ' + str(t_end / (3600 * 24)) + '};\n')
            else:   # if not a line we want to change, just write it to new file
                file.write(line)
    
    
def runScripts(fname, ind):
    '''
    Runs all GMAT scripts in a given directory

    NOTE: Just tries to run anything that is a file, so make sure there are
    only GMAT .script files in the directory (folders are OK)
    '''

    name = 'station_keeping/gmat_scripts/' + fname + str(ind) + '_MoonOrb.script'

    gmat.LoadScript(name)
    gmat.RunScript()


def readGMATReport(fname, ind, epoch_dt):
    with open('station_keeping\gmat_reports\\' + fname + str(ind) + '_MoonOrb.txt', 'r') as file:
        allLines = file.readlines()
        numLines = len(allLines) - 1

        ts_gmat = np.zeros(numLines)
        sv_gmat = np.zeros((numLines, 6))
        
        line = file.readline()

        lineNum = -1
        while True:
            line = file.readline()
            lineNum += 1
            if not line:
                break

            data = line.split()

            sv_gmat[lineNum, 0] = float(data[4])
            sv_gmat[lineNum, 1] = float(data[5])
            sv_gmat[lineNum, 2] = float(data[6])
            sv_gmat[lineNum, 3] = float(data[7])
            sv_gmat[lineNum, 4] = float(data[8])
            sv_gmat[lineNum, 5] = float(data[9])

            t_str = data[0] + ' ' + data[1] + ' ' + data[2] + ' ' + data[3] + '000'
            t_i = dt.datetime.strptime(t_str, '%d %b %Y %H:%M:%S.%f')
            t_diff = t_i - epoch_dt
            ts_gmat[lineNum] = t_diff.total_seconds()
    
    return ts_gmat, sv_gmat


if __name__ == '__main__':
    print('Running')

    # User inputs
    a = 3000 # km
    ecc = 0.1
    inc = 40 # deg
    raan = 80 # deg
    arg_p = 50 # deg
    f0 = 0 # deg

    t_step = 10 # s
    t_end = 86400 # s

    epoch = '26 May 2024 12:16:34.044'

    pos_error_tol = 100 # km
    time_to_return = 2 * pi * sqrt((a**3)/mu) * 0.25 # s 
    # time_to_return is currently set at a quarter rev as this may make the calculations close better,
    # though further testing could be done to better optimize station-keeping, as this code is ultimately a crude estimation
    if time_to_return < t_step:
        time_to_return = t_step

    fname = 'skTest'


    # Derived inputs
    period = 2 * pi * sqrt((a**3)/mu) # s
    epoch_dt = dt.datetime.strptime(epoch, '%d %b %Y %H:%M:%S.%f')



    # Get propagated 2-body states (reference trajectory)
    ts_2b = np.arange(0, t_end, t_step)
    r0, v0 = orb2st(a, ecc, inc, raan, arg_p, f0)
    sv_2b = propagate2Bod(r0, v0, ts_2b)

    # Generate and run the GMAT script
    ind = 0
    generateScript(a, ecc, inc, raan, arg_p, f0, epoch, t_step, t_end, fname, ind)
    runScripts(fname, ind)

    # Read in GMAT data
    ts_gmat, sv_gmat = readGMATReport(fname, ind, epoch_dt)


    # Loop through ephemeris and check position error
    i = 0
    dv_total = 0
    while i < (t_end/t_step):
        # Find current time
        curr_time = ts_2b[i]

        # Get current position vectors for reference and GMAT trajectories
        pos_2b = sv_2b[i, :3]
        pos_gmat = sv_gmat[i, :3]

        # Find position error
        pos_error = np.linalg.norm(sv_gmat - sv_2b)

        # If position error is above a specified threshold, maneuver
        if pos_error > pos_error_tol:
            # Round the "time to return" to the nearest timestamp
            ttr_rounded = (np.round(time_to_return/t_step) * t_step)
            t_return = curr_time + ttr_rounded
            # Cap the return time to the end of the simulation
            if t_return > t_end:
                t_return = t_end

            # Find the reference trajectory position at the time of return
            i_return = t_return / t_step
            pos_2b_return = sv_2b[i_return, :3]

            # Calculate starting velocity before first maneuver and final velocity after second maneuver
            v1_gmat = sv_gmat[i, 3:]
            v2_2b = sv_2b[i, 3:]

            # Solve Lambert's problem to find transfer starting and ending velocity vectors
            dm = getDM(ttr_rounded, period)
            num_revs = np.floor(ttr_rounded/period)
            v1_l, v2_l = lambert_z(pos_gmat, pos_2b_return, ttr_rounded, dm, num_revs, mu)
            
            # Check that new velocity is still going in the same direction; if not, fix
            if np.dot(v1_gmat, v1_l) < 0:
                v1_l, v2_l = lambert_z(pos_gmat, pos_2b_return, ttr_rounded, -dm, num_revs, mu)

            # Calculate delta-Vs for both maneuvers
            dv1_vec = v1_l - v1_gmat
            dv1 = np.linalg.norm(dv1_vec)

            dv2_vec = v2_2b - v2_l
            dv2 = np.linalg.norm(dv2_vec)

            # Add to total delta-V
            dv_total = dv_total + dv1 + dv2

            # Update loop index
            i = i_return


            if t_return < t_end:
                # Calculate new true anomaly
                rmag_return = np.linalg.norm(pos_2b_return)
                vmag_return = np.linalg.norm(v2_2b)
                
                eng = (vmag_return**2)/2 - (mu/rmag_return)
                hvec = np.cross(pos_2b_return, v2_2b)
                evec = (1/mu) * (np.cross(v2_2b, hvec) - ((mu * pos_2b_return)/rmag_return))
    
                f_new = acos(np.dot(evec, pos_2b_return) / (ecc * rmag_return))
                if np.dot(pos_2b_return, v2_2b) < 0:
                    f_new = (2*pi) - f_new
                
                # Generate and run the new GMAT script
                ind += 1
                new_epoch_dt = epoch_dt + dt.timedelta(seconds=t_return)
                generateScript(a, ecc, inc, raan, arg_p, f_new, new_epoch_dt.strftime('%d %b %Y %H:%M:%S.%f'), t_step, t_end, fname, ind)
                runScripts(fname, ind)

                # Read in new GMAT data and update times
                ts_gmat, sv_gmat = readGMATReport(fname, ind, new_epoch_dt)
                ts_gmat = ts_gmat + t_return

                # Pad new results to keep indexing consistent
                ts_gmat = np.concatenate((np.zeros(i), ts_gmat))
                sv_gmat = np.concatenate((np.zeros((i,6)), sv_gmat), axis=0)

        else:
            i += 1


    print('Total delta-V: ' + dv_total)         

