from Vector import Vector
import re
import numpy as np
import matplotlib.pyplot as plt

def parseGmatData(filename, gmatReport=False):
    '''
    Parse position data and sort into Satellite objects
     filename - full relative path to file
    '''
    
    head = True
    sats = []
    with open(filename) as data:
        if gmatReport:
            for line in data:
                if head:
                    head = False
                    tspan = [0, 0]
                    n = len(line.split())
                    posarray = np.zeros([1,n])
                else:
                    row = np.array([[float(x.strip()) for x in line.split()]])
                    posarray = np.concatenate((posarray, row), axis=0)
                    tspan[1] += 1

            tspan[1] -= 1
            for i in range(0,int(n/3)):
                pos = posarray[1:,3*i:3*i+3].tolist()
                sats.append(Vector(pos, tspan))
            if n % 3 != 0:  # throw excess data in if there's more
                sats.append(posarray[1:,3*int(n/3):])
        else:
            for line in data:
                if head:
                    # parse first line for time stamps
                    head = False
                    time = [int(x) for x in line.strip().split()]
                    tspan = [time[0], time[-1]]
                else:
                    # capture anything between "[...]" excl. brackets and quotes
                    vecs = re.findall(r'"\[(.*?)\]"', line)
                    pos = []
                    # convert list of position strings to list of 3-element positions
                    [pos.append([float(x) for x in vec.split(', ')]) for vec in vecs]
                    # create satellite from data
                    sats.append(Vector(pos, tspan))

    return sats

def findVisibleSats(pos, sats, time, elev=5):
    '''
    Finds satellites visible to lunar surface user

     pos  - satellite position vector (km)
     sats - list of satellites to consider
     time - point in time to evaluate
     elev - minimum elevation of satellites (deg), aka mask angle
    '''
    elev = elev * np.pi/180     # elevation angle (5 deg)
    r = 1737.4                  # km, volumetric mean radius of moon
    zcone = lambda x,y: -r - np.tan(elev) * np.sqrt(x**2 + y**2)

    ## Compute rotation matrix
    rv = np.array([0, 0, -r])
    # Euler axis and angle for rotation
    if np.linalg.norm(np.cross(pos,rv)) != 0.:  # make sure vectors aren't already aligned
        ax = np.cross(pos, rv) / np.linalg.norm(np.cross(pos, rv))
        ang = np.arccos(np.dot(pos, rv) / (np.linalg.norm(pos) * np.linalg.norm(rv)))
    else:   # set to arbitrary axis and zero rotation
        ax = np.array([1,0,0])
        ang = 0.
    U = np.array([[0, -ax[2], ax[1]],[ax[2], 0, -ax[0]],[-ax[1], ax[0], 0]])
    R = np.cos(ang) * np.eye(3) + np.sin(ang) * U + (1 - np.cos(ang)) * np.outer(ax, ax)
    # print(R @ pos)

    visible = []
    for sat in sats:
        sx,sy,sz = sat.vec(time)
        loc = R @ np.array([sx, sy, sz])    # rotate satellite position
        if loc[2] < zcone(loc[0],loc[1]):
            visible.append(sat)

    return visible, R

def computeDOP(user, sats, time):
    '''
    Computes the Dilution of Precision (DOP) matrix H = inv(G'G)
     user - inertial user position as numpy array
     sats - satellites under consideration
     time - time to evaluate DOP at
    '''
    num = len(sats)
    G = np.ones((num,4))
    # Get latitude and longitude for ECI->ENU rotation
    lat = np.pi/2 - np.arctan2(user[2], np.sqrt(user[1]**2 + user[2]**2))
    lon = np.pi/2 + np.arctan2(user[1], user[0])

    R1 = np.array([[1,0,0],[0,np.cos(lat),np.sin(lat)],[0,-np.sin(lat),np.cos(lat)]])
    R2 = np.array([[np.cos(lon),np.sin(lon),0],[-np.sin(lon),np.cos(lon),0],[0,0,1]])
    R = np.concatenate((R1@R2, np.zeros((3,1))), axis=1)  
    R = np.concatenate((R, np.array([[0,0,0,1]])), axis=0)  # ECI->ENU rotation matrix

    for i in range(0,num):
        ik = (sats[i].vec(time) - user)
        G[i,0:3] = -ik / np.linalg.norm(ik)
    
    G = G@np.transpose(R)
    try:
        return np.linalg.inv(np.transpose(G)@G)
    except np.linalg.LinAlgError: 
        return float("inf")*np.ones((4,4))


def set_axes_equal(ax):
    '''
    Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

def goodness_of_fit(xdata, ydata, xlabel, ylabel):
    '''
    Scatter data of y vs. x, where x and y and 3D vectors, then the r^2
    value to y = x
     xdata  - (3,n) matrix, truth data
     ydata  - (3,n) matrix, estimator data
     xlabel - str, what to call this in plots
     ylabel - str, what to call this in plots
    '''
    fig = plt.figure()
    ax = [None] * 3
    ax[0] = plt.subplot(131)
    ax[1] = plt.subplot(132)
    ax[2] = plt.subplot(133)
    vavg = 0

    for i in range(0,3):
        # Compute coefficient of determination
        x = xdata[i]
        y = ydata[i]
        mean = np.mean(y)
        sse = np.dot(y - x, y - x)
        sst = np.dot(y - mean, y - mean)
        vavg += sse/len(y)

        # Plot data
        ax[i].scatter(x, y, s=5)
        ax[i].set_title("$\sigma^2$ = {:.4f} $cm^2/s^4$".format(sse/len(y)*1e10))
        axis= ['x', 'y', 'z']
        #ax[i].set_title("Force, " + axis[i])
        ax[i].grid()

    ax[0].set_ylabel(ylabel)
    ax[0].set_xlabel(xlabel+ ", x")
    ax[1].set_xlabel(xlabel+ ", y")
    ax[2].set_xlabel(xlabel+ ", z")
    print(vavg)

                    

if __name__ == "__main__":
    sats = parseGmatData("data/const_eph/0_MoonOrb.txt", gmatReport=True)
    t = 100;        # evaluation time
    r = 1737.4;     # moon radius, km
    pos = np.array([np.sin(np.pi/12)*r, 0, -np.cos(np.pi/12)*r])
    visible, R = findVisibleSats(pos, sats, t)
    H = computeDOP(pos, visible, t)
    print(f'PDOP: {np.sqrt(H[0,0]**2 + H[1,1]**2 + H[2,2]**2)}')

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_aspect('auto')

    # plot cone
    elev = 5 * np.pi/180;   # elevation angle (5 deg)
    r = 1737.4              # km, volumetric mean radius of moon
    x,y = np.mgrid[-23000:23000:1000j, -23000:23000:1000j]    
    z = -r - np.tan(elev) * np.sqrt(x**2 + y**2)
    ax.plot_surface(x, y, z)

    # plot moon
    u, v = np.mgrid[0:2 * np.pi:30j, 0:np.pi:20j]
    x = r * np.cos(u) * np.sin(v)
    y = r * np.sin(u) * np.sin(v)
    z = r * np.cos(v)
    ax.plot_surface(x, y, z)

    # # plot satellites
    # for sat in sats:
    #     sat.plot(ax)
    for sat in visible:
        x,y,z = sat.vec(t)
        ax.scatter(x,y,z)
        #loc = R @ np.array([x, y, z])    # rotate satellite position
        #ax.scatter(loc[0], loc[1], loc[2])
    set_axes_equal(ax)
    plt.show()

    pass