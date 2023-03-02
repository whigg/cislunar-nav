import numpy as np

class Vector:
    '''
    time history of a 3D vector
     vec   - list of 3-element lists
     tspan - 2-element list of the form [start time; end time]
    '''
    def __init__(self, data, tspan):
        self._vec = np.array(data)
        self.start = tspan[0]
        self.end = tspan[1]

    def vec(self, time):
        return self._vec[int(time - self.start)]

    @property
    def x(self):
        ''' get time series of all x values '''
        return self._vec[:,0]

    @property
    def y(self):
        ''' get time series of all y values '''
        return self._vec[:,1]

    @property
    def z(self):
        ''' get time series of all z values '''
        return self._vec[:,2]

    def plot(self, ax):
        ''' plot values on given pyplot axes '''
        ax.plot3D(self.x, self.y, self.z)