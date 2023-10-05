# library imports
import spiceypy as spice

class LunarDynamics:

    def __init__(self):
        '''
        Initialize necessary SPICE kernels for dynamics function
        '''
        spice.furnsh(('data/kernels/de430.bsp','data/kernels/naif0012.tls'))


    def __enter__(self):
        ''' Support for with statements '''
        return self

    def __exit__(self, ex_type, ex_val, ex_traceback):
        ''' Catch exit errors '''
        spice.kclear()              # Unload kernels

        if ex_type is not None:
            print("\nExecution type:", ex_type)
            print("\nExecution value:", ex_val)
            print("\nTraceback:", ex_traceback)

    def orbdyn(self, t, x):
        '''
        compute state vector derivatives for a lunar orbiter
        '''
        merc  = spice.spkezr('MERCURY', t, 'J2000', 'NONE', 'MOON')[0][0:3]
        venus = spice.spkezr('VENUS', t, 'J2000', 'NONE', 'MOON')[0][0:3]
        earth = spice.spkezr('EARTH', t, 'J2000', 'NONE', 'MOON')[0][0:3]
        moon  = spice.spkezr('MOON' , t, 'J2000', 'NONE', 'MOON')[0][0:3]
        jupi  = spice.spkezr('JUPITER BARYCENTER', t, 'J2000', 'NONE', 'MOON')[0][0:3]
        saturn= spice.spkezr('SATURN BARYCENTER', t, 'J2000', 'NONE', 'MOON')[0][0:3]
        uranus= spice.spkezr('URANUS BARYCENTER', t, 'J2000', 'NONE', 'MOON')[0][0:3]
        nept  = spice.spkezr('NEPTUNE BARYCENTER', t, 'J2000', 'NONE', 'MOON')[0][0:3]
        sun   = spice.spkezr('SUN'  , t, 'J2000', 'NONE', 'MOON')[0][0:3]