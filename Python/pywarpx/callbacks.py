# Copyright 2017-2022 The WarpX Community
#
# This file is part of WarpX.
#
# Authors: David Grote, Roelof Groenewald
#
# License: BSD-3-Clause-LBNL

"""call back operations
=====================

These are the functions which allow installing user created functions so that
they are called at various places along the time step.

The following three functions allow the user to install, uninstall and verify
the different call back types.
 - installcallback: Installs a function to be called at that specified time
 - uninstallcallback: Uninstalls the function (so it won't be called anymore)
 - isinstalled: Checks if the function is installed

These functions all take a callback location name (string) and function or
instance method as an argument. Note that if an instance method is used, an
extra reference to the method's object is saved.

The install can be done using a decorator, which has the prefix "callfrom". See
example below.

Functions can be called at the following times:
 - beforeInitEsolve: before the initial solve for the E fields (i.e. before the PIC loop starts)
 - afterinit: immediately after the init is complete
 - beforeEsolve: before the solve for E fields
 - poissonsolver: In place of the computePhi call but only in an electrostatic simulation
 - afterEsolve: after the solve for E fields
 - beforedeposition: before the particle deposition (for charge and/or current)
 - afterdeposition: after particle deposition (for charge and/or current)
 - beforestep: before the time step
 - afterstep: after the time step
 - afterdiagnostics: after diagnostic output
 - oncheckpointsignal: on a checkpoint signal
 - onbreaksignal: on a break signal. These callbacks will be the last ones executed before the simulation ends.
 - particlescraper: just after the particle boundary conditions are applied
                    but before lost particles are processed
 - particleloader: at the time that the standard particle loader is called
 - particleinjection: called when particle injection happens, after the position
                      advance and before deposition is called, allowing a user
                      defined particle distribution to be injected each time step
 - appliedfields: allows directly specifying any fields to be applied to the particles
                  during the advance

To use a decorator, the syntax is as follows. This will install the function
``myplots`` to be called after each step.

@callfromafterstep
def myplots():
  ppzx()

This is equivalent to the following:

def myplots():
  ppzx()
installcallback('afterstep', myplots)

"""
import copy
import sys
import time
import types

import numpy

from ._libwarpx import libwarpx


class CallbackFunctions(object):
    """
    Class to handle call back function lists.

    Note that for functions passed in that are methods of a class instance,
    a full reference of the instance is saved. This extra reference means
    that the object will not actually be deleted if the user deletes the
    original reference. This is good since the user does not need to keep
    the reference to the object (for example it can be created using a local
    variable in a function). It may be bad if the user thinks an object was
    deleted, but it actually isn't since it had (unkown to the user)
    installed a method in one of the call back lists.
    """

    def __init__(self,name=None,lcallonce=0,singlefunconly=False):
        self.funcs = []
        self.time = 0.
        self.timers = {}
        self.name = name
        self.lcallonce = lcallonce
        self.singlefunconly = singlefunconly

    def __call__(self,*args,**kw):
        """Call all of the functions in the list"""
        tt = self.callfuncsinlist(*args,**kw)
        self.time = self.time + tt
        if self.lcallonce: self.funcs = []

    def clearlist(self):
        """Unregister/clear out all registered C callbacks"""
        self.funcs = []
        libwarpx.libwarpx_so.remove_python_callback(self.name)

    def __bool__(self):
        """Returns True if functions are installed, otherwise False"""
        return self.hasfuncsinstalled()

    def __len__(self):
        """Returns number of functions installed"""
        return len(self.funcs)

    def hasfuncsinstalled(self):
        """Checks if there are any functions installed"""
        return len(self.funcs) > 0

    def _getmethodobject(self,func):
        """For call backs that are methods, returns the method's instance"""
        return func[0]

    def callbackfunclist(self):
        """Generator returning callable functions from the list"""
        funclistcopy = copy.copy(self.funcs)
        for f in funclistcopy:
            if isinstance(f,list):
                object = self._getmethodobject(f)
                if object is None:
                    self.funcs.remove(f)
                    continue
                result = getattr(object,f[1])
            elif isinstance(f,str):
                import __main__
                if f in __main__.__dict__:
                    result = __main__.__dict__[f]
                    # --- If the function with the name is found, then replace the
                    # --- name in the list with the function.
                    self.funcs[self.funcs.index(f)] = result
                else:
                    continue
            else:
                result = f
            if not callable(result):
                print("\n\nWarning: a call back was found that is not callable.")
                if self.name is not None:
                    print("For %s"%self.name)
                print("Only callable objects can be installed.")
                print("It is possible that the callable's name has been overwritten")
                print("by something not callable. This can happen during restart")
                print("if a function name had later been used as a variable name.")
                print(self.name)
                if isinstance(f,str):
                    print("The name of the call back is %s"%f)
                print("\n\n")
                continue
            yield result

    def installfuncinlist(self,f):
        """Check if the specified function is installed"""
        if self.singlefunconly and self.hasfuncsinstalled():
            raise RuntimeError(
                f"Only one function can be installed for callback {self.name}."
            )

        if len(self.funcs) == 0:
            # If this is the first function installed, set the callback in the C++
            # to call this class instance.
            libwarpx.libwarpx_so.add_python_callback(self.name, self)
        if isinstance(f,types.MethodType):
            # --- If the function is a method of a class instance, then save a full
            # --- reference to that instance and the method name.
            finstance = f.__self__
            fname = f.__name__
            self.funcs.append([finstance,fname])
        elif callable(f):
            # --- If a function had already been installed by name, then skip the install.
            # --- This is problematic, since no warning message is given, but it is unlikely
            # --- to arise under normal circumstances.
            # --- The purpose of this check is to avoid redundant installation of functions
            # --- during a restore from a dump file. Without the check, functions that had been
            # --- installed via a decorator would be installed an extra time since the source
            # --- of the function contains the decoration (which is activated when the source
            # --- is exec'd).
            if f.__name__ not in self.funcs:
                self.funcs.append(f)
        else:
            self.funcs.append(f)

    def uninstallfuncinlist(self,f):
        """Uninstall the specified function"""
        # --- An element by element search is needed
        # --- f can be a function or method object, or a name (string).
        # --- Note that method objects can not be removed by name.
        funclistcopy = copy.copy(self.funcs)
        for func in funclistcopy:
            if f == func:
                self.funcs.remove(f)
                break
            elif isinstance(func,list) and isinstance(f,types.MethodType):
                object = self._getmethodobject(func)
                if f.__self__ is object and f.__name__ == func[1]:
                    self.funcs.remove(func)
                    break
            elif isinstance(func,str):
                if f.__name__ == func:
                    self.funcs.remove(func)
                    break
            elif isinstance(f,str):
                if isinstance(func,str): funcname = func
                elif isinstance(func,list): funcname = None
                else:                        funcname = func.__name__
                if f == funcname:
                    self.funcs.remove(func)
                    break

        # check that a function was removed
        if len(self.funcs) == len(funclistcopy):
            raise Exception(f'Warning: no function, {f}, had been installed')

        # if there are no functions left, remove the C callback
        if not self.hasfuncsinstalled():
            self.clearlist()

    def isinstalledfuncinlist(self,f):
        """Checks if the specified function is installed"""
        # --- An element by element search is needed
        funclistcopy = copy.copy(self.funcs)
        for func in funclistcopy:
            if f == func:
                return 1
            elif isinstance(func,list) and isinstance(f,types.MethodType):
                object = self._getmethodobject(func)
                if f.__self__ is object and f.__name__ == func[1]:
                    return 1
            elif isinstance(func,str):
                if f.__name__ == func:
                    return 1
        return 0

    def callfuncsinlist(self,*args,**kw):
        """Call the functions in the list"""
        bb = time.time()
        for f in self.callbackfunclist():
            #barrier()
            t1 = time.time()
            f(*args,**kw)
            #barrier()
            t2 = time.time()
            # --- For the timers, use the function (or method) name as the key.
            self.timers[f.__name__] = self.timers.get(f.__name__,0.) + (t2 - t1)
        aa = time.time()
        return aa - bb

#=============================================================================

callback_instances = {
    "beforeInitEsolve": {},
    "afterInitEsolve": {},
    "afterinit": {},
    "beforecollisions": {},
    "aftercollisions": {},
    "beforeEsolve": {},
    "poissonsolver": {'singlefunconly': True}, # external Poisson solver
    "afterEsolve": {},
    "beforedeposition": {},
    "afterdeposition": {},
    "particlescraper": {},
    "particleloader": {},
    "beforestep": {},
    "afterstep": {},
    "afterdiagnostics": {},
    "afterrestart": {},
    "oncheckpointsignal": {},
    "onbreaksignal": {},
    "particleinjection": {},
    "appliedfields": {}
}

# --- Now create the actual instances.
for key, val in callback_instances.items():
    callback_instances[key] = CallbackFunctions(name=key, **val)

def installcallback(name, f):
    "Adds a function to the list of functions called by this callback"
    callback_instances[name].installfuncinlist(f)

def uninstallcallback(name, f):
    "Removes the function from the list of functions called by this callback"
    callback_instances[name].uninstallfuncinlist(f)

def isinstalled(name, f):
    "Checks if the function is called by this callback"
    return callback_instances[name].isinstalledfuncinlist(f)

#=============================================================================

def printcallbacktimers(tmin=1.,lminmax=False,ff=None):
    """Prints timings of installed functions.
    - tmin=1.: only functions with time greater than tmin will be printed
    - lminmax=False: If True, prints the min and max times over all processors
    - ff=None: If given, timings will be written to the file object instead of stdout
    """
    if ff is None: ff = sys.stdout
    for c in callback_instances.values():
        for fname, time in c.timers.items():
            #vlist = numpy.array(gather(time))
            vlist = numpy.array([time])
            #if me > 0: continue
            vsum = numpy.sum(vlist)
            if vsum <= tmin: continue
            vrms = numpy.sqrt(max(0.,numpy.sum(vlist**2)/len(vlist) - (numpy.sum(vlist)/len(vlist))**2))
            npes = 1. # Only works for one processor
            ff.write('%20s %s %10.4f  %10.4f %10.4f'%(c.name,fname,vsum,vsum/npes,vrms))
            if lminmax:
                vmin = numpy.min(vlist)
                vmax = numpy.max(vlist)
                ff.write('  %10.4f  %10.4f'%(vmin,vmax))
            it = libwarpx.libwarpx_so.warpx_getistep(0)
            if it > 0:
                ff.write('   %10.4f'%(vsum/npes/(it)))
            ff.write('\n')

#=============================================================================

# ----------------------------------------------------------------------------
def callfrombeforeInitEsolve(f):
    installbeforeInitEsolve(f)
    return f
def installbeforeInitEsolve(f):
    "Adds a function to the list of functions called before the initial E solve"
    _beforeInitEsolve.installfuncinlist(f)
def uninstallbeforeInitEsolve(f):
    "Removes the function from the list of functions called before the initial E solve"
    _beforeInitEsolve.uninstallfuncinlist(f)
def isinstalledbeforeInitEsolve(f):
    "Checks if the function is called before the initial E solve"
    return _beforeInitEsolve.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromafterInitEsolve(f):
    installafterInitEsolve(f)
    return f
def installafterInitEsolve(f):
    "Adds a function to the list of functions called after the initial E solve"
    _afterInitEsolve.installfuncinlist(f)
def uninstallafterInitEsolve(f):
    "Removes the function from the list of functions called after the initial E solve"
    _afterInitEsolve.uninstallfuncinlist(f)
def isinstalledafterInitEsolve(f):
    "Checks if the function is called after the initial E solve"
    return _afterInitEsolve.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromafterinit(f):
    installcallback('afterinit', f)
    return f

# ----------------------------------------------------------------------------
def callfrombeforecollisions(f):
    installcallback('beforecollisions', f)
    return f

# ----------------------------------------------------------------------------
def callfromaftercollisions(f):
    installcallback('aftercollisions', f)
    return f

# ----------------------------------------------------------------------------
def callfrombeforeEsolve(f):
    installcallback('beforeEsolve', f)
    return f

# ----------------------------------------------------------------------------
def callfrompoissonsolver(f):
    installcallback('poissonsolver', f)
    return f

# ----------------------------------------------------------------------------
def callfromafterEsolve(f):
    installcallback('afterEsolve', f)
    return f

# ----------------------------------------------------------------------------
def callfrombeforedeposition(f):
    installcallback('beforedeposition', f)
    return f

# ----------------------------------------------------------------------------
def callfromafterdeposition(f):
    installcallback('afterdeposition', f)
    return f

# ----------------------------------------------------------------------------
def callfromparticlescraper(f):
    installcallback('particlescraper', f)
    return f

# ----------------------------------------------------------------------------
def callfromparticleloader(f):
    installcallback('particleloader', f)
    return f

# ----------------------------------------------------------------------------
def callfrombeforestep(f):
    installcallback('beforestep', f)
    return f

# ----------------------------------------------------------------------------
def callfromafterstep(f):
    installcallback('afterstep', f)
    return f

# ----------------------------------------------------------------------------
def callfromafterdiagnostics(f):
    installcallback('afterdiagnostics', f)
    return f

# ----------------------------------------------------------------------------
def oncheckpointsignal(f):
    installcallback('oncheckpointsignal', f)
    return f

# ----------------------------------------------------------------------------
def onbreaksignal(f):
    installonbreaksignal(f)
    return f
def installonbreaksignal(f):
    "Adds a function to the list of functions called on a break signal"
    _onbreaksignal.installfuncinlist(f)
def uninstallonbreaksignal(f):
    "Removes the function from the list of functions called on a break signal"
    _onbreaksignal.uninstallfuncinlist(f)
def isinstalledonbreaksignal(f):
    "Checks if the function is called on a break signal"
    return _onbreaksignal.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromparticleinjection(f):
    installcallback('particleinjection', f)
    return f
