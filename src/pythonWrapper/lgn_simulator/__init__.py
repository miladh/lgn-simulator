print "hello world!"
def __bootstrap__():
    global __bootstrap__, __loader__, __file__
    import sys, pkg_resources, imp
    __file__ = pkg_resources.resource_filename("lgn_simulator", 'lgn_simulator.so')
    __loader__ = None; del __bootstrap__, __loader__
    imp.load_dynamic("lgn_simulator",__file__)
__bootstrap__()
