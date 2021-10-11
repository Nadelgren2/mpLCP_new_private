###############################################################################
#
#   Author:         Nathan Adelgren
#   Affiliation:    Andlinger Center For Energy and the Environment
#                   Princeton University
#
#   Purpose:        Declare a counter class to be used to keep track of the
#                   number of invariancy regions left to process.
#
################################################################################

import multiprocessing

class Counter(object):
    def __init__(self, initval=0):
        self.val = multiprocessing.Value('i', initval)
        self.lock = multiprocessing.Lock()

    def Increment(self):
        with self.lock:
            self.val.value += 1
            
    def Decrement(self):
        with self.lock:
            self.val.value -= 1

    def Value(self):
        with self.lock:
            return self.val.value
