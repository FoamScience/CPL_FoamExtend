#!/usr/bin/env python2
import sys, os
sys.path.append(os.path.abspath('..'))
from test_forces import compare_forces 

if __name__ == "__main__":
    compare_forces(1e-6, steps=2, dt=2, mode="debug")
