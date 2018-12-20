#!/usr/bin/env python
<<<<<<< HEAD

import numpy as np
=======
>>>>>>> add_parser

import numpy as np
import win32pipe as wp
import win32file as wf
import matplotlib.pyplot as plt
import time
import datetime
import subprocess

from sys import executable
from subprocess import Popen, CREATE_NEW_CONSOLE


np.set_printoptions(threshold=np.nan) # Allow printing of large string

# Returns buffer size int
def get_buffer_size(pipe_handle, num_bytes):
    buffer_size = int.from_bytes(wf.ReadFile(pipe_handle, num_bytes)[1], byteorder='big')
    return buffer_size

def trans(inp,trans):
    return (np.square(np.real(np.matmul(trans,inp)))*1).astype(np.uint8)

################################### INITIALIZE #############################################

NUM_BYTES_BUFFER = 4

k,l = 32*6, 24*6 # dimensions of input masks (optimal ~ 33,33)
m,n = 32, 24 # dimensions of the output field. (optimal ~ 33,33)
N = k*l # input matrix pixels (modes).
M = m*n # output matrix pixels (modes).


S = 30 # population of generated input phase masks matrices. (optimal ~ 30)
I = 1/np.sqrt(N) # intensity coefficient.
I0 = 0. # average intensity of output field before optimization.
Im = 0. # intensity of target output mode.

phsteps = 256 # number of discrete phase values
phdist = np.arange(0,2*np.pi+.00000001,2*np.pi/phsteps) # distribution of discrete phase values

tmat = np.empty((M,N),complex) # Transmission matrix. Size MxN.
tshift = np.random.randint(0,255,size=(M,N)) # Transmission matrix phase value indices.
tmat = np.exp(1j*phdist[tshift])


pipe_in = wp.CreateNamedPipe("\\\\.\\pipe\\LABVIEW_IN",
                       wp.PIPE_ACCESS_DUPLEX,
                       wp.PIPE_TYPE_BYTE | wp.PIPE_WAIT,
                       1,65536, 65536,300,None)

pipe_out = wp.CreateNamedPipe("\\\\.\\pipe\\LABVIEW_OUT",
                       wp.PIPE_ACCESS_DUPLEX,
                       wp.PIPE_TYPE_BYTE | wp.PIPE_WAIT,
                       1,65536, 65536,300,None)

print('ready....')
print('calling waveopt_interface.py....')

p = subprocess.Popen('python waveopt_interface.py --SLM_WIDTH 192 --SLM_HEIGHT 144 --PLOT True',
                     creationflags=CREATE_NEW_CONSOLE)

print('...done')
for i in range(300000):
    wp.ConnectNamedPipe(pipe_in, None)
    wp.ConnectNamedPipe(pipe_out, None)

    read_pipe = wf.ReadFile(pipe_in, get_buffer_size(pipe_in, NUM_BYTES_BUFFER))
    read_array = list(read_pipe[1])
    print('pipe read...')

    ### Plot Masks after reading ####
##    plt.matshow(np.asarray(read_array).reshape(l,k))
##    plt.show()

    out_field = trans(I*np.exp(1j*phdist[read_array]),tmat)[:5]
##    print(out_field)
##    out_field = read_array
    data = bytearray(len(bytearray(out_field)).to_bytes(NUM_BYTES_BUFFER, byteorder='big')) + bytearray(out_field)
    wf.WriteFile(pipe_out,data)

    









    
