#!/usr/bin/env python
"""Genetic Algorithm optimizer for phase modulated wavefront.

This algorithm is approximately optimized when number of segments ~ 1000, number of generations pipe_out ~ 2500, and population of
input masks POP ~ 30. Adjust fitness function and PARAMETERS under MAIN.
"""

import argparse
import numpy as np
import win32pipe as wp
import win32file as wf
import matplotlib.pyplot as plt
import time
import datetime
import sys

__author__ = "Jesse Weller"
__copyright__ = "Copyright 2018, Jesse Weller, All rights reserved"
__version__ = "1.1"
__email__ = "wellerj@oregonstate.edu"

######################### FUNCTIONS ###################################

def flatten_mask(mask):
    """Flatten mask. Return 1d array."""
    return np.stack(np.stack(mask,axis=1),axis=2).flatten()

def get_buffer_size(pipe_handle, num_bytes):
    """Return buffer size as int."""
    buffer_size = int.from_bytes(wf.ReadFile(pipe_handle, num_bytes)[1], byteorder='little', signed=False)
    return buffer_size

def encode_mask(mask, num_bytes_buffer):
    """Return buffer size and mask as bytearray."""
    return bytearray(len(mask).to_bytes(num_bytes_buffer, byteorder='little', signed=False)) + bytearray(mask)


def plot_mask(mask,mask_width,mask_height):
    """Plots the flattened and then reshaped phase values of the mask"""
    plt.matshow(flatten_mask(mask).reshape(mask_height,mask_width))
    plt.show()

def get_output_fields(input_masks,
                      segment_width,
                      segment_height,
                      segment_rows,
                      segment_cols,
                      pipe_handle_in,
                      pipe_handle_out,
                      num_bytes_buffer):
    """Transmit mask pixel data through pipe to apparatus. Return list of ROI field arrays."""
    t0=time.time()
    roi_placeholder = []
    ###### PLOT MASK ####
    # plot_mask(input_masks[0],segment_width*segment_cols,segment_height*segment_rows)
    
    for mask in input_masks:
        wf.WriteFile(pipe_handle_out,
                     encode_mask(flatten_mask(mask), num_bytes_buffer))

        read_pipe = wf.ReadFile(pipe_handle_in, get_buffer_size(pipe_handle_in, num_bytes_buffer))
        read_array = list(read_pipe[1])
        roi_placeholder.append(fitness(read_array))
    print(time.time()-t0)
    return roi_placeholder

def ranksort(output_fields, input_masks, generation, prev_time):
    """Return list of input masks sorted by fitness value"""
    ranks=[]
    sorted_input=[]

    ranks=np.argsort(output_fields)
    for i in range(len(input_masks)):
        sorted_input.append(input_masks[ranks[i]])
    if generation % 1 == 0:
        print("Generation ",generation,": output focus=",output_fields[ranks[-1]])
        print("time: ",str(datetime.timedelta(seconds=time.time()-prev_time)))
    return sorted_input 

def breed(parent1, parent2, mutate, phase_vals):
    """Breed two "parent" masks and return new mutated "child" input mask array."""
    shape = parent1.shape
    rand_bool = np.random.choice((0,1),size=shape[0]*shape[1]).reshape(shape[0],shape[1])
    child = parent1.copy()
    for row in range(shape[0]):
        for col in range(shape[1]):
            if rand_bool[row,col]>0:
                child[row,col,:,:] = parent2[row,col,:,:]
    for i in range(mutate):
        child[np.random.randint(0,shape[0]),np.random.randint(0,shape[1]),:,:] = phase_vals[np.random.randint(0,len(phase_vals))]
    return child

def select_parents(parent_index_dist, probabilities):
    """Return an array of two randomly chosen ints."""
    return np.random.choice(parent_index_dist, size=2, p=probabilities, replace=False)

def init_masks(num_masks, segment_rows, segment_cols, segment_height, segment_width, phase_vals):
    """Return list of initialized input masks."""
    input_masks = []
    for mask in range(num_masks):
        newmask = np.empty((segment_rows, segment_cols, segment_height, segment_width), np.uint8)
        for row in range(segment_rows):
            for col in range(segment_cols):
                newmask[row,col,:,:] = phase_vals[np.random.randint(0,len(phase_vals))]
        input_masks.append(newmask)
        
    return input_masks

def get_num_segments(slm_width, slm_height, segment_width, segment_height):
    """Return number of segments in each mask as int."""
    mask_width = int(slm_width/segment_width)
    mask_height = int(slm_height/segment_height)
    return int(mask_width*mask_height)

def wavefront_phase_optimizer(pipe_in_handle,
                              pipe_out_handle,
                              bytes_buffer_size,
                              num_segments,
                              num_masks,
                              segment_width,
                              segment_height,
                              segment_rows,
                              segment_cols,
                              generations,
                              num_phase_vals,
                              mutate_initial_rate,
                              mutate_final_rate,
                              mutate_decay_factor,
                              prev_time,
                              plot_bool):
    """Optimize wavefront over generations. Return final mask array, list of output values

    Note: iterates over generations. For each generation: generates output intensity field matrices,
    sorts input array based on rank. Generates "G" new child masks. For each child, randomly chooses
    two parents for breeding based on probability distribution "prob", and breeds them. Culls lowest
    ranked input masks and replaces them with child masks.
    """
    
    phase_vals = []
    for val in range(num_phase_vals):
        phase_vals.append(np.ones((segment_height,segment_width))*val) # Distribution of phase values
    num_childs = int(round(num_masks/2)) # Number of new children per generation
    parent_index_dist = np.arange(num_childs,num_masks,1) # Distribution of index values for choosing parents.
    pstep = 1/sum(np.arange(0,num_masks-num_childs+1,1)) # Increment for generating probability distribution.
    parent_probability_dist = np.arange(pstep,(num_masks-num_childs+1)*pstep,pstep) # Probability distribution for parent selection.

    input_masks = init_masks(num_masks, segment_rows, segment_cols, segment_height, segment_width, phase_vals)
    
    pipe_in = wf.CreateFile(pipe_in_handle,
                 wf.GENERIC_READ | wf.GENERIC_WRITE,
                 0, None,
                 wf.OPEN_EXISTING,
                 0,None)
                             
    pipe_out = wf.CreateFile(pipe_out_handle,
                 wf.GENERIC_READ | wf.GENERIC_WRITE,
                 0, None,
                 wf.OPEN_EXISTING,
                 0,None)
                              
    max_output_vals = []
    for gen in range(generations):
        slm_outputs = []
        num_mutations = int(round(num_segments * ((mutate_initial_rate - mutate_final_rate)
                                    * np.exp(-gen / mutate_decay_factor)
                                    + mutate_final_rate)))
        child_placeholder = []

        slm_outputs = get_output_fields(input_masks,
                                        segment_width,
                                        segment_height,
                                        segment_rows,
                                        segment_cols,
                                        pipe_in,
                                        pipe_out,
                                        bytes_buffer_size)
        max_output_vals.append(np.max(slm_outputs)) # record max output values for visualization
        input_masks = ranksort(slm_outputs,input_masks,gen, prev_time) # sort input masks by rank
        if plot_bool is True:
            plt.plot(max_output_vals)
            plt.pause(0.05)    
        for child in range(num_childs): # make "G" new masks through breeding
            parents = select_parents(parent_index_dist, parent_probability_dist)
            input_masks[child] = breed(input_masks[parents[0]],input_masks[parents[1]],num_mutations,phase_vals)
    final_mask = input_masks[-1]
    return final_mask, max_output_vals


def fitness(output_field):
    """Return fitness value of mask as float.

    Note: Adjust fitness function to suit your optimization process.
    """
    return np.mean(output_field)


################################# MAIN ##############################################

if __name__ == '__main__':   

    ############## Parse Command Line Arguments #############################   
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--PIPE_IN_HANDLE',
        type=str,
        default='\\\\.\\pipe\\LABVIEW_OUT',
        help="""\Input Pipe handle.\
      """
    )
    parser.add_argument(
        '--PIPE_OUT_HANDLE',
        type=str,
        default='\\\\.\\pipe\\LABVIEW_IN',
        help="""\Output Pipe handle.\
      """
    )        
    parser.add_argument(
        '--BYTES_BUFFER_SIZE',
        type=int,
        default=4,
        help="""\Pipe inputoutput buffer size.\
      """
    )  
    parser.add_argument(
        '--PLOT',
        type=bool,
        default=False,
        help="""\Turn on/off visualization of optimization.\
      """
    )
    parser.add_argument(
        '--SLM_WIDTH',
        type=int,
        default=1024,
        help="""\Pixel width of SLM.\
      """
    )
    parser.add_argument(
        '--SLM_HEIGHT',
        type=int,
        default=768,
        help="""\Pixel height of SLM.\
      """
    )
    parser.add_argument(
        '--SEGMENT_WIDTH',
        type=int,
        default=32,
        help="""\Pixel width of each segment (group of pixels on SLM).\
      """
    )
    parser.add_argument(
        '--SEGMENT_HEIGHT',
        type=int,
        default=24,
        help="""\Pixel height of each segment (group of pixels on SLM).\
      """
    )
       
    parser.add_argument(
        '--POP',
        type=int,
        default=30,
        help="""\Initial population of randomly generated phase masks in genetic algorithm.\
      """
    )
    parser.add_argument(
        '--GENS',
        type=int,
        default=1000,
        help="""\Number of generations to run genetic algorithm.\
      """
    )
    parser.add_argument(
        '--MUTATE_INITIAL_RATE',
        type=float,
        default=.1,
        help="""\Initial mutation rate for genetic algorithm.\
      """
    )
    parser.add_argument(
        '--MUTATE_FINAL_RATE',
        type=float,
        default=.013,
        help="""\Final mutation rate for genetic algorithm.\
      """
    )
    parser.add_argument(
        '--MUTATE_DECAY_FACTOR',
        type=float,
        default=650,
        help="""\Final mutation rate for genetic algorithm.\
      """
    )
    parser.add_argument(
        '--NUM_PHASE_VALS',
        type=int,
        default=256,
        help="""\Number of discrete phase values to be passed to SLM.\
      """
    )
    parser.add_argument(
        '--SAVE_PATH',
        type=str,
        default='/tmp/optimized_masks/optimized_mask.txt',
        help="""\Path of text file to save optimized mask.\
      """
    )
    args = parser.parse_args()
    
    
    PIPE_IN_HANDLE = args.PIPE_IN_HANDLE
    PIPE_OUT_HANDLE = args.PIPE_OUT_HANDLE
    BYTES_BUFFER_SIZE = args.BYTES_BUFFER_SIZE
    PLOT = args.PLOT # plots fitness values for each generation if True, no plot if False
                            
    SLM_WIDTH = args.SLM_WIDTH
    SLM_HEIGHT = args.SLM_HEIGHT
    SEGMENT_WIDTH = args.SEGMENT_WIDTH # SLM_WIDTH % SEGMENT_WIDTH must be 0
    SEGMENT_HEIGHT = args.SEGMENT_HEIGHT # SLM_HEIGHT % SEGMENT_HEIGHT must be 0
    POP = args.POP # Population of generated input phase masks. (optimal ~ 30)
    GENS = args.GENS # Number of generations to run algorithm. (optimal ~ 2000)

    MUTATE_INITIAL_RATE = args.MUTATE_INITIAL_RATE # (optimal ~ .1)
    MUTATE_FINAL_RATE = .args.MUTATE_FINAL_RATE # (optimal ~ .013)
    MUTATE_DECAY_FACTOR = args.MUTATE_DECAY_FACTOR # (optimal ~ 650)

    NUM_PHASE_VALS = args.NUM_PHASE_VALS
    num_segments = get_num_segments(SLM_WIDTH, SLM_HEIGHT, SEGMENT_WIDTH, SEGMENT_HEIGHT)
    segment_rows = int(SLM_HEIGHT/SEGMENT_HEIGHT)
    segment_cols = int(SLM_WIDTH/SEGMENT_WIDTH)
    
    time_start = time.time()

    optimized_mask, max_output_vals = wavefront_phase_optimizer(PIPE_IN_HANDLE,
                                                                PIPE_OUT_HANDLE,
                                                                BYTES_BUFFER_SIZE,
                                                                num_segments,
                                                                POP,
                                                                SEGMENT_WIDTH,
                                                                SEGMENT_HEIGHT,
                                                                segment_rows,
                                                                segment_cols,
                                                                GENS,
                                                                NUM_PHASE_VALS,
                                                                MUTATE_INITIAL_RATE,
                                                                MUTATE_FINAL_RATE,
                                                                MUTATE_DECAY_FACTOR,
                                                                time_start,
                                                                PLOT)
    time_end = time.time()
    print("Optimization Time: ",str(datetime.timedelta(seconds=time_end-time_start)))
    
    np.savetxt(args.SAVE_PATH, optimized_mask)
    
    plt.plot(max_output_vals)
    plt.show()
