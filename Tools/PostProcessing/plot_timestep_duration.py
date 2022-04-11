#!/usr/bin/env python3

import re
import matplotlib.pyplot as plt
import numpy as np
import sys

def extract_data(filename):
    
    regex_step = re.compile(
        r"STEP [0-9]* ends.*\n.* Avg\. per step = ([0-9]*[.])?[0-9]+ s", re.MULTILINE)
    
    string_data = []
    
    print("Processing " + output_log_file_name + " ...", end='')
    with open(output_log_file_name) as f:
        text = f.read()
   
        string_data = [s.group(0) for s in regex_step.finditer(text)]
    
    regex_real = re.compile(
        r" -?[\d.]+(?:e-?\d+)?", re.MULTILINE)
    
    time_data = np.zeros([len(string_data), 6])
    for i, ss in enumerate(string_data):
        numbers = regex_real.findall(ss)
        time_data[i,:] = np.array(numbers) 
    
    print("...done!")
    return time_data

def plot_timestep_duration(time_data, name):
    fig_name = name + "_ts_duration.png"
    print("Generating " + fig_name + "...", end='')

    plt.rcParams.update({'font.size': 20})
    plt.rcParams['axes.linewidth'] = 3

    f, ax = plt.subplots(figsize=(12,6))

    ax.set_ylabel("timestep duration [s]")
    ax.set_xlabel("step [#]")

    ax.semilogy(time_data[:,0], time_data[:,4])

    ax.spines['bottom'].set_color('gray')
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_color('gray')
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    plt.savefig(fig_name, transparent=False, dpi=300)
    print("...done!")
    
    
def plot_cumulative_duration(time_data, name):
    fig_name = name + "_cumulative_duration.png"
    print("Generating " + fig_name + "...", end='')
   
    plt.rcParams.update({'font.size': 20})
    plt.rcParams['axes.linewidth'] = 3

    f, ax = plt.subplots(figsize=(12,6))

    ax.set_ylabel("cumulative duration [s]")
    ax.set_xlabel("step [#]")

    ax.plot(time_data[:,0], np.cumsum(time_data[:,4]))

    ax.spines['bottom'].set_color('gray')
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_color('gray')
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    plt.savefig(name + "_cumulative_duration.png", transparent=False, dpi=300)
    print("...done!")    

output_log_file_name = sys.argv[1]
time_data = extract_data(output_log_file_name)

plot_timestep_duration(time_data, output_log_file_name)
plot_cumulative_duration(time_data, output_log_file_name)
