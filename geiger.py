# -*- coding: utf-8 -*-
"""
Created on Tue May 16 23:49:48 2023

@author: Charlie Perkins
"""


import numpy as np
from scipy.io import wavfile
import matplotlib.pyplot as plt
import getopt
import sys


def get_args():
    argumentList = sys.argv[1:]
    options = "hi:o:d:Dan:t:w:I:f:s:b:"
    long_options = ["help", "in_dir=", "out_dir=", "disp_name=", "Debug", "analyse", "num=", "time=", "width=", "Isotopes=", "filename=", "significance=", "background="]
    in_dir = "./"
    out_dir = "./"
    disp_name = "samplename"
    debug = False
    analyse = False
    use_num = True
    num = 20
    time = 3
    width = 0.00041
    isotopes = 1
    filename = "sample.wav"
    significance = 0.75
    use_background = False
    background = "background.wav"
    try:
        arguments, values = getopt.getopt(argumentList, options, long_options)
        if ("-n" or "--name") in arguments:
            if ("-t" or "--time") in arguments:
                print("ERROR: Must give ONLY -n (--num) OR -t (--time) not both")
                raise sys.exit() # SWAP FOR OPTION TO DEFINE WITHOUT QUITTING
        if ("-a" or "--analyse") in arguments:
            if ("-I" or "--Isotopes") in arguments:
                pass
            else:
                print("ERROR: Must give BOTH -a (--analyse) AND -I (--Isotopes) or NEITHER")
                raise sys.exit() # SWAP FOR OPTION TO DEFINE WITHOUT QUITTING
        for argument, value in arguments:
            # THESE ALL NEED VALIDATING
            if argument in ("-h", "help"):
                disp_help()
            elif argument in ("-i", "--in_dir"):
                in_dir = value
            elif argument in ("-o", "--out_dir"):
                out_dir = value
            elif argument in ("-d", "--disp_name"):
                disp_name = value
            elif argument in ("-D", "--Debug"):
                debug = True
            elif argument in ("-a", "--analyse"):
                analyse = True
            elif argument in ("-n", "--num"):
                num = int(value)
            elif argument in ("-t", "--time"):
                time = float(value)
                use_num = False
            elif argument in ("-w", "--width"):
                width = float(value)
            elif argument in ("-I", "--Isotopes"):
                isotopes = int(value)
            elif argument in ("-f", "--filename"):
                filename = value
            elif argument in ("-s", "--significance"):
                significance = float(value)
            elif argument in ("-b", "--background"):
                use_background = True
                background = value
    except getopt.error as err:
        print("ERROR: Argument error")
        print(str(err))
    return (in_dir, out_dir, disp_name, debug, analyse, use_num, num, time, width, isotopes, filename, significance, use_background, background)


def disp_help():
    print("""
_______________________________________          
  ___     _                           
 / __|___(_)__ _ ___ _ _              
| (_ / -_) / _` / -_) '_|             
 \___\___|_\__, \___|_|               
    /_\  _ |___/ _| |_  _ ___ ___ _ _ 
   / _ \| ' \/ _` | | || (_-</ -_) '_|
  /_/ \_\_||_\__,_|_|\_, /__/\___|_|  
                     |__/             
_______________________________________

Arguments:
    -h, --help         : Print help (but you already knew that...)
    -f, --filename     : Set the filename, default - sample.wav
    -b, --background   : Set the background filename, default - background.wav
    -i, --in_dir       : Set the input file directory, default - ./
    -o, --out_dir      : Set the output file directory, default - ./
    -d, --disp_name    : Set the filename and display names, default - sample
    -D, --Debug        : Select debug (full) mode (dump all data)
    -a, --analyse      : Guess the isotope based on decay curve (requires -I)
    -I, --Isotopes     : Define number active of isotopes in sample
    -n, --num          : Define the number of bins in file time
    OR
    -t, --time         : Define the duration of one bin (seconds)
    -w, --width        : Define the width of one count (seconds)
    -s, --significance : Define the significance for a count (0-1)
        """)
    raise sys.exit()
        
        
def get_data(in_dir, filename):
    # NEEDS DATA VALIDATION
    rate, data = wavfile.read(in_dir + filename)
    return rate, data


def process(bulk, fore_rate, fore_data, *args):
    in_dir, out_dir, disp_name, debug, analyse, use_num, num, time, width, isotopes, filename, significance, use_background, background = bulk
    noise = 0
    back_cpm = 0
    if use_num:
        time = len(fore_data) / (fore_rate * num)
    else:
        num = len(fore_data) / round(time * fore_rate, 0)
        
    if use_background:
        back_rate, back_data = args
        noise, back_cpm = get_noise(bulk, back_rate, back_data)
        
    counts = count_peaks(width, fore_rate, fore_data, significance, noise)
    
    discretised_counts = np.histogram(counts, num)[0]
    activity = discretised_counts * 60 / time
    cpm = countrate(counts, fore_rate) - back_cpm
    return counts, discretised_counts, activity, cpm

def get_noise(bulk, rate, data):
    in_dir, out_dir, disp_name, debug, analyse, use_num, num, time, width, isotopes, filename, significance, use_background, background = bulk
    noise = np.mean(np.abs(data))
    counts = count_peaks(width, rate, data, significance, noise)
    cpm = countrate(counts, rate)
    print("The backround activity was:", round(cpm, 0), "cpm with", round(noise, 0), "noise.")
    return noise, cpm
	
def countrate(counts, rate):
	return np.mean(counts) * rate * 60


def count_peaks(width, rate, data, significance, noise):
    data = np.abs(data) - np.abs(noise)
    c_max = np.max(data)
    crit_vol = significance * c_max
    samples = int(round(width * rate, 0))
    counts = []
    i = 0
    while i < len(data)-samples:
        if data[i] > crit_vol:
            counts.append(1)
            counts.extend([0] * (samples - 1))
            i = i + samples
        else:
            counts.append(0)
            i += 1
    return counts

def plot(bulk, activity):
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))
    # ax2 = ax.twinx()
    ax.plot(activity, color="red", linestyle='-')
    ax.grid()
    plt.tight_layout()
    plt.show()
    plt.clf()
            
            
if __name__ == "__main__":
    args = get_args()
    in_dir, out_dir, disp_name, debug, analyse, use_num, num, time, width, isotopes, filename, significance, use_background, background = args
    fore_rate, fore_data = get_data(in_dir, filename)
    if use_background:
        back_rate, back_data = get_data(in_dir, background)
        counts, discretised_counts, activity, cpm = process(args, fore_rate, fore_data, back_rate, back_data)
    else:
        counts, discretised_counts, activity, cpm = process(args, fore_rate, fore_data)
    # plot(args, activity)
    print("The sample activity was:", round(cpm, 0), "cpm (", round(cpm / 60, 0), "Bq )")
    
    