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
from matplotlib.ticker import AutoMinorLocator
from pathlib import Path


def get_args():
    def value_error(argument, value):
        print("WARNING:", argument, "cannot have value", value, "... using default instead...")
    argumentList = sys.argv[1:]
    options = "hi:o:d:Dan:t:w:I:f:s:b:H:"
    long_options = ["help", "in_dir=", "out_dir=", "disp_name=", "Debug", "analyse", "num=", "time=", "width=", "Isotopes=", "filename=", "significance=", "background=", "Histogram="]
    in_dir = "./"
    out_dir = "./"
    disp_name = "defaultsamplename"
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
    histogram = 10
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
                if int(value) == 0:
                    value_error(argument, value)
                else:
                    num = np.abs(int(value))
            elif argument in ("-t", "--time"):
                if float(value) == 0:
                    value_error(argument, value)
                else:
                    time = np.abs(float(value))
                    use_num = False
            elif argument in ("-w", "--width"):
                if float(value) == 0:
                    value_error(argument, value)
                else:
                    width = np.abs(float(value))
            elif argument in ("-I", "--Isotopes"):
                if float(value) == 0:
                    value_error(argument, value)
                else:
                    isotopes = np.abs(int(value))
            elif argument in ("-f", "--filename"):
                if value == "":
                    value_error(argument, value)
                elif ".wav" in value:
                    filename = value
                else:
                    filename = value + ".wav"
            elif argument in ("-s", "--significance"):
                if float(value) >= 1:
                    print("WARNING: Setting", argument, "=", value, "may cause unexpected behaviour. Continuing...")
                significance = float(value)
            elif argument in ("-b", "--background"):
                use_background = True
                if value == "":
                    value_error(argument, value)
                elif ".wav" in value:
                    background = value
                else:
                    background = value + ".wav"   
            elif argument in ("-H", "--Histogram"):
                if float(value) <= 0:
                    value_error(argument, value)
                else:
                    histogram = int(value)
    except getopt.error as err:
        print("ERROR: Argument error. Attempting to continue...")
        print(str(err))
    except ValueError as err:
        print("ERROR: Argument error. Attempting to continue...")
        print(str(err))
    return (in_dir, out_dir, disp_name, debug, analyse, use_num, num, time, width, isotopes, filename, significance, use_background, background, histogram)


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
    -n, --num          : Define the number of bins for activity plot
    OR
    -t, --time         : Define the duration of one bin in activity plot (seconds)
    -w, --width        : Define the width of one count (seconds)
    -s, --significance : Define the significance for a count (0-1)
    -H, --Histogram    : Define the number of bins for the histogram
        """)
    raise sys.exit()
        
        
def get_data(in_dir, filename):
    # NEEDS DATA VALIDATION
    rate, data = wavfile.read(in_dir + filename)
    return rate, data


def process(bulk, fore_rate, fore_data, *args):
    in_dir, out_dir, disp_name, debug, analyse, use_num, num, time, width, isotopes, filename, significance, use_background, background, histogram = bulk
    noise = 0
    back_cpm = 0        
    if use_background:
        back_rate, back_data = args
        noise, back_cpm = get_noise(bulk, back_rate, back_data)
        
    counts = count_peaks(width, fore_rate, fore_data, significance, noise)
    
    discretised_counts = discretise(bulk, fore_rate, counts)
    dist = np.histogram(discretised_counts)
    dist_x = dist[1]
    dist_y = dist[0]
    activity = discretised_counts * 60 / time
    cpm = countrate(counts, fore_rate) - back_cpm
    times = []
    for i in range(len(activity)):
        times.append(i*time)
    return counts, discretised_counts, activity, cpm, dist_x, dist_y, times

def get_noise(bulk, rate, data):
    in_dir, out_dir, disp_name, debug, analyse, use_num, num, time, width, isotopes, filename, significance, use_background, background, histogram = bulk
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

def plot(bulk, activity, cpm, dist_x, dist_y, times):
    in_dir, out_dir, disp_name, debug, analyse, use_num, num, time, width, isotopes, filename, significance, use_background, background, histogram = bulk
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(12, 10))
    
    ax1.plot(times, activity, color="crimson", linestyle='-', label="Sample Activity")
    ax1.axhline(y=cpm, color='mediumaquamarine', linestyle='dotted', label="Mean Activity")
    ax1.set_title("Activity", fontsize=18)
    ax1.set_xlabel("Time, s")
    ax1.set_ylabel("Activity, CPM")
    ax1.grid(which='major')
    ax1.grid(which='minor', linewidth='0.3')
    ax1.minorticks_on()
    ax1.tick_params(which='minor', bottom=False, left=False)
    ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax1.legend()
    
    y, x, _ = ax2.hist(activity, bins=histogram, color="green")
    y = -1 * np.max(y)
    x = np.max(x)
    ax2.set_title("Activity Distribution", fontsize=18)
    ax2.set_xlabel("Activity, CPM")
    ax2.set_ylabel("Frequency")
    ax2.grid(which='major')
    ax2.grid(which='minor', linewidth='0.3')
    ax2.minorticks_on()
    ax2.tick_params(which='minor', bottom=False, left=False)
    ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    
    ax3.text(0.03, 1.1, "Sample: " + disp_name, fontsize=12)
    
    ax3.text(0.03, 0.93, "Count Width:", fontsize=10)
    ax3.text(0.16, 0.93, str(round(width, 6)) + " s", fontsize=10)
    ax3.text(0.03, 0.86, "CPM:", fontsize=10)
    ax3.text(0.16, 0.86, str(round(cpm, 0)), fontsize=10)
    ax3.text(0.03, 0.79, "Estimated Activity:", fontsize=10)
    ax3.text(0.16, 0.79, str(round(cpm / 60, 3)) + " Bq", fontsize=10)
    ax3.text(0.03, 0.72, "# Distribution Bins:", fontsize=10)
    ax3.text(0.16, 0.72, str(histogram), fontsize=10)
    
    ax3.text(0.3, 0.93, "# Isotopes:", fontsize=10)
    ax3.text(0.45, 0.93, str(isotopes), fontsize=10)
    ax3.text(0.3, 0.86, "Activity Bin Time: ", fontsize=10)
    ax3.text(0.45, 0.86, str(round(time, 3)) + " s", fontsize=10)
    ax3.text(0.3, 0.79, "Background Removed:", fontsize=10)
    ax3.text(0.45, 0.79, str(use_background), fontsize=10)
    ax3.text(0.3, 0.72, "Count Significance:", fontsize=10)
    ax3.text(0.45, 0.72, str(significance * 100) + " %", fontsize=10)
    
    ax3.text(0, 0.03, "Geiger Analyser on GitHub from https://github.com/Chaddyfynn/geiger-analyser",
             fontsize=10)
    ax3.text(0.9, 0.03, "v0.2-beta", fontsize=10)
    
    ax3.axis('off')

    
    filepath = path_checker(disp_name, ".png", out_dir)
    plt.tight_layout()
    plt.savefig(filepath, dpi=800)
    plt.show()
    plt.clf()
    
def discretise(bulk, rate, counts):
    in_dir, out_dir, disp_name, debug, analyse, use_num, num, time, width, isotopes, filename, significance, use_background, background, histogram = bulk
    split_counts = np.array_split(counts, num)
    discretised = []
    for i in range(len(split_counts)):
        discretised.append(np.sum(split_counts[i]))
    discretised = np.asarray(discretised)
    return discretised

def path_checker(ideal_filename, extension, address):
    filename = ideal_filename
    number = int(0)
    path = Path(address + filename + extension)  # Initial Path

    while path.is_file():
        number = number + int(1)
        number_append = str(number)
        filename = ideal_filename + '_' + number_append
        path = Path(address + filename + extension)
    return address + filename + extension
            
            
if __name__ == "__main__":
    args = get_args()
    in_dir, out_dir, disp_name, debug, analyse, use_num, num, time, width, isotopes, filename, significance, use_background, background, histogram = args
    fore_rate, fore_data = get_data(in_dir, filename)
    if use_num:
        time = len(fore_data) / (fore_rate * num)
    else:
        num = len(fore_data) / round(time * fore_rate, 0)
    if use_background:
        back_rate, back_data = get_data(in_dir, background)
        counts, discretised_counts, activity, cpm, dist_x, dist_y, times = process(args, fore_rate, fore_data, back_rate, back_data)
    else:
        counts, discretised_counts, activity, cpm, dist_x, dist_y, times = process(args, fore_rate, fore_data)
    plot(args, activity, cpm, dist_x, dist_y, times)
    print("The sample activity was:", round(cpm, 0), "cpm (", round(cpm / 60, 0), "Bq )")
    
    