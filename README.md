# Gieger Counter WAV File Analyser
This is a command line python file to analyse activity of radioactive samples by taking the audio of a geiger counter.
Currently this is limited to giving the cpm and estimated activity, but it may in the future graph activity and more.

## File Types
It is important to note that 16-bit Integer WAV files work best. The sample rate is not limited.

## Usage
Navigate to the directory of geiger.py.

For a list of commands, use:
'''bash
Python geiger.py -h
'''

The default foreground filename is sample.wav, and the background background.wav.