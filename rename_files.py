# rename_files.py
# Rename files for ease of plotting
#
# Author Z. Wallace
# Created: 4.12.19

from subprocess import run
#from netCDF4 import Dataset, Variable

import numpy as np 
import matplotlib.pyplot as plt

# acquire list of file names
path_to_dir = "/Users/zach/Documents/Patagonia_Iron/data/NASA/chlorophyll/SeaWiFS/yearly/mapped"
proc_result = run(["ls", path_to_dir], capture_output=True)
file_names = proc_result.stdout.splitlines()

# Define string components to which files will be reneamed
month_suffix = ["01", "02", "03", "04", "05", "06", "07", "08", "09", 
                "10", "11", "12"]
year_label = ["1997", "1998", "1999", "2000", "2001", "2002", "2003",
              "2004", "2005", "2006", "2007", "2008", "2009", "2010"]

file_names_new = []
curr_year = 0  # First data are in 1997
curr_month = 9  # October, indexed from zero
m_flag = 1

for fname in file_names:
    # Account for case that 2008 is missing Jan - April data
    if b"2008" in fname and m_flag:
        curr_month = 4
        m_flag = 0  # prevent curr_month being set to 4 every iteration in 2008
    # Account for case that 2009 is missing May, September, and October data
    elif b"2009" in fname and curr_month == 4:
        curr_month += 1
    elif b"2009" in fname and curr_month == 8:
        curr_month += 2  # skip October and November

    # Create the new string to which the file will be renamed.
    iterable = [year_label[curr_year], month_suffix[curr_month]]
    new_fname = "_".join(iterable)
    new_fname = ".".join([new_fname, "nc"])

    # Save list of new file names for debugging.
    file_names_new.append(new_fname)

    # Rename current file.
    run(["mv", fname, new_fname])  
    
    curr_month += 1

    # Reset month counter to January, and increment year counter.
    if curr_month >= 12:
        curr_month = 0
        curr_year += 1