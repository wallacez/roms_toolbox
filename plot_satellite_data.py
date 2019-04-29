# plot_satellite_data.py
# Process and plot variable of interest from satellite data (e.g. SeaWiFS, Aqua)
#
# Author Z. Wallace
# Created: 4.15.19

#from argparse import ArgumentParser
from sys import argv
from subprocess import run
from netCDF4 import Dataset
from numpy.ma import mean as mmean
#from numpy import mean as mmean
from numpy import log10, zeros

from matplotlib.pyplot import pcolormesh, show, title, colorbar, tight_layout

def generate_data(month_string, start_year, num_years):
    # Acquire list of file names.
    path_to_dir = "/Users/zach/Documents/Patagonia_Iron/data/NASA/chlorophyll/SeaWiFS/yearly/mapped"
    proc_result = run(["ls", path_to_dir], capture_output=True)
    file_names = proc_result.stdout.splitlines()
    
    file_list = []
    # Filter file list by month of interest
    for f in file_names:
        if month_string in f.decode("utf-8", "strict"):
            file_list.append("/".join([path_to_dir, f.decode("utf-8", "strict")]))
    
    # Filter file list by years of interest
    for f in file_list:
        if start_year in f:  # find index of file from which to start data collection
            findex = file_list.index(f)

    filtered_list = []
    [filtered_list.append(file_list[f]) for f in range(findex, findex+num_years)]

    # Read grid data (same across all files, so can read from any).
    var_list = Dataset(file_list[0]).variables
    lon = var_list["lon"][:]
    lat = var_list["lat"][:]

    chla_records = zeros((len(lat), len(lon), num_years))
    
    # Acquire desired number of chlorophyll time records for month of
    # interest and take the average.
    for record in range(len(filtered_list)):
        chla_records[:,:,record] = Dataset(filtered_list[record]).variables['chlor_a'][:]

    # Save data in a dictionary rather than returning and passing multiple arguments.
    chla_avg = {
        "lon" : lon,
        "lat" : lat,
        "chla" : chla_records
    }

    return chla_avg


def plot_data(av_chla, start_year, num_years, month_string):
    pcolormesh(av_chla["lon"],
               av_chla["lat"], 
               log10(mmean(av_chla["chla"],axis=2)),
               vmin=-1, vmax=0.7, cmap="jet", shading='gouraud')
    colorbar()
    if num_years != 1:
        title("SeaWiFS {} Chl-a | {} - {} Mean".format(month_string,
                                                       start_year, 
                                                       str(int(start_year)+num_years-1)))
    else:
        title("SeaWiFS {} Chl-a | {}".format(month_string,
                                             start_year))
    tight_layout()
    show()


def main():
    # arg_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
    #              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    # parser = ArgumentParser(description="Select month to plot.")
    # parser.add_argument(arg_names)
    month_string = argv[1]
    start_year = argv[2]
    num_years = int(argv[3])

    month_dict = { 
        "Jan" : "01",
        "Feb" : "02",
        "Mar" : "03",
        "Apr" : "04",
        "May" : "05",
        "Jun" : "06",
        "Jul" : "07",
        "Aug" : "08",
        "Sep" : "09",
        "Oct" : "10",
        "Nov" : "11",
        "Dec" : "12",
    }
    # Set appropriate month suffix based on user input    
    month_suffix = "_".join(["",month_dict[month_string]])
    
    # Generate data from appropriate files
    data_to_plot = generate_data(month_suffix, start_year, num_years)

    # Plot data
    plot_data(data_to_plot, start_year, num_years, month_string)

if __name__ == "__main__":
    main()