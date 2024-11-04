import sys
sys.path.append('C:/Users/flapet/OneDrive - NOC/Documents/utils_python')
from functions.profiling import *
import numpy as np
import pandas as pd
import os
from tqdm import tqdm

#Set relative path
os.chdir("C:/Users/flapet/OneDrive - NOC/Documents/IDAPro/lib/db_building/")

#Paths of my datasets
doombar_path = "data/glider/OG1/Doombar_648_R.nc"
cabot_path = "data/glider/OG1/Cabot_645_R.nc"
churchill_path = "data/glider/OG1/Churchill_647_R.nc"
nelson_path = "data/glider/OG1/Nelson_646_R.nc"

paths_list = [doombar_path, cabot_path, churchill_path, nelson_path]

# Function to extract variables from netCDF to DataFrame
def from_nc_to_pd(dat, vars, prof_indexes):
    print("Building the core data...")

    # Initialize lists for columns
    prof_col = []
    time_col = []
    lon_col = []
    lat_col = []
    pres_col = []
    direction_col = []

    # Loop through each profile index
    for prof_ii in tqdm(range(len(prof_indexes)), desc="Processing profiles"):
        prof = prof_indexes[prof_ii]
        
        #Extract the vars values of this profile
        time_temp = dat['TIME'][prof].values
        lon_temp = dat['LONGITUDE_GPS'][prof].values
        lat_temp = dat['LATITUDE_GPS'][prof].values
        pres_temp = dat['PRES'][prof].values

        if any(pres_temp):
            pres_temp = interp_nan(pres_temp)
            if pres_temp[0] > pres_temp[len(pres_temp) - 1]:
                direction = 'asc'
            else:
                direction = 'desc'
        else:
            direction = 'NULL'
        
        # Repeat profile index for matching dimensions
        prof_col_temp = np.repeat(prof_ii, len(time_temp))
        direction_temp = np.repeat(direction, len(time_temp))
        
        # Append extracted data to lists
        prof_col.extend(prof_col_temp)
        time_col.extend(time_temp)
        lon_col.extend(lon_temp)
        lat_col.extend(lat_temp)
        pres_col.extend(pres_temp)
        direction_col.extend(direction_temp)

    # Create the initial DataFrame from core columns
    data_df = pd.DataFrame({
        'profile': prof_col,
        'direction': direction_col,
        'time': time_col,
        'lon': lon_col,
        'lat': lat_col,
        'pres': pres_col
    })

    # Process each variable in vars
    for var in vars:
        if var in dat.variables:
            var_array = dat[var].values
            col = []

            # Extract data for each profile
            for prof_ii in tqdm(range(len(prof_indexes)), desc=f"Processing {var}"):
                prof = prof_indexes[prof_ii]
                col_temp = var_array[prof]
                col.extend(col_temp)
            
            # Add the flattened variable data to the DataFrame
            data_df[var] = col
        else:
            print(f"{var} doesn't exist in the dataset")

    return data_df

# Test on doombar_path

for path in paths_list:

    glider_name = path[path.rfind('/') + 1 : path.find('_', path.rfind('/'))]

    print(f"Processing {glider_name} data")
    profiles = find_profiles_by_depth(path)
    glider_nc = xr.open_dataset(path)
    glider_df = from_nc_to_pd(glider_nc, [
    'TEMP',
    'CNDC',
    'CHLA',
    'BBP700',
    'DPHASE_DOXY',
    'MOLAR_DOXY',
    'TPHASE_DOXY',
    'TEMP_DOXY',
    'OXYSAT_DOXY',
    'DOWNWELLING_PAR'],
    profiles)

    glider_df['glider_name'] = glider_name

    glider_df.to_csv('data/glider/csv/' + glider_name + '_raw_profile.csv')