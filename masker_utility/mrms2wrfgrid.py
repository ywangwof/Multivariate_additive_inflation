# DESCRIPTION: Interpolates MRMS composite reflectivity to
# wrf grid using scipy.interpolate.RegularGridInterpolator
# with linearstod interpolation method
#
# This script is created by Dr. Yongming Wang based on MRMS2fv3grid.py, 
#      crediting to Dr. Nick Gasperoni
#
import netCDF4 as nc
import numpy as np
from scipy.interpolate import griddata
import time
import argparse

start_time = time.time()

# Create a parser object to handle command-line arguments
parser = argparse.ArgumentParser(description='''Create mask field based on observed MRMS composite reflectivity
                                 with linearstod interpolation method''')

# Add an argument for the user-defined value with a default
parser.add_argument('dbzthresh', type=float, nargs='?', default=25.0,
                    help='''threshold for masking composite reflectivity (default=35.0)''')

# Parse the command-line arguments
args      = parser.parse_args()
dbzthresh = args.dbzthresh

if dbzthresh == 25.0:
    print(f"Using the default value for dbzthresh ({dbzthresh})")
else:
    print(f"User-defined value provided for dbzthresh: {dbzthresh}")

# Replace 'your_file.nc' with the path to your NetCDF file
file_path = 'dbzobs.nc'

# Open the NetCDF file for reading
nc_file = nc.Dataset(file_path, 'r')

# Access dimensions, variables, and attributes
nlat_mrms=len(nc_file.dimensions['index'])
print(f"\nDimensions of mrms reflectivity in {file_path}:")
for dim in nc_file.dimensions:
    print(f" - {dim}: {len(nc_file.dimensions[dim])}")

lev_mrms=nc_file.variables["height"][:]
lat_mrms=nc_file.variables["lat"][:]
lon_mrms=nc_file.variables["lon"][:]
lon_mrms=np.where(lon_mrms > 180.0, lon_mrms-360.0, lon_mrms)
print(f"\nlat_mrms range: {np.amin(lat_mrms)}  {np.amax(lat_mrms)}")
print(f"lon_mrms range: {np.amin(lon_mrms)}  {np.amax(lon_mrms)}")
print(f"lev_mrms range: {np.amin(lev_mrms)}  {np.amax(lev_mrms)}")

mrms_refl3d = nc_file.variables["value"][:]
nc_file.close()

print(f"\nDimensions (nlat, nlon, nlev) of mrms_refl3d: {mrms_refl3d.shape}")
print(f"Min and Max values for mrms_refl3d: {np.amin(mrms_refl3d)}  {np.amax(mrms_refl3d)}")

# Read in gridspec file for latitude, longitude arrays
nc_grid = nc.Dataset("wrf_mask","r")
lat2d_wrfgrid = nc_grid.variables["XLAT"][:]
lon2d_wrfgrid = nc_grid.variables["XLONG"][:]
lon2d_wrfgrid = np.where(lon2d_wrfgrid > 180.0, lon2d_wrfgrid-360.0, lon2d_wrfgrid)

ph_wrfgrid    = nc_grid.variables["PH"][:]
phb_wrfgrid   = nc_grid.variables["PHB"][:]
lev_wrfgrid   = (ph_wrfgrid + phb_wrfgrid)/9.81

print(f"\nDimensions of wrf_mask for interpolation destination: {lat2d_wrfgrid.shape}")
print(f"lat2d_wrfgrid range: {np.amin(lat2d_wrfgrid)}  {np.amax(lat2d_wrfgrid)}")
print(f"lon2d_wrfgrid range: {np.amin(lon2d_wrfgrid)}  {np.amax(lon2d_wrfgrid)}")
print(f"lev_wrfgrid range: {np.amin(lev_wrfgrid)}  {np.amax(lev_wrfgrid)}")
nc_grid.close()

dims = lev_wrfgrid.shape
lon2d_wrfgrid_3D = np.empty([dims[0],dims[1],dims[2],dims[3]])
lat2d_wrfgrid_3D = np.empty([dims[0],dims[1],dims[2],dims[3]])

for iz in range(dims[1]-1):
    lon2d_wrfgrid_3D[:,iz+1,:,:] = lon2d_wrfgrid
    lat2d_wrfgrid_3D[:,iz+1,:,:] = lat2d_wrfgrid

# Interpolate mrms cref to wrf grid locations with bilinear interpolation
refl_wrfgrid = griddata((lon_mrms,lat_mrms,lev_mrms),mrms_refl3d,
                        (lon2d_wrfgrid_3D,lat2d_wrfgrid_3D,lev_wrfgrid),
                        method='linear',fill_value=0.0)
print(f"After interpolation, refl_wrfgrid dimensions: {refl_wrfgrid.shape}")
print(f"Min Max of interpolated refl_wrfgrid: {np.amin(refl_wrfgrid)} {np.amax(refl_wrfgrid)}")

refl_wrfgrid = np.where(refl_wrfgrid >= dbzthresh, 200, 0.0 )
cref_wrfgrid = np.amax(refl_wrfgrid, axis=1)


# Write to netcdf file
foutname="mask"
print(f'\nWriting to output netcdf file {foutname}')
xaxis = lon2d_wrfgrid.shape[2]
yaxis = lat2d_wrfgrid.shape[1]
zaxis = lev_wrfgrid.shape[1]-1
ncout = nc.Dataset(foutname, 'w', format='NETCDF4_CLASSIC')
ncout.createDimension('west_east', xaxis)
ncout.createDimension('south_north', yaxis)
ncout.createDimension('bottom_top',zaxis)
ncout.createDimension('Time', None)

# Create a variable for time (unlimited)
time_var = ncout.createVariable('Time', np.float32, ('Time',))

# Create a netcdf variable for output interpolated cref
cref_ncvar = ncout.createVariable('QRAIN', np.float32, ('Time','bottom_top','south_north', 'west_east'), zlib=True, complevel=4)

# Assign the data from your array to the NetCDF variable
cref_ncvar[0,:,:,:] = np.tile(cref_wrfgrid[:,np.newaxis,:,:], (1,zaxis,1,1))

# Optionally, add attributes to the variable or the NetCDF file
cref_ncvar.units = "dBZ"
ncout.description = "MRMS Composite reflectivity interpolated to WRF-ARW grid"
ncout.intmethod = "scipy.interpolate.griddata(method=linear)"

ncout.close()

print("\n--- mrms2fv3grid.py completed! Took %s seconds ---" % (time.time() - start_time))
