"""
A simple script to reformat EMIT netCDFs to alternate formats.

Author: Philip G. Brodrick, philip.brodrick@jpl.nasa.gov
"""
#import argparse
import netCDF4
import numpy as np
from spectral.io import envi
#from emit_utils.file_checks import envi_header
import os
from osgeo import gdal
import subprocess


envi_typemap = {
    'uint8': 1,
    'int16': 2,
    'int32': 3,
    'float32': 4,
    'float64': 5,
    'complex64': 6,
    'complex128': 9,
    'uint16': 12,
    'uint32': 13,
    'int64': 14,
    'uint64': 15
}



def envi_header(inputpath):
    """
    Convert a envi binary/header path to a header, handling extensions
    Args:
        inputpath: path to envi binary file
    Returns:
        str: the header file associated with the input reference.

    """
    if os.path.splitext(inputpath)[-1] == '.img' or os.path.splitext(inputpath)[-1] == '.dat' or os.path.splitext(inputpath)[-1] == '.raw':
        # headers could be at either filename.img.hdr or filename.hdr.  Check both, return the one that exists if it
        # does, if not return the latter (new file creation presumed).
        hdrfile = os.path.splitext(inputpath)[0] + '.hdr'
        if os.path.isfile(hdrfile):
            return hdrfile
        elif os.path.isfile(inputpath + '.hdr'):
            return inputpath + '.hdr'
        return hdrfile
    elif os.path.splitext(inputpath)[-1] == '.hdr':
        return inputpath
    else:
        return inputpath + '.hdr'
        
def single_image_ortho(img_dat, glt, glt_nodata_value=0):
    """Orthorectify a single image

    Args:
        img_dat (array like): raw input image
        glt (array like): glt - 2 band 1-based indexing for output file(x, y)
        glt_nodata_value (int, optional): Value from glt to ignore. Defaults to 0.

    Returns:
        array like: orthorectified version of img_dat
    """
    outdat = np.zeros((glt.shape[0], glt.shape[1], img_dat.shape[-1]))
    valid_glt = np.all(glt != glt_nodata_value, axis=-1)
    glt[valid_glt] -= 1 # account for 1-based indexing
    outdat[valid_glt, :] = img_dat[glt[valid_glt, 1], glt[valid_glt, 0], :]
    return outdat

# parser = argparse.ArgumentParser(description="Resave EMIT NetCDF file.")
# parser.add_argument('input_netcdf', type=str, help='File to convert.')
# parser.add_argument('output_dir', type=str, help='Base directory for output ENVI files')
# parser.add_argument('-ot', '--output_type', type=str, default='ENVI', choices=['ENVI'], help='Output format')
# parser.add_argument('--interleave', type=str, default='BIL', choices=['BIL','BIP','BSQ'], help='Interleave of ENVI file to write')
# parser.add_argument('--overwrite', action='store_true', help='Overwrite existing file')
# parser.add_argument('--orthorectify', action='store_true', help='Orthorectify data')
# args = parser.parse_args()

def reformat(input_netcdf,output_dir,output_type='ENVI',interleave='BIL',\
             overwrite=False,orthorectify=False):

    nc_ds = netCDF4.Dataset(input_netcdf, 'r', format='NETCDF4')

    if os.path.isdir(output_dir) is False:
        err_str = f'Output directory {output_dir} does not exist - please create or try again'
        raise AttributeError(err_str)

    if orthorectify:
        glt = np.zeros(list(nc_ds.groups['location']['glt_x'].shape) + [2], dtype=np.int32)
        glt[...,0] = np.array(nc_ds.groups['location']['glt_x'])
        glt[...,1] = np.array(nc_ds.groups['location']['glt_y'])

    if output_type == 'ENVI':
        dataset_names = list(nc_ds.variables.keys())
        for ds in dataset_names:
            output_name = os.path.join(output_dir, os.path.splitext(os.path.basename(input_netcdf))[0] + '_' + ds)
            if os.path.isfile(output_name) and overwrite is False:
                err_str = f'File {output_name} already exists. Please use --overwrite to replace'
                raise AttributeError(err_str)
            nbands = 1
            if len(nc_ds[ds].shape) > 2:
                nbands = nc_ds[ds].shape[2]

            metadata = {
                'lines': nc_ds[ds].shape[0],
                'samples': nc_ds[ds].shape[1],
                'bands': nbands,
                'interleave': interleave,
                'header offset' : 0,
                'file type' : 'ENVI Standard',
                'data type' : envi_typemap[str(nc_ds[ds].dtype)],
                'byte order' : 0
            }

            for key in list(nc_ds.__dict__.keys()):
                if key == 'summary':
                    metadata['description'] = nc_ds.__dict__[key]
                elif key not in ['geotransform','spatial_ref' ]:
                    metadata[key] = f'{{ {nc_ds.__dict__[key]} }}'

            if orthorectify:
                metadata['lines'] = glt.shape[0]
                metadata['samples'] = glt.shape[1]
                gt = np.array(nc_ds.__dict__["geotransform"])
                metadata['map info'] = f'{{Geographic Lat/Lon, 1, 1, {gt[0]}, {gt[3]}, {gt[1]}, {gt[5]*-1},WGS-84}}'

                metadata['coordinate system string'] = f'{{ {nc_ds.__dict__["spatial_ref"]} }}' 

            if 'sensor_band_parameters' in nc_ds.__dict__.keys():
                band_parameters = nc_ds['sensor_band_parameters'].variables.keys() 
                for bp in band_parameters:
                    if bp == 'wavelengths' or bp == 'radiance_wl':
                        metadata['wavelength'] = np.array(nc_ds['sensor_band_parameters'].variables[bp]).astype(str).tolist()
                    elif bp == 'radiance_fwhm':
                        metadata['fwhm'] = np.array(nc_ds['sensor_band_parameters'].variables[bp]).astype(str).tolist()
                    elif bp == 'observation_bands':
                        metadata['band names'] = np.array(nc_ds['sensor_band_parameters'].variables[bp]).astype(str).tolist()
                    else:
                        metadata[bp] = np.array(nc_ds['sensor_band_parameters'].variables[bp]).astype(str).tolist()
            
            if 'wavelength' in list(metadata.keys()) and 'band names' not in list(metadata.keys()):
                metadata['band names'] = metadata['wavelength']

            envi_ds = envi.create_image(envi_header(output_name), metadata, ext='', force=overwrite) 
            mm = envi_ds.open_memmap(interleave='bip',writable=True)

            dat = np.array(nc_ds[ds])
            if len(dat.shape) == 2:
                dat = dat.reshape((dat.shape[0],dat.shape[1],1))

            if orthorectify:
                mm[...] = single_image_ortho(dat, glt)
            else:
                mm[...] = np.array(dat)
            del mm, envi_ds


# See gdalwarp documentation for more details on options flags:
# https://gdal.org/programs/gdalwarp.html
# -r : resampling method (near,bilinear,cubic,etc...)
# -of: output file format (geotiff, etc)
# -co: output file creation options
# -t_srs: target spatial reference (will use info from input dataset by default)
# -srcnodata: the "no data" value for the source data, to ignore areas without data
# -dstnodata: the "no data" value to record for the desitnation file in areas without data

# parser = argparse.ArgumentParser(description="Check output of obs")
# parser.add_argument('source_file', type=str)
# parser.add_argument('match_file', type=str)
# parser.add_argument('dst_file', type=str)
# parser.add_argument('-r', type=str, default='near')
# parser.add_argument('-of', type=str, default='GTiff')
# parser.add_argument('-co', type=str, default='')
# parser.add_argument('-t_srs', type=str, default='')
# parser.add_argument('-srcnodata', type=str, default='')
# parser.add_argument('-dstnodata', type=str, default='')

# args = parser.parse_args()


def match_extent(source_file,match_file,dst_file,resample_method='near',output_format='GTiff',\
                creation_options=None,target_srs=None,srcnodata=None,dstnodata=None):

    match_ds = gdal.Open(match_file,gdal.GA_ReadOnly)
    trans = match_ds.GetGeoTransform()
    x_min = trans[0]
    x_max = trans[0] + trans[1]*match_ds.RasterXSize
    y_max = trans[3]
    y_min = trans[3] + trans[5]*match_ds.RasterYSize

    cmd = f'gdalwarp {source_file} {dst_file}'
    cmd += f' -te {x_min} {y_min} {x_max} {y_max}'
    cmd += f' -tr {trans[1]} {trans[5]}'
    cmd += f' -r {resample_method}'
    cmd += f' -of {output_format}'
    if creation_options is not None:
        cmd += f' -co {creation_options}'

    if srcnodata is not None:
        cmd += f' -srcnodata {srcnodata}'

    if dstnodata is not None:
        cmd += f' -dstnodata {dstnodata}'

    if target_srs is not None:
        cmd += f' -t_srs {target_srs}'

    print(cmd)
    subprocess.call(cmd,shell=True)
    
