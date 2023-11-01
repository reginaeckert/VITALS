"""
Wrapper on gdalwarp to warp the source file to the geolocation of the match file

Author: Philip G. Brodrick, philip.brodrick@jpl.nasa.gov
"""
import argparse
from osgeo import gdal
import subprocess

# See gdalwarp documentation for more details on options flags:
# https://gdal.org/programs/gdalwarp.html
# -r : resampling method (near,bilinear,cubic,etc...)
# -of: output file format (geotiff, etc)
# -co: output file creation options
# -t_srs: target spatial reference (will use info from input dataset by default)
# -srcnodata: the "no data" value for the source data, to ignore areas without data
# -dstnodata: the "no data" value to record for the desitnation file in areas without data

parser = argparse.ArgumentParser(description="Check output of obs")
parser.add_argument('source_file', type=str)
parser.add_argument('match_file', type=str)
parser.add_argument('dst_file', type=str)
parser.add_argument('-r', type=str, default='near')
parser.add_argument('-of', type=str, default='GTiff')
parser.add_argument('-co', type=str, default='')
parser.add_argument('-t_srs', type=str, default='')
parser.add_argument('-srcnodata', type=str, default='')
parser.add_argument('-dstnodata', type=str, default='')

args = parser.parse_args()


def main(args):

    match_ds = gdal.Open(args.match_file,gdal.GA_ReadOnly)
    trans = match_ds.GetGeoTransform()
    x_min = trans[0]
    x_max = trans[0] + trans[1]*match_ds.RasterXSize
    y_max = trans[3]
    y_min = trans[3] + trans[5]*match_ds.RasterYSize

    cmd = f'gdalwarp {args.source_file} {args.dst_file}'
    cmd += f' -te {x_min} {y_min} {x_max} {y_max}'
    cmd += f' -tr {trans[1]} {trans[5]}'
    cmd += f' -r {args.r}'
    cmd += f' -of {args.of}'
    if args.co != '':
        cmd += f' -co {args.co}'

    if args.srcnodata != '':
        cmd += f' -srcnodata {args.srcnodata}'

    if args.dstnodata != '':
        cmd += f' -dstnodata {args.dstnodata}'

    if args.t_srs != '':
        cmd += f' -t_srs {args.t_srs}'

    print(cmd)
    subprocess.call(cmd,shell=True)
    
if __name__ == '__main__':
    main(args)
