# Outling gaps in Canopy Height Models
# authors: Dirk Pflugmacher (HU Berlin) and Kirsten Krüger (TU München)
# last revised: 20.10.2023


import os
import numpy as np
from osgeo import gdal
from skimage.measure import label
from skimage.morphology import binary_opening, disk


def gap_identify(in_file, out_file=None, mask_file=None, connectivity=2, background=0,
                 height=5, close_radius=0):
    """""
    close_radius:   Performs morphological closing operation using a circular kernel if greater than 0.
                    Size of kernel is 2 * close_radius + 1.
    connectivity:   1 for 4 neighbors and 2 for eight neighbors.
                    Two pixels are connected when they are neighbors and have the same value.
                    In 2D, they can be neighbors either in a 1- or 2-connected sense.
                    The value refers to the maximum number of orthogonal hops to consider a
                    pixel/voxel a neighbor::
    """

    if out_file is None:
        out_file = in_file[:-4] + '_patchid_cn%scr%s.tif' % (connectivity, close_radius)

    print('Mapping gaps --> %s' % out_file)

    if not os.path.isfile(in_file):
        print('File does not exist %s' % in_file)
        return

    if os.path.exists(out_file):
        print('File exists: %s' % out_file)
        return

    src_ds = gdal.Open(in_file)
    src_band = src_ds.GetRasterBand(1)
    src_img = src_band.ReadAsArray()

    # create binary image using height threshold
    bin_img = src_img <= height

    if mask_file is not None:
        if not os.path.isfile(mask_file):
            print('File does not exist %s' % in_file)
            return
        mask_ds = gdal.Open(mask_file)
        mask_img = mask_ds.ReadAsArray()
        bin_img = mask_img.astype(np.uint8) * bin_img.astype(np.uint8)
    else:
        bin_img = bin_img.astype(np.uint8)

    if close_radius > 0:
        kernel = disk(close_radius, dtype=np.uint8)
        bin_img = binary_opening(bin_img, footprint=kernel)

    result = label(bin_img, background=background, return_num=False, connectivity=connectivity)

    drv = gdal.GetDriverByName('GTiff')
    dst_ds = drv.Create(out_file, src_ds.RasterXSize, src_ds.RasterYSize, 1,
                        gdal.GDT_Int32, options=['COMPRESS=LZW'])

    dst_ds.SetGeoTransform(src_ds.GetGeoTransform())
    dst_ds.SetProjection(src_ds.GetProjectionRef())

    dst_band = dst_ds.GetRasterBand(1)
    dst_band.SetNoDataValue(0)

    dst_band.WriteArray(result, 0, 0)

    dst_ds = None
    dst_band = None

    return out_file


def sieve(src_file, out_file=None, mmu=6, neighbors=8):

    src_ds = gdal.Open(src_file)
    src_band = src_ds.GetRasterBand(1)

    if out_file is None:
        out_file = src_file[:-4] + '_mmu%sn%s.tif' % (mmu, neighbors)

    drv = gdal.GetDriverByName('GTiff')
    dst_ds = drv.Create(out_file, src_ds.RasterXSize, src_ds.RasterYSize, 1,
                        src_band.DataType, options=['COMPRESS=LZW'])

    dst_ds.SetGeoTransform(src_ds.GetGeoTransform())
    dst_ds.SetProjection(src_ds.GetProjectionRef())

    dst_band = dst_ds.GetRasterBand(1)

    gdal.SieveFilter(src_band, None, dst_band, mmu, neighbors, callback=gdal.TermProgress_nocb)

    dst_ds = None
    src_ds = None


def sieve_no_holes(src_file, out_file=None, mmu=6, neighbors=8, nodata_value=0):

    # this sieve preserves nodata pixels, i.e., it does not fill holes

    if out_file is None:
        out_file = src_file[:-4] + '_mmu%sn%s.tif' % (mmu, neighbors)

    outpath = os.path.join(os.path.dirname(out_file))

    if os.path.exists(out_file):
        print('File exists: %s' % out_file)
        return

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    # masking

    src_ds = gdal.Open(src_file)
    src_band = src_ds.GetRasterBand(1)
    src_img = src_ds.ReadAsArray()

    ##################################################
    drv = gdal.GetDriverByName('GTiff')
    dst_ds = drv.Create(out_file, src_ds.RasterXSize, src_ds.RasterYSize, 1,
                        src_band.DataType, options=['COMPRESS=LZW'])

    dst_ds.SetGeoTransform(src_ds.GetGeoTransform())
    dst_ds.SetProjection(src_ds.GetProjectionRef())

    dst_band = dst_ds.GetRasterBand(1)
    dst_band.SetNoDataValue(nodata_value)

    # we don't want undisturbed holes (pixel value = 0) to be filled in.
    # To ignore 0 pixels we need to do some trick because ignoring 0 pixels
    # by setting the mask also leaves small islands of disturbed pixels
    gdal.SieveFilter(src_band, None, dst_band, mmu, neighbors) #, callback=gdal.TermProgress_nocb)

    # last step filled-in undisturbed holes
    # remove filled-in zero holes
    dst_img = dst_band.ReadAsArray()
    dst_img = dst_img * (src_img != 0)
    dst_band.WriteArray(dst_img, 0, 0)
    # dst_band.FlushCache()

    dst_band = None
    dst_ds = None
    src_ds = None
    mem_ds = None

    # print('Created %s' % os.path.basename(out_file))

if __name__ == '__main__':

    chm_files = [r"f:\spring\berchtesgaden\lidar\2009\berchtesgaden_2009_chm_1m.tif",
                 r"f:\spring\berchtesgaden\lidar\2017\berchtesgaden_2017_chm_1m.tif",
                 r"f:\spring\berchtesgaden\lidar\2021\berchtesgaden_2021_chm_1m.tif"]

    mask_files = [r"f:\spring\berchtesgaden\lidar\2009\berchtesgaden_2009_forestmask_1m.tif",
                  r"f:\spring\berchtesgaden\lidar\2017\berchtesgaden_2017_forestmask_1m.tif",
                  r"f:\spring\berchtesgaden\lidar\2021\berchtesgaden_2021_forestmask_1m.tif"]

    for chm_file in zip(chm_files, mask_files):
        for j in (0, 2):
            patch_file = gap_identify(chm_file[0], mask_file=chm_file[1], connectivity=2, background=0, height=5,
                                      close_radius=j)
            sieve_no_holes(patch_file, mmu=400, neighbors=8)
            os.remove(patch_file)


#create gap layers with height (2,5,10), min size (100, 400, 900) 

    for chm_file in zip(chm_files, mask_files):
        for j in (0, 2):
            for h_thres in (2,5,10):
                for min_size in (100, 400, 900):
            patch_file = gap_identify(chm_file[0], mask_file=chm_file[1], connectivity=2, background=0, height=h_thres,
                                      close_radius=j)
            sieve_no_holes(patch_file, mmu=min_size, neighbors=8)
            os.remove(patch_file)