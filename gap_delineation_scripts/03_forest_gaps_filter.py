# Filtering outlined gaps for gaps extent
# (removing elongated gap shapes)
# authors: Dirk Pflugmacher (HU Berlin) and Kirsten Krüger (TU München)
# last revised: 20.10.2023


import os
import numpy as np
import time
from osgeo import gdal
from skimage.morphology import binary_dilation, disk
from multiprocessing import Pool
# from matplotlib import pyplot


def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]


def combine_chunks(file_basename, out_file, nodata_value=0):

    chunk_dir = os.path.dirname(file_basename)
    file_list = [os.path.join(chunk_dir, x) for x in os.listdir(chunk_dir) if os.path.basename(file_basename) in x]
    tif_list = [x for x in file_list if x[-3:].lower() == "tif"]
    csv_list = [x for x in file_list if x[-3:].lower() == "csv"]

    patch_ds = gdal.Open(tif_list[0])
    patch_band = patch_ds.GetRasterBand(1)
    patch_img = patch_ds.ReadAsArray()

    for i in range(1, len(tif_list)):
        this_ds = gdal.Open(tif_list[i])
        this_img = this_ds.ReadAsArray()
        patch_img += this_img
        this_ds = None

    # write output dataset
    drv = gdal.GetDriverByName('GTiff')
    dst_ds = drv.Create(out_file, patch_ds.RasterXSize, patch_ds.RasterYSize, 1,
                        patch_band.DataType, options=['COMPRESS=LZW'])

    dst_ds.SetGeoTransform(patch_ds.GetGeoTransform())
    dst_ds.SetProjection(patch_ds.GetProjectionRef())
    dst_band = dst_ds.GetRasterBand(1)
    dst_band.SetNoDataValue(nodata_value)

    dst_band.WriteArray(patch_img, 0, 0)
    dst_ds = None
    patch_ds = None
    patch_band = None

    for fn in tif_list:
        os.remove(fn)

    with open(out_file[:-3] + 'csv', "w") as f_out:
        f_out.writelines('pid,count,x_extent,y_extent,buffer_height\n')
        for fn in csv_list:
            with open(fn) as infile:
                contents = infile.read()
                f_out.write(contents)
            os.remove(fn)

    os.rmdir(chunk_dir)


def gap_get_ids(patch_file, return_counts=False):

    patch_ds = gdal.Open(patch_file)
    patch_img = patch_ds.ReadAsArray()
    result = np.unique(patch_img, return_counts=return_counts)
    return result


def gap_filter(patch_file, out_file, chm_file=None, buffer_radius=20,
               buffer_min_height=10, min_extent=20, max_extent=1000, pid_range=None):

    # open patch file
    patch_ds = gdal.Open(patch_file)
    patch_band = patch_ds.GetRasterBand(1)
    patch_img = patch_ds.ReadAsArray()

    # determine pixel size
    gt = patch_ds.GetGeoTransform()
    pixel_size = gt[1]

    # open chm file
    if chm_file is not None and buffer_radius > 0:
        chm_ds = gdal.Open(chm_file)
        chm_band = chm_ds.GetRasterBand(1)
        chm_img = chm_ds.ReadAsArray()

    if pid_range is not None:
        patch_img[patch_img < pid_range[0]] = 0
        patch_img[patch_img > pid_range[1]] = 0

    # make a list of patch ids
    pids, counts = np.unique(patch_img, return_counts=True)
    # hist = np.bincount(patch_img.flatten())

    # print("Found %s patches." % len(pids))
    n = len(pids)

    csv_file = out_file[:-3] + "csv"
    f_out = open(csv_file, "w")

    # loop through patches. erase patches that do not meet requirements
    for pid in pids:

        # pid = 33135 # 2592
        this_patch_img = patch_img == pid

        # get pixel coordinates
        x, y = np.where(this_patch_img)

        # get patch extent
        x_xtent = (np.max(x) - np.min(x)) * pixel_size
        y_xtent = (np.max(y) - np.min(y)) * pixel_size

        # erase patch if minimum extent is not met
        if (x_xtent < min_extent) or (y_xtent < min_extent):
            patch_img = patch_img * (this_patch_img == 0)
            f_out.writelines('%s,%s,%s,%s,%s\n' % (pid, len(x), x_xtent, y_xtent, "NA"))
            continue

        # erase patch if maximum extent is not met
        if (x_xtent > max_extent) or (y_xtent > max_extent):
            patch_img = patch_img * (this_patch_img == 0)
            f_out.writelines('%s,%s,%s,%s,%s\n' % (pid, len(x), x_xtent, y_xtent, "NA"))
            continue
        
        # apply a minimun canop height filter surrounding the gaps, to ensure break in ecosystem
        # ! for this study, filter was not applied, as it distorts the results !
        
        if chm_file is not None and buffer_radius > 0:
            # cut to extent of patch with some buffer
            buffer_radius_px = int(round(float(buffer_radius) / pixel_size))
            minx = np.max([np.min(x) - 2 * buffer_radius_px, 0])
            miny = np.max([np.min(y) - 2 * buffer_radius_px, 0])
            maxx = np.min([np.max(x) + 2 * buffer_radius_px, this_patch_img.shape[0]])
            maxy = np.min([np.max(y) + 2 * buffer_radius_px, this_patch_img.shape[1]])

            this_patch_cut = this_patch_img[minx:maxx, miny:maxy]
            this_chm_cut = chm_img[minx:maxx, miny:maxy]

            # make buffer
            kernel = disk(buffer_radius_px, dtype=np.uint8)
            this_dil_img = binary_dilation(this_patch_cut, footprint=kernel)
            this_buf_img = this_dil_img * (this_patch_cut == 0)

            # set non-buffer pixels to NA
            tmp = np.where(this_buf_img == 0, this_chm_cut, np.nan)
            buffer_height = np.nanpercentile(tmp, 75)

            # erase patch if 75th percentile height in buffer is less than threshold
            if buffer_height < buffer_min_height:
                patch_img = patch_img * (this_patch_img == 0)

            f_out.writelines('%s,%s,%s,%s,%s\n' % (pid, len(x), x_xtent, y_xtent, buffer_height))
        else:
            f_out.writelines('%s,%s,%s,%s,%s\n' % (pid, len(x), x_xtent, y_xtent, "NA"))

    # pyplot.imshow(this_buf_img)
    # pyplot.show()
    #
    # pyplot.imshow(tmp)
    # pyplot.show()
    #
    # pyplot.imshow(this_chm_cut)
    # pyplot.show()

    # write output dataset
    drv = gdal.GetDriverByName('GTiff')
    dst_ds = drv.Create(out_file, patch_ds.RasterXSize, patch_ds.RasterYSize, 1,
                        patch_band.DataType, options=['COMPRESS=LZW'])

    dst_ds.SetGeoTransform(patch_ds.GetGeoTransform())
    dst_ds.SetProjection(patch_ds.GetProjectionRef())
    dst_band = dst_ds.GetRasterBand(1)

    dst_band.WriteArray(patch_img, 0, 0)
    dst_ds = None
    dst_band = None


def gap_filter_parallel(patch_file, out_file, chm_file=None, buffer_radius=20, buffer_min_height=10,
                        min_extent=20, max_extent=1000, cores=1, overwrite=False):

    if os.path.exists(out_file) and not overwrite:
        print("%s exists. Use overwrite keyword." % out_file)
        return

    pids = gap_get_ids(patch_file)
    pid_iter = batch(list(pids), round(np.ceil(len(pids) / cores)))

    chunk_path = os.path.join(os.path.dirname(out_file), "_tmp000")
    chunk_filebase = os.path.join(chunk_path, os.path.basename(out_file)[:-4])
    if not os.path.exists(chunk_path):
        os.mkdir(chunk_path)

    pool = Pool(processes=cores)
    for i, pids in enumerate(pid_iter):
        pid_range = (min(pids), max(pids))
        out_chunk_file = chunk_filebase + "_chunk_%s.tif" % ('{0:02d}'.format(i))
        pool.apply_async(gap_filter,
                         (patch_file, out_chunk_file),
                          kwds={'chm_file': chm_file,
                                'pid_range': pid_range,
                                'buffer_radius': buffer_radius,
                                'buffer_min_height': buffer_min_height,
                                'min_extent': min_extent,
                                'max_extent': max_extent})
    pool.close()
    pool.join()

    print(" Combining chunks..")
    combine_chunks(chunk_filebase, out_file)


if __name__ == '__main__':

    cores = 40
    buffer_radius = 0  # in meter - skip this filter if set to 0 -> We did not apply the canopy height filter in study!
    buffer_min_height = 10  # in meter
    min_extent = 20  # in meter
    max_extent = 2000  # in meter

    start = time.time()

    for year in (2021, 2017, 2009):
        for cr in (0, 2):
            for mmu in (100, 400):
                chm_file = r"f:\Spring\Berchtesgaden\lidar\%s\berchtesgaden_%s_chm_1m_artefacts_masked.tif" % (year, year)
                patch_file = r"f:\Spring\Berchtesgaden\lidar\%s\berchtesgaden_%s_chm_1m_patchid_cn2cr%s_mmu%sn8.tif" % (
                year, year, cr, mmu)
                out_file = patch_file[:-4] + "_filtered.tif"
                print("Creating %s" % out_file)
                gap_filter_parallel(patch_file, out_file, chm_file=None, buffer_radius=buffer_radius, cores=cores,
                                    buffer_min_height=buffer_min_height, min_extent=min_extent, max_extent=max_extent)

    end = time.time()
    hours, rem = divmod(end - start, 3600)
    minutes, seconds = divmod(rem, 60)
    print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds))


