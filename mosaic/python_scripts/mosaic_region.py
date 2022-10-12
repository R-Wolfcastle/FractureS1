#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from osgeo import gdal
import numpy as np
sys.path.insert(1, os.environ['MODULE_HOME'])
from gdal_hacks import add_merge_fct_to_vrt, add_blur_to_vrt
from pathlib import Path
from odd_geo_fcts import array_to_geotiff
import os
import subprocess
import random
#from tqdm import tqdm
import itertools
from datetime import datetime, timedelta
from multiprocessing import Pool
#import cv2 as cv
from scipy import signal
import shutil

os.system("export OMP_NUM_THREADS=1")
num_procs_py = 8


def warp_tiff(damage_fp, vx_fp, vy_fp, date:str, ref_date:str, out_tiff_fp, res):
    vx = gdal.Open(vx_fp, gdal.GA_ReadOnly).ReadAsArray()
    vy = gdal.Open(vy_fp, gdal.GA_ReadOnly).ReadAsArray()

    vx[np.isnan(vx)]=0
    vy[np.isnan(vy)]=0

    #vx = ndi.zoom(vx_big, 20)
    #vy = ndi.zoom(vy_big, 20)

    dt = days_from(date, ref_date)

    dx_real = vx*dt/(res*365.25) #metres per day times days between date and ref_date
    dy_real = -vy*dt/(res*365.25)

    #HACK ALERT!!!!!!
    dx = np.zeros((dx_real.shape[0], dx_real.shape[1]+1))
    dy = np.zeros((dx_real.shape[0], dx_real.shape[1]+1))
    dx[:,:-1] = dx_real
    dy[:,:-1] = dy_real

    damage_tiff = gdal.Open(damage_fp, gdal.GA_ReadOnly)
    damage_data = damage_tiff.ReadAsArray()

    yy, xx = np.indices(damage_data.shape)
    xmap = (xx-dx).astype(np.float32)
    ymap = (yy-dy).astype(np.float32)
    warped = cv.remap(damage_data, xmap, ymap ,cv.INTER_LINEAR)

    array_to_geotiff(warped, damage_tiff, out_tiff_fp, compression="DEFLATE")



def arbitrary_upper_round(x, base):
    return base * np.ceil(x/base)

def chunk_domain_fct(tlx, tly, brx, bry, num_imgs, chunk_size=250000):
    #500Mb of float 64 data is maximum 2000*2000 (=100000x100000m) pixels, so break roughly into that size?
    #About right for 500 images, so:
    chunk_size = round(chunk_size / np.sqrt(num_imgs/500))
    lx = brx-tlx
    ly = tly-bry
    new_lx = arbitrary_upper_round(lx, chunk_size)
    new_ly = arbitrary_upper_round(ly, chunk_size)
    print(lx, new_lx)
    print(ly, new_ly)
    n_chunks_x = int(new_lx/chunk_size)
    n_chunks_y = int(new_ly/chunk_size)
    print(n_chunks_x, n_chunks_y)
    chunk_coord_lists = []
    for i in range(n_chunks_x):
        for j in range(n_chunks_y):
            ntlx = tlx + i*chunk_size
            nbrx = ntlx + chunk_size
            ntly = tly - j*chunk_size
            nbry = ntly - chunk_size
            chunk_coord_lists.append([ntlx, ntly, nbrx, nbry])
    return chunk_coord_lists


def dates_between(date1, date2):
    """
    Get list of string dates between date1:str and date2:str.
    """
    dt1 = datetime.strptime(date1, "%Y%m%d")
    dt2 = datetime.strptime(date2, "%Y%m%d")
    delta = dt2 - dt1
    dates = [dt1 + timedelta(days=i) for i in range(delta.days + 1)]
    datesstr = [datetime.strftime(date, "%Y%m%d") for date in dates]
    return datesstr


def concatenate_txt_files(fp_list:list, outfilepath:str, superdir_prefix:str):
    with open(outfilepath, 'w') as outfile:
        for fp in fp_list:
            with open(fp) as infile:
                for line in infile:
                    outfile.write(superdir_prefix+line)

def find_files(damage_type:str, date1:str, date2:str, data_superdir:str, outdir:str,\
               overwrite_all_dates_txt:bool=True, merge_txt_files:bool=True):

    dates = dates_between(date1, date2)
    if damage_type=="1a":
        file_basename="damage_probs_1a_sqrt.tif"
    elif damage_type=="1b":
        file_basename="damage_probs_1b_filtered.tif"
    
    date_fp_txt_fps = []
    for date in dates:
        date_txt_fp = "{}/date_fp_txt_files/{}_{}_fps.txt".format(os.environ['DATA_SUPERDIR'], date, damage_type)
        if not os.path.isfile(date_txt_fp):
            print("No filepaths text file for {} ¯\_(ツ)_/¯").format(date)
        date_fp_txt_fps.append(date_txt_fp)

    if not merge_txt_files:
        return date_fp_txt_fps

    all_dates_fp_list_fp = outdir+"{}_{}_{}_fp_list.txt".format(date1, date2, damage_type)
    if os.path.isfile(all_dates_fp_list_fp) and overwrite_all_dates_txt:
        os.remove(all_dates_fp_list_fp)

    if not os.path.isfile(all_dates_fp_list_fp):
        concatenate_txt_files(date_fp_txt_fps, all_dates_fp_list_fp, data_superdir)
    
    return all_dates_fp_list_fp



def build_vrt(vrt_fp, res, ulx, uly, lrx, lry, srcnodata, vrtnodata, input_files_txt):
    subprocess.call("gdalbuildvrt {} -resolution user -tr {} {} -te {} {} {} {}\
                    -srcnodata {} -vrtnodata {} -r {} -input_file_list {}".format(
                    vrt_fp,
                    res, res, ulx,
                    lry, lrx, uly,
                    srcnodata, vrtnodata,
                    "bilinear", input_files_txt), shell=True)


def function_tiff_from_vrt(vrt_fp, res, ulx, uly, brx, bry, tiff_fp):
    gdal.SetConfigOption('GDAL_VRT_ENABLE_PYTHON', 'YES')
    print(ulx, uly, brx, bry)
    os.system("gdal_translate --config GDAL_VRT_ENABLE_PYTHON YES -ot {}\
               -co BIGTIFF=YES -co TILED=YES -co COMPRESS=DEFLATE\
               -projwin {} {} {} {} -tr {} {} {} {}".format(
               "Float64", ulx, uly,
               brx, bry, res, res,
               vrt_fp, tiff_fp))
    gdal.SetConfigOption('GDAL_VRT_ENABLE_PYTHON', None)



def process_chunk(chunk_domain, vrt_fp, tiff_fp, res, fps_txt, vrt_src_nodata, vrt_target_nodata, function="median"):
    if not os.path.isfile(tiff_fp) or os.path.getsize(tiff_fp)<1e6: #hack to check if filesize is over 1Mb
        ntlx, ntly, nbrx, nbry = chunk_domain
        if not os.path.isfile(vrt_fp):
            build_vrt(vrt_fp, res, ntlx, ntly, nbrx, nbry, vrt_src_nodata, vrt_target_nodata, fps_txt)
            if function is not None:
                add_merge_fct_to_vrt(vrt_fp, function)
        function_tiff_from_vrt(vrt_fp, res, ntlx, ntly, nbrx, nbry, tiff_fp)

def process_chunk_p_var(chunk_domain, vrt_fp, vrt_2_fp, tiff_fp, res, fps_txt):
    if not os.path.isfile(tiff_fp) or os.path.getsize(tiff_fp)<1e6: #hack to check if filesize is over 1Mb
        ntlx, ntly, nbrx, nbry = chunk_domain
        if not os.path.isfile(vrt_fp):
            build_vrt(vrt_fp, res, ntlx, ntly, nbrx, nbry, 0, 0, fps_txt)
            add_merge_fct_to_vrt(vrt_fp, "variance")
        if not os.path.isfile(vrt_2_fp):
            make_blurred_vrt(vrt_fp, vrt_2_fp)

        function_tiff_from_vrt(vrt_2_fp, res, ntlx, ntly, nbrx, nbry, tiff_fp)


def unpack_args_and_process_chunk(args):
    process_chunk(*args)

def unpack_args_and_process_chunk_p_var(args):
    process_chunk_p_var(*args)

def zip_up_args(chunk_domains, vrt_fps, tiff_fps, res, fps_txt, vrt_src_nodata, vrt_target_nodata, function="median"):
    return zip(chunk_domains, vrt_fps, tiff_fps, itertools.repeat(res), itertools.repeat(fps_txt),
               itertools.repeat(vrt_src_nodata),
               itertools.repeat(vrt_target_nodata), itertools.repeat(function))

def zip_up_args_p_var(chunk_domains, vrt_fps, vrt_2_fps, tiff_fps, res, fps_txt):
    return zip(chunk_domains, vrt_fps, vrt_2_fps, tiff_fps, itertools.repeat(res), itertools.repeat(fps_txt))

def merge_chunks(chunk_fps, out_fp):
    os.system("gdal_merge.py -o {} -co TILED=YES -co COMPRESS=DEFLATE\
               -n {} -ot Float32 -tap {}".format(
               out_fp, 0, " ".join(chunk_fps)))

def merge_chunks_fp_txt(chunk_fps_txt, out_fp):
    os.system("gdal_merge.py -o {} -co TILED=YES -co COMPRESS=DEFLATE\
               -n {} -ot Float32 -tap -file {}".format(
               out_fp, 0, chunk_fps_txt))

def merge_chunks_warp(chunk_fps, out_fp):
    os.system("gdalwarp {} {}".format(out_fp, 0, chunk_fps_txt))


def mosaic_1a(date1, date2, data_superdir, outdir, aoi, aoi_name, res, chunk_domain=False, merge_function=None, lagcor=False):
    ulx, uly, lrx, lry = aoi
    
    if lagcor:
        fps_txt_fps = find_files(damage_type, date1, date2, data_superdir=data_superdir, outdir=outdir, merge_txt_files=False)

        for fp_fp in fps_txt_fps:
            pass        
    else:
        fps_txt_fp = find_files(damage_type, date1, date2, data_superdir=data_superdir, outdir=outdir, merge_txt_files=True)
        
        if chunk_domain:
            num_imgs = 0
            with open(fps_txt_fp) as f:
                for line in f:
                    num_imgs += 1
            
            if num_imgs == 0:
                return None
            
            chunk_domains = chunk_domain_fct(ulx, uly, lrx, lry, num_imgs)

            vrt_fps = [outdir+"{}_{}_{}_1a_{}_{}m_dt{}.vrt".format(date1, date2, aoi_name, merge_function, res, domain_tag) for domain_tag in list(np.arange(len(chunk_domains)))]
            tiff_fps = [outdir+"{}_{}_{}_1a_{}_{}m_dt{}.tiff".format(date1, date2, aoi_name, merge_function, res, domain_tag) for domain_tag in list(np.arange(len(chunk_domains)))]
        
            argzip = zip_up_args(chunk_domains, vrt_fps, tiff_fps, res, fps_txt_fp, 0, np.nan, merge_function)
        
            if __name__ == '__main__':
                pool = Pool(num_procs_py)
                pool.map(unpack_args_and_process_chunk, argzip)
                pool.close()
                pool.join()
        
            completed_tiff_fps = [fp for fp in tiff_fps if os.path.isfile(fp)]
        
            merged_tiff_outpath = outdir+"{}_{}_{}_1a_{}_{}m.tiff".format(date1, date2, aoi_name, merge_function, res)
            merge_chunks(completed_tiff_fps, merged_tiff_outpath)

        else:

            vrt_fp = outdir+"{}_{}_{}_1a_{}_{}m.vrt".format(date1, date2, aoi_name, merge_function, res)
            
            print("building vrt file: {}".format(vrt_fp))
            build_vrt(vrt_fp, res, ulx, uly, lrx, lry, np.nan, np.nan, fps_txt_fp)
            print("adding fct: {} to vrt_file: {}".format(merge_function, vrt_fp))
            if merge_function is not None:
                add_merge_fct_to_vrt(vrt_fp, merge_function)

            merged_tiff_outpath = outdir+"{}_{}_{}_1a_{}_{}m.tiff".format(date1, date2, aoi_name, merge_function, res)
            print("translating {} to {}".format(vrt_fp, tiff_fp))
            function_tiff_from_vrt(vrt_fp, res, ulx, uly, lrx, lry, merged_tiff_outpath)

    return merged_tiff_outpath

def gkern_coefs(kernlen=21, std=5):
    gkern1d = signal.gaussian(kernlen, std=std).reshape(kernlen, 1)
    gkern2d = np.outer(gkern1d, gkern1d)
    return gkern2d

def make_blurred_vrt(vrt_fp, vrt_2_fp):
    #gdal.SetConfigOption('GDAL_VRT_ENABLE_PYTHON', 'YES')
    if not os.path.isfile(vrt_2_fp):
        subprocess.call("gdalbuildvrt {} {}".format(vrt_2_fp, vrt_fp), shell=True)

    coef_list = gkern_coefs(35, 10)
    add_blur_to_vrt(vrt_2_fp, coef_list, 35)

    #os.system("gdal_translate --config GDAL_VRT_ENABLE_PYTHON YES -co TILED=YES {} {}".format(vrt_2_fp, outdir+"bs_filt_var_blur_2021.tif"))
    #gdal.SetConfigOption('GDAL_VRT_ENABLE_PYTHON', None)


def mosaic_1b(date1, date2, data_superdir, outdir, aoi, aoi_name, res, chunk_domain=False, merge_function=None, p_var_calculation=False):
    ulx, uly, lrx, lry = aoi

    fps_txt_fp = find_files("1b", date1, date2, data_superdir=data_superdir, outdir=outdir, merge_txt_files=True)
   
    if p_var_calculation:
        merge_function = "variance"

        num_imgs = 0
        with open(fps_txt_fp) as f:
            for line in f:
                num_imgs += 1 

        if num_imgs == 0:
            return None

        chunk_domains = chunk_domain_fct(ulx, uly, lrx, lry, num_imgs)

        vrt_fps = [outdir+"{}_{}_{}_1b_{}_{}m_dt{}.vrt".format(date1, date2, aoi_name, "variance", res, domain_tag) for domain_tag in list(np.arange(len(chunk_domains)))]
        vrt_2_fps = [outdir+"{}_{}_{}_1b_{}_{}m_dt{}.vrt".format(date1, date2, aoi_name, "blur", res, domain_tag) for domain_tag in list(np.arange(len(chunk_domains)))]
        tiff_fps = [outdir+"{}_{}_{}_1b_{}_{}m_dt{}.tiff".format(date1, date2, aoi_name, "p_var", res, domain_tag) for domain_tag in list(np.arange(len(chunk_domains)))]

        argzip = zip_up_args_p_var(chunk_domains, vrt_fps, vrt_2_fps, tiff_fps, res, fps_txt_fp)
        
        #if __name__ == '__main__':
        #    pool = Pool(num_procs_py)
        #    pool.map(unpack_args_and_process_chunk_p_var, argzip)
        #    pool.close()
        #    pool.join()

        completed_tiff_fps = [fp for fp in tiff_fps if os.path.isfile(fp)]
        #ctfp_fp_txt = outdir+"{}_{}_{}_1b_{}_{}m_domain_fps.txt"
        #with open(ctfp_fp_txt, 'w') as file_:
        #for ctfp in completed_tiff_fps:
        #    file_.write(ctfp)

        print("MERGING {} TIFFS".format(len(completed_tiff_fps)))

        merged_tiff_outpath = outdir+"{}_{}_{}_1b_{}_{}m.tiff".format(date1, date2, aoi_name, "p_var", res)
        merge_chunks(completed_tiff_fps, merged_tiff_outpath)


#        vrt_fp = outdir+"{}_{}_{}_1b_{}_{}m.vrt".format(date1, date2, aoi_name, merge_function, res)
#
##        make_vrt(vrt_fp, ulx, uly, lrx, lry, res, fps_txt_fp, 0, 0, "variance")
#
#        if not os.path.isfile(vrt_fp):
#            build_vrt(vrt_fp, res, ulx, uly, lrx, lry, 0, 0, fps_txt_fp)
#            add_merge_fct_to_vrt(vrt_fp, "variance")
#
#        vrt_2_fp = outdir+"{}_{}_{}_1b_{}_{}m.vrt".format(date1, date2, aoi_name, "p_var", res)
#        #vrt_2_fp = outdir+"{}_var_blur_{}.vrt".format(damage_type, time_)
#       
#        if not os.path.isfile(vrt_2_fp):
#            make_blurred_vrt(vrt_fp, vrt_2_fp)
#        
#        num_imgs = 0
#        with open(fps_txt_fp) as f:
#            for line in f:
#                num_imgs += 1 
#
#        chunk_domains = chunk_domain_fct(ulx, uly, lrx, lry, num_imgs)
#        
#        vrt_fps = [vrt_2_fp for i in range(num_imgs)]
#        tiff_fps = [outdir+"{}_{}_{}_1b_{}_{}m_dt{}.tiff".format(date1, date2, aoi_name, "p_var", res, domain_tag) for domain_tag in list(np.arange(len(chunk_domains)))]
#
#        #It'll sort it out even though I've given it None as the fp input, as I've already made the vrt.
#        argzip = zip_up_args(chunk_domains, vrt_fps, tiff_fps, res, None, 0, np.nan, function=None)
#        
#        if __name__ == '__main__':
#            pool = Pool(num_procs_py)
#            pool.map(unpack_args_and_process_chunk, argzip)
#            pool.close()
#            pool.join()
#
#        completed_tiff_fps = [fp for fp in tiff_fps if os.path.isfile(fp)]



    elif chunk_domain:
        num_imgs = 0
        with open(fps_txt_fp) as f:
            for line in f:
                num_imgs += 1

        chunk_domains = chunk_domain_fct(ulx, uly, lrx, lry, num_imgs)

        vrt_fps = [outdir+"{}_{}_{}_1b_{}_{}m_dt{}.vrt".format(date1, date2, aoi_name, merge_function, res, domain_tag) for domain_tag in list(np.arange(len(chunk_domains)))]
        tiff_fps = [outdir+"{}_{}_{}_1b_{}_{}m_dt{}.tiff".format(date1, date2, aoi_name, merge_function, res, domain_tag) for domain_tag in list(np.arange(len(chunk_domains)))]

        argzip = zip_up_args(chunk_domains, vrt_fps, tiff_fps, res, fps_txt_fp, 0, np.nan, merge_function)

        if __name__ == '__main__':
            pool = Pool(num_procs_py)
            pool.map(unpack_args_and_process_chunk, argzip)
            pool.close()
            pool.join()

        completed_tiff_fps = [fp for fp in tiff_fps if os.path.isfile(fp)]

        merged_tiff_outpath = outdir+"{}_{}_{}_1b_{}_{}m.tiff".format(date1, date2, aoi_name, merge_function, res)
        merge_chunks(completed_tiff_fps, merged_tiff_outpath)
    else:
        pass

    return merged_tiff_outpath

def mosaic(damage_type:str, date1:str, date2:str,\
           data_superdir:str, outdir:str, savedir:str,\
           aoi:list, aoi_name:str, res:int,\
           chunk_domain:bool=False, merge_function=None,\
           lagcor:bool=False, p_var_calulcation:bool=False):

    Path(outdir).mkdir(parents=True, exist_ok=True)
    Path(savedir).mkdir(parents=True, exist_ok=True)

    if damage_type not in {"1a","1b"}:
        raise ValueError('Argument `damage_type` should be in {1a, 1b}')
    if merge_function not in {None}.union({"max", "mean", "median"}):
        raise ValueError('Argument `merge_function` should be in {None, `mean`, `median`, `max`}')

    if not os.path.isfile(savedir+"{}_{}_{}_{}_mosaic.tiff".format(date1, res, aoi_name, merge_function)):

        if damage_type=="1a":
            outfp = mosaic_1a(date1, date2, data_superdir, outdir, aoi, aoi_name, res, chunk_domain, merge_function, lagcor)
        else:
            outfp = mosaic_1b(date1, date2, data_superdir, outdir, aoi, aoi_name, res, chunk_domain, merge_function, p_var_calculation)
    
        if outfp is not None:
            shutil.copy(outfp, savedir+"{}_{}_{}_{}_mosaic.tiff".format(date1, res, aoi_name, merge_function))
            os.remove(outfp)
            files = glob.glob(outpath+"/*")
            for f in files:
                os.remove(f)

    

damage_type="1a"
date1 = str(sys.argv[1])
date2 = date1
data_superdir = os.environ['DATA_SUPERDIR']
merge_function = "median"
lagcor = False
savedir = "{}/{}/".format(os.environ['DATA_SUPERDIR'], date1)
if os.environ['TEMP_STORAGE']:
    outdir = "{}/daily_mosaics/1a/{}/".format(os.environ['TEMP_STORAGE'], date1)
else:
    outdir = savedir.copy()
aoi = [-2615000, 2251000, 3003000, -2200000] ##Antarctica
aoi_name = "AIS"
res = 250
chunk_domain=True
p_var_calculation=False

mosaic(damage_type, date1, date2, data_superdir, outdir, savedir, aoi, aoi_name, res, chunk_domain, merge_function, lagcor, p_var_calculation)






