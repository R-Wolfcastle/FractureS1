#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 16:03:34 2021

@author: eetss
"""
from osgeo import gdal, osr
import numpy as np
from skimage.morphology import skeletonize
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from pyproj import Proj, transform
# from scipy.ndimage import binary_dilation as bd
import rasterio
import rasterio.mask
#import cartopy

def subset_geotiff(geotiff_data, im_dim):
    # data = gdal.Open(file, gdal.GA_ReadOnly)#.ReadAsArray()
    data_array = geotiff_data.ReadAsArray()
    
    # clipping to aoi
    gt = geotiff_data.GetGeoTransform()
    
    width = geotiff_data.RasterXSize
    height = geotiff_data.RasterYSize
    
    sizex = gt[1]
    sizey = gt[5]
 
    minx = gt[0]
    miny = gt[3] + width*gt[4] + height*sizey 
    maxx = gt[0] + width*sizex + height*gt[2]
    maxy = gt[3]
 
    aoi_tl_pixel_col = int(((im_dim[0]-minx)/(maxx-minx))*width)
    aoi_tl_pixel_row = int(((im_dim[1]-maxy)/(miny-maxy))*height)
    aoi_br_pixel_col = int(((im_dim[2]-minx)/(maxx-minx))*width)
    aoi_br_pixel_row = int(((im_dim[3]-maxy)/(miny-maxy))*height)
 
    # aoi_pixel_coords = [aoi_tl_pixel_col, aoi_tl_pixel_row, aoi_br_pixel_col, aoi_br_pixel_row]
    
    data_array_aoi = data_array[aoi_tl_pixel_row:aoi_br_pixel_row, aoi_tl_pixel_col:aoi_br_pixel_col]
    
    return data_array_aoi
  
def subset_geotiff_by_smaller_tiff(geotiff_data, smaller_tiff):
    # data = gdal.Open(file, gdal.GA_ReadOnly)#.ReadAsArray()
    
    gt = smaller_tiff.GetGeoTransform()
    width = smaller_tiff.RasterXSize
    height = smaller_tiff.RasterYSize
    
    sizex = gt[1]
    sizey = gt[5]
 
    minx = gt[0]
    miny = gt[3] + width*gt[4] + height*sizey 
    maxx = gt[0] + width*sizex + height*gt[2]
    maxy = gt[3]
    
    im_dim = [minx, maxy, maxx, miny]
    
    print(im_dim)
    
    out = subset_geotiff(geotiff_data, im_dim)
    
    return out
  
def coords_in_mask(mask_tif):
    mask_array = mask_tif.ReadAsArray()
    gt = mask_tif.GetGeoTransform()
    
    width = mask_tif.RasterXSize
    height = mask_tif.RasterYSize
    
    sizex = gt[1]
    sizey = gt[5]
 
    minx = gt[0]
    miny = gt[3] + width*gt[4] + height*sizey 
    maxx = gt[0] + width*sizex + height*gt[2]
    maxy = gt[3]
    
    print(minx, miny, maxx, maxy)
    
    nzs = np.nonzero(mask_array)
    
    coords = []
    
    for (row,col) in list(zip(nzs[0], nzs[1])):
      x = int(minx + (maxx-minx)*col/width)
      y = int(maxy - (maxy-miny)*row/height)
      coords.append((x,y))
      
    return coords
  
def ll_in_grid(tiff):
    mask_array = tiff.ReadAsArray()
    gt = tiff.GetGeoTransform()
    
    width = tiff.RasterXSize
    height = tiff.RasterYSize
    
    sizex = gt[1]
    sizey = gt[5]
 
    minx = gt[0]
    miny = gt[3] + width*gt[4] + height*sizey 
    maxx = gt[0] + width*sizex + height*gt[2]
    maxy = gt[3]
    
    print(minx, miny, maxx, maxy)
    
    nzs = np.nonzero(mask_array)
    
    coords = np.zeros((2, mask_array.shape[0], mask_array.shape[1]))
    
    for row in range(mask_array.shape[0]):
      for col in range(mask_array.shape[1]):
        x = int(minx + (maxx-minx)*col/width)
        y = int(maxy - (maxy-miny)*row/height)
        x, y  = aps_to_lat_lon([x,y])
        coords[0,row,col] = x
        coords[1,row,col] = y
      
    return coords
  

def array_coords_from_crs(coords, aoi, geotiff_array):
    minx = aoi[0]
    maxy = aoi[1]
    maxx = aoi[2]
    miny = aoi[3]
    
    height = geotiff_array.shape[0]
    width = geotiff_array.shape[1]
    
    rcs = []
    for coord_pair in coords:
      col = int((coord_pair[0]-minx)/(maxx-minx)*width)
      row = int((1-(coord_pair[1]-miny)/(maxy-miny))*height)
      rcs.append([row,col])
    
    return rcs
        
  
def tiff_data_at_coords(tiff, coords, zero_to_nan=False, operation=None):
    data_array = tiff.ReadAsArray()
    
    if operation == 'skeletonize':
      data_array = skeletonize(data_array)
    
    # clipping to aoi
    gt = tiff.GetGeoTransform()
    
    width = tiff.RasterXSize
    height = tiff.RasterYSize
    
    sizex = gt[1]
    sizey = gt[5]
 
    minx = gt[0]
    miny = gt[3] + width*gt[4] + height*sizey 
    maxx = gt[0] + width*sizex + height*gt[2]
    maxy = gt[3]
    
    coord_data = []
    
    for (x,y) in coords:
      pixel_col = int(((x-minx)/(maxx-minx))*width)
      pixel_row = int(((y-maxy)/(miny-maxy))*height)
      
      datum = data_array[pixel_row, pixel_col]
      
      if zero_to_nan:
        if datum==0:
          datum=np.nan
        # else:
          # print(datum, x, y, pixel_row, pixel_col)
          
      
      coord_data.append(datum)
      
    return coord_data
  
# def mean_tiff_data_in_coords(tiff, coords, zero_to_nan=False):
#     data_array = tiff.ReadAsArray()
    
#     # clipping to aoi
#     gt = tiff.GetGeoTransform()
    
#     width = tiff.RasterXSize
#     height = tiff.RasterYSize
    
#     sizex = gt[1]
#     sizey = gt[5]
 
#     minx = gt[0]
#     miny = gt[3] + width*gt[4] + height*sizey 
#     maxx = gt[0] + width*sizex + height*gt[2]
#     maxy = gt[3]
    
#     coord_data = []
    
#     for (x,y) in coords:
#       pixel_col = int(((x-minx)/(maxx-minx))*width)
#       pixel_row = int(((y-maxy)/(miny-maxy))*height)
      
#       datum = data_array[pixel_col, pixel_row]
      
#       if zero_to_nan:
#         if datum==0:
#           datum=np.nan
      
#       coord_data.append(datum)
      
#     return coord_data
  


def mask_from_coords(orgininal_tif, mask_coords):
  data_array = orgininal_tif.ReadAsArray()
  gt = orgininal_tif.GetGeoTransform()
    
  width = orgininal_tif.RasterXSize
  height = orgininal_tif.RasterYSize
  
  sizex = gt[1]
  sizey = gt[5]
 
  minx = gt[0]
  miny = gt[3] + width*gt[4] + height*sizey 
  maxx = gt[0] + width*sizex + height*gt[2]
  maxy = gt[3]
  
  mask_coords_tl = mask_coords[0]
  mask_coords_br = mask_coords[-1]
  
  mask = np.zeros_like(data_array)*np.nan
  
  row_tl = int(height*(mask_coords_tl[1]-maxy)/(miny-maxy))
  row_br = int(height*(mask_coords_br[1]-maxy)/(miny-maxy))
  
  col_tl = int(width*(mask_coords_tl[0]-minx)/(maxx-minx))
  col_br = int(width*(mask_coords_br[0]-minx)/(maxx-minx))
  
  mask[row_tl:row_br,col_tl:col_br]=1
    
  return mask

  
def array_to_geotiff(array, original_tif, new_tif_dir, compression=None):
    driver = gdal.GetDriverByName('GTiff')
    ny = array.shape[0]
    nx = array.shape[1]
    if compression:
        new_data = driver.Create(new_tif_dir, nx, ny, 1, gdal.GDT_Float32, options=['COMPRESS={}'.format(compression)])
    else:
        new_data = driver.Create(new_tif_dir, nx, ny, 1, gdal.GDT_Float32)
    geo_transform = original_tif.GetGeoTransform()  #get GeoTranform from existing dataset
    projection = original_tif.GetProjection() #similarly get from orignal tifD
    new_data.SetGeoTransform(geo_transform)
    new_data.SetProjection(projection)
    
    new_data.GetRasterBand(1).WriteArray(array)
    
    new_data.FlushCache() #write to disk
    return new_data

def array_to_multiband_geotiff(array, original_tif, new_tif_dir, compression=None):
    assert array.ndim==3, "This function only takes 3D arrays, if yours is 2D consider using array_to_geotiff"

    driver = gdal.GetDriverByName('GTiff')
    n_bands = array.shape[0]
    ny = array.shape[1]
    nx = array.shape[2]
    if compression:
        new_data = driver.Create(new_tif_dir, nx, ny, n_bands, gdal.GDT_Float32, options=['COMPRESS={}'.format(compression)])
    else:
        new_data = driver.Create(new_tif_dir, nx, ny, n_bands, gdal.GDT_Float32)
    geo_transform = original_tif.GetGeoTransform()  #get GeoTranform from existing dataset
    projection = original_tif.GetProjection() #similarly get from orignal tifD
    new_data.SetGeoTransform(geo_transform)
    new_data.SetProjection(projection)
	
    for band in range(n_bands):
        new_data.GetRasterBand(band+1).WriteArray(array[band,:,:])

    new_data.FlushCache() #write to disk
    return new_data
  
def geocode_array(array, bounding_crs_coords, resolution, filename='/Users/eetss/Documents/calving_fronts/thwaites_investigation/vulnerability.tif'):
  array = np.flip(array,0)
  drv = gdal.GetDriverByName("GTiff")
  ds = drv.Create(filename, array.shape[1], array.shape[0], 1, gdal.GDT_Float32)
  ds.SetGeoTransform([bounding_crs_coords[0], resolution, 0, bounding_crs_coords[3], 0, resolution])
  srs = osr.SpatialReference()
  srs.ImportFromEPSG(3031)
  ds.SetProjection(srs.ExportToWkt())
  # ds.SetProjection('WGS 84')
  ds.GetRasterBand(1).WriteArray(array)
 
  
def geocode_array_1(array, bounding_crs_coords, resolution, filename, compression=None):
  # array = np.flip(array,0)
  drv = gdal.GetDriverByName("GTiff")
  if compression:
    ds = drv.Create(filename, array.shape[1], array.shape[0], 1, gdal.GDT_Float32, options=['COMPRESS={}'.format(compression)])
  else:
    ds = drv.Create(filename, array.shape[1], array.shape[0], 1, gdal.GDT_Float32)
  ds.SetGeoTransform([bounding_crs_coords[0], resolution, 0, bounding_crs_coords[1], 0, resolution])
  srs = osr.SpatialReference()
  srs.ImportFromEPSG(3031)
  ds.SetProjection(srs.ExportToWkt())
  # ds.SetProjection('WGS 84')
  ds.GetRasterBand(1).WriteArray(array)
  # ds.FlushCache()
  # return ds 

 
def get_all_dates(data_dates, date_string_format='%Y%m%d', days=6):
  from_date_time = datetime.strptime(data_dates[0], date_string_format)
  to_date_time = datetime.strptime(data_dates[-1], date_string_format)
  
  all_dates = [from_date_time.strftime(date_string_format)]
  date_time = from_date_time
  while date_time < to_date_time:
      date_time += timedelta(days=6)
      all_dates.append(date_time.strftime(date_string_format))
      
  return all_dates

def dates_between(date1, date2):
  date_string_format='%Y%m%d'
  from_date_time = datetime.strptime(date1, date_string_format)
  to_date_time = datetime.strptime(date2, date_string_format)

  all_dates = [from_date_time.strftime(date_string_format)]
  date_time = from_date_time
  while date_time < to_date_time:
      date_time += timedelta(days=1)
      all_dates.append(date_time.strftime(date_string_format))

  return all_dates



def coordinate_grid(tiff):
    tiff_array = tiff.ReadAsArray()
    
    gt = tiff.GetGeoTransform()
    
    width = tiff.RasterXSize
    height = tiff.RasterYSize
    
    
    sizex = gt[1]
    sizey = gt[5]
 
    minx = gt[0]
    miny = gt[3] + width*gt[4] + height*sizey 
    maxx = gt[0] + width*sizex + height*gt[2]
    maxy = gt[3]
    
    coords = np.zeros((2, height, width))
    
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))

    
    xs = (minx + (maxx-minx)*cols/width).astype(int)
    ys = (maxy - (maxy-miny)*rows/height).astype(int)
    
    coords[0,:,:] = xs
    coords[1,:,:] = ys
    
    return coords
  

def true_north_to_grid_north_angles(coordinate_grid):
    x_coords = coordinate_grid[0,:,:]
    y_coords = coordinate_grid[1,:,:]
    
    # gives in the range [-pi, pi] with positive angle defined anticlockwise from x axis in 1st quadrant
    angles = np.arctan2(y_coords, x_coords)
    
    angles[angles<0] = -angles[angles<0]
    angles[angles>0] = (2*np.pi-angles[angles>0])
    
    return angles
  
def northing_easting_to_x_y(northing_data_tiff, easting_data_tiff):
    
  coord_grid = coordinate_grid(northing_data_tiff)
  
  angles = true_north_to_grid_north_angles(coord_grid)
  
  n_data = northing_data_tiff.ReadAsArray()
  e_data = easting_data_tiff.ReadAsArray()
  
  #rotate with 2d rotation matrix:
  
  x_data = e_data*np.cos(angles) + n_data*np.sin(angles)
  y_data = e_data*np.sin(angles) - n_data*np.cos(angles)
  
  return x_data, y_data









