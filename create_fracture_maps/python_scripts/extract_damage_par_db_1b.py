#!/usr/bin/env python
# -*- coding: utf-8 -*-


from osgeo import gdal
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
# from skimage.feature import canny
# from PIL import Image
from scipy import ndimage
from scipy.ndimage import binary_dilation as bd
import torch
import torchvision.transforms.functional as tf
import torch.nn.functional as fct
import os
from pathlib import Path
import imageio

from scipy.ndimage import binary_dilation as bd


import sys
sys.path.insert(1, os.environ['MODULE_HOME'])
from new_unet import UNet
from odd_geo_fcts import array_to_geotiff


from multiprocessing import Pool
# os.environ['OMP_NUM_THREADS'] = '4'

# from scipy.stats import mode

sar_filepath= str(sys.argv[1])
print("processing 1b damage for: {}".format(sar_filepath))
num_procs_py = int(sys.argv[2])
print("number of available processes for damage mapping: {}".format(num_procs_py))
num_omp_threads = int(sys.argv[3])
torch.set_num_threads(num_omp_threads)
print("number of available openMP threads for torch: {}".format(torch.get_num_threads()))

device = torch.device("cuda:0" if (torch.cuda.is_available()) else "cpu")

# kernel = gkern()
kernel = 1

# model = CalveNet(1,2)
model = UNet(1,2)
model.load_state_dict(torch.load('{}/unet_26_t2_only_manual_t2_epoch_90_sd'.format(os.environ['MODELB_DICT_DIR']),
                    										map_location=torch.device(device))['model_state_dict'])
model.eval().to(device)


def gkern(kernlen=256):
    x, y = np.meshgrid(np.linspace(-1,1,256), np.linspace(-1,1,256))
    d = np.sqrt(x*x+y*y)
    sigma, mu = 1.0, 0.0
    g = np.exp(-( (d-mu)**2 / ( 2.0 * sigma**2 ) ) )
    g /= 4
    return g


def process_patch(patch, top_left_x, top_left_y):
    if np.sum(patch) != 0 and patch.shape[-1]==patch.shape[-2] and patch.shape[-1]==256:
        sar = tf.to_tensor(patch)
        sar = tf.normalize(sar,(0.5,), (0.5,))
        sar = sar.unsqueeze(0)

        with torch.no_grad():
            output_activations = model(sar).squeeze(0)
            output_activations = fct.softmax(output_activations).numpy()

            #Up for modifying this fudge
            output_activations[5:,5:] = 0
            output_activations[:-5,:-5] = 0
            
            if kernel is not None:
                output_activations *= kernel

            if np.count_nonzero(patch==0)>=10:
                patch_mask = np.zeros_like(patch)
                patch_mask[patch==0] = 1
                for bdn in range(10):
                    patch_mask=bd(patch_mask)
                patch_mask = (1-patch_mask).astype("float")
                patch_mask[patch_mask==0]=np.nan

                output_activations *= patch_mask

    else:
        output_activations = np.zeros((2, patch.shape[0], patch.shape[1]))*np.nan

    return [output_activations, top_left_x, top_left_y]


def unpack_args_proc_patch(args):
    # map only supports calling functions with one arg. We have 4 args so use this fct to expand zipped list of all our args.
    return process_patch(*args)


def process_tiff(sar_data, tile_size, stride, model, kernel):
    model.train(False)
    height, width = sar_data.shape

    patches = []
    tlxs = []
    tlys = []
    
    prob_map = np.zeros((2, height, width))

    row=1
    column=1
    patch_number=1
    while True:

        if ((column-1)*stride+tile_size)>width:
            top_left_x = width-tile_size
            reached_right = True
        else:
            reached_right = False
            top_left_x = (column-1)*stride

        if ((row-1)*stride+tile_size)>height:
            top_left_y = height-tile_size
            reached_bottom = True
        else:
            reached_bottom = False
            top_left_y = (row-1)*stride

        pch1_image = sar_data[top_left_y:(top_left_y+tile_size), top_left_x:(top_left_x+tile_size)]
        
        patches.append(pch1_image)
        tlxs.append(top_left_x)
        tlys.append(top_left_y)

        patch_number+=1
        if reached_right:
            column = 1
            row += 1
        else:
            column += 1

        if reached_right and reached_bottom:
            break


    if __name__ == '__main__':
        pool = Pool(num_procs_py)
        output_activations_and_tags = pool.map(unpack_args_proc_patch, zip(patches, tlxs, tlys))
        pool.close()
        pool.join()

    for patch_and_tag in output_activations_and_tags:
        output_activations = patch_and_tag[0]
        top_left_x = patch_and_tag[1]
        top_left_y = patch_and_tag[2]

        prob_map[:, (top_left_y):(top_left_y+tile_size), (top_left_x):(top_left_x+tile_size)] = prob_map[:, (top_left_y):(top_left_y+tile_size), (top_left_x):(top_left_x+tile_size)]+output_activations


    prob_tensor = torch.from_numpy(prob_map)
    
   # prob_map = fct.softmax(prob_tensor).numpy()
    prob_map = prob_tensor.numpy()/2
    
    
    
    prob_map = fct.softmax(prob_tensor).numpy()
    prob_map = fct.softmax(prob_tensor).numpy()
    

    segment_data = np.argmax(prob_map, 0)
    
    return segment_data, prob_map



def generate_damage(sar_tiff_filepath, outpath, seg_path, model, aoi=None, kernel=None, device="cpu"):
    date = sar_tiff_filepath.split('/')[-1].split('.')[0]

    dirpath = '/'.join(sar_tiff_filepath.split('/')[:-1])

    sar_data = gdal.Open(sar_tiff_filepath, gdal.GA_ReadOnly)

    # clip_region = aoi
    # clip_out_filepath = outpath+"/{}.damage_clip.bmp.tif".format(date)
    # if clip_region is not None:
    #     os.system("gdal_translate -projwin {} {} {} {} -of GTiff {} {}".format(clip_region[0],clip_region[1],clip_region[2],clip_region[3],sar_tiff_filepath,clip_out_filepath))
    # else:
    #     os.system("cp {} {}".format(sar_tiff_filepath, clip_out_filepath))

    # sar_data = gdal.Open(clip_out_filepath, gdal.GA_ReadOnly)


    sar_data_array = sar_data.ReadAsArray().astype(np.float32)

    damage_map, probs = process_tiff(sar_data_array, 256, 128, model.to(device), kernel)
    damage_map[np.isnan(damage_map)]=0
    probs[np.isnan(probs)]=0

    array_to_geotiff(probs[1,:,:], sar_data, outpath+"/{}.damage_probs_1b.tif".format(date), compression="DEFLATE")

    if seg_path is not None:
        mask = gdal.Open(seg_path, gdal.GA_ReadOnly).ReadAsArray()
        masked_damage = damage_map*mask
        array_to_geotiff(masked_damage, sar_data, outpath+"/{}.masked_damage.tif".format(date), compression="DEFLATE")



outpath = os.getcwd()

#IF you want to mask the damage maps with the calving front segmentations, supply a:
cf_seg_path=None

#If you want to clip to a certain region in the image, soecify aoi in APS:
aoi=None

generate_damage(sar_filepath, outpath, cf_seg_path, model, aoi=aoi, kernel=kernel, device=device)
