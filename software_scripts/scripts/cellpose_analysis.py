import cellpose
import geopandas as gpd
import numpy as np
import pandas as pd
import scipy
import shapely
import skimage as ski
import tifffile
import zarr

from cellpose import core, utils, io, models, metrics, plot
from collections import Counter
from matplotlib import pyplot as plt
from shapely.geometry import Polygon
from tqdm import tqdm
from matplotlib.backends.backend_pdf import PdfPages

# Allow cellpose to print stdout progress
io.logger_setup()

# Use cyto3 cellpose model
model = models.Cellpose(gpu=False, model_type='cyto3')

# Parameter to use single-image mode.
channels = np.array([[0,0]])

# data location
data_dir = '/share/workshop/Spatial_scRNA_workshop/Data/'
project_dir = "/share/workshop/Spatial_scRNA_workshop/username/"
filename = 'morphology_focus/morphology_focus_0003.ome.tif'
fileloc = f'{data_dir}Xenium_Prime_Human_Lymph_Node_Reactive_FFPE_xe_outs/{filename}'

# output dir
outdir = f'{project_dir}cellpose/'

# Levels: pixelsize in µm https://kb.10xgenomics.com/hc/en-us/articles/11636252598925
# Xenium image results are pyramidal TIFFs that store different resolution images as separate series or levels. The pixel size to micron (µm) conversion factors differ for each series and are given in the link. 

scalefactors = {
  0: 0.2125,
  1: 0.4250,
  2: 0.85,
  3: 1.7,
  4: 3.4,
  5: 6.8,
  6: 13.6,
  7: 27.2,
}


# load the stain image and visualize
level = 0
pixelsize = scalefactors[level]
a = tifffile.imread(fileloc, is_ome=False, level=level)
plt.imshow(a, cmap='binary')
plt.axis('scaled')

with PdfPages('stain_image.pdf') as pdf:
  pdf.savefig(plt.show())

pdf.close()

# subset a small test region for parameter sweeping
ystart = 16500
yend = 18500
xstart = 33050
xend = 35050

# subset the small test region and shift coordinates
subset = a[ystart:yend, xstart:xend]

## load the ground truth
gtgeoj = f'{data_dir}a.geojson'

gdf = pd.read_json(gtgeoj)
gj = gpd.read_file(gtgeoj)

numtruthmasks = len(gdf)
numtruthmasks

## Wrangle GeoDataFrame
sgeoms = []
for i,row in gj.iterrows():
  sgeoms.append(shapely.affinity.translate(row['geometry'], xoff=-xstart, yoff=-ystart))
gj['geometry'] = sgeoms

## Wrangle DataFrame
sgeoms2 = []
for i,row in gdf.iterrows():
  pcoors = row.features['geometry']['coordinates'][0]
  v = np.array([xstart, ystart])
  result = [pcoors - v[:]]
  item = {'geometry': {'coordinates': result}}
  sgeoms2.append(item)
gdf['sgeoms'] = sgeoms2


# create a mask of the polygons
mm = np.zeros(subset.shape)
for i,row in gdf.iterrows():
  p = np.array(row['sgeoms']['geometry']['coordinates']).squeeze()
  r,c = ski.draw.polygon(p[:,1],p[:,0], subset.shape)
  mm[r,c] = 1

# define downsample levels and diameters in micros to test
factors = [0,1,2,3]
diameters = [50, 75, 100] 

# create downsampled versions of the mask and stain image
dmasks = {}
for factor in factors:
  dm = scipy.ndimage.zoom(mm, pixelsize/scalefactors[factor])
  dmasks[factor] = dm

# print the image factors and shape

images = {}
for factor in factors:
  img = scipy.ndimage.zoom(subset, pixelsize/scalefactors[factor])
  images[factor] = img
  print(factor, img.shape)

# run cellpose 3 on the test area
results = {}
for diameter in diameters:
  results[diameter] = {}

  for factor,img in images.items():
  # Rescale micron diameter to match image pixels
    diam = round(diameter/(scalefactors[factor]), 0)
    print(factor, diameter, diam)

  # Run cellpose3
    masks,flows,styles,diams = model.eval(img, diameter=diam, channels=channels)
    results[diameter][factor] = [masks,flows,styles,diams]

# Visualize overlap
with PdfPages('check_overlap_image.pdf') as pdf:
  for diameter,result in results.items():
    for factor,r in result.items():
      masks,flows,styles,diams = r
      fig = plt.figure(figsize=(12,5))
      plot.show_segmentation(fig, images[factor], masks, flows[0], channels=channels[0])
      plt.tight_layout()
      plt.title(f'factor={factor} µm-diameter={diameter} ({diams})')
      pdf.savefig(plt.show())

pdf.close()

# calculate fraction of overlap and compare the number of polygons against the ground truth mask

with PdfPages('check_overlap_fraction.pdf') as pdf:
  for diameter,result in results.items():
    for factor,r in result.items():
      masks,flows,styles,diams = r
      truthmask = dmasks[factor]
      intersection_mask = np.logical_and(masks, truthmask)
  
      # Fraction overlap
      score = round(np.sum(intersection_mask)/np.sum(truthmask),3)
      # Number of polygons
      nummasks = len(Counter(masks.flatten())) - 1
      print(f'factor={factor} µm-diameter={diameter} score={score} num={nummasks}')
  
      # Configure cutoffs and plot those passing muster
      if ((score < 1.6) and (score > 1.4) and
        (nummasks <=numtruthmasks*1) and (nummasks >=numtruthmasks/1.6)):
        fig, axs = plt.subplots(1, 4, figsize=(12, 3))
        axs[0].imshow(masks)
        axs[0].set_title(f'cellpose3 mask n={nummasks}')
        axs[1].imshow(utils.masks_to_outlines(masks), cmap='binary')
        axs[1].set_title('cellpose3 mask outlines')
        axs[2].imshow(subset, cmap='binary')
        axs[2].set_title('stain image')
        axs[3].imshow(truthmask)
        axs[3].set_title(f'truth mask n={numtruthmasks}')
        plt.tight_layout()
        pdf.savefig(plt.show())

pdf.close()

# define chosen parameters
factor = 2
diameter = 75

# Run cellpose3 on the full image
img = scipy.ndimage.zoom(a, pixelsize/scalefactors[factor])
diam = round(diameter/(scalefactors[factor]), 0)
masks2,flows2,styles2,diams2 = model.eval(img, diameter=diam, channels=channels)

# check the number of cells
numcells = len(set(masks2.flatten())) - 1
numcells

# visualize
with PdfPages('cellpose_full_image.pdf') as pdf:
  plt.rcParams["figure.figsize"] = (20,20)
  plt.imshow(utils.masks_to_outlines(masks2), cmap='binary')

pdf.close()

# scale and save as NPY
resampled_mask = ski.transform.resize(masks2, a.shape, order=0, preserve_range=True, anti_aliasing=False)
resampled_mask.shape
np.save(f'{outdir}a.npy', resampled_mask, allow_pickle=True)

## the saved npy file can be used in xeniumranger to improve segmentation
#xeniumranger import-segmentation \
#             --id xr30is-a-npy \
#             --xenium-bundle /Xenium_Prime_Human_Lymph_Node_Reactive_FFPE_outs/ \
#             --cells a.npy \
#             --localcores 32 \
#             --localmem 128
