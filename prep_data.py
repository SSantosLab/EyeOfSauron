import numpy as np
import configparser
import os
import healpy as hp
from astropy.cosmology import Planck15 as cosmo
from astropy.io import fits
from math import erf

# read config.ini file

config=configparser.ConfigParser()

config.read('config.ini')

# set directory structure

directory=os.path.abspath(config['path']['data'])

subdirs=['galaxies','search','templates','log']

dirs=[]

for subdir in subdirs:
    d=os.path.join(directory,subdir))                           
    if not os.path.exists(d): os.makedirs(d)
    dirs.append(d)

# read skymap

print "Reading skymap for event ", event

event=config['event']['id']

d=config['path']['event_maps']

f=os.path.join(d,event+config['map']['name']+config['map']['extension'])

m=hp.read_map(f,field=range(4),verbose=False)

nside=config['catalog']['nside']

prob=hp.ud_grade(m[0],nside,power=-2)

distmu,distsigma,distnorm=hp.ud_grade(m[1:],nside)

cl=0.9

idx_sort = np.argsort(prob)
idx_sort_up = list(reversed(idx_sort))
sum = 0.
idx = 0
while prob_sum<cl:
    this_idx = idx_sort_up[idx]
    prob_sum = prob_sum+prob[this_idx]
     idx = idx+1

idx_sort_cut = idx_sort_up[:idx]

roi_area = hp.nside2pixarea(nside,degrees=True) * idx

print "ROI total area: ", roi_area, "sq-deg (npix: "+str(idx)+")"

# read galaxy cat

d=config['path']['galaxy_cat']

nsigma=erf(cl/np.sqrt(2.))

deventmin=distmu-nsigma*distsigma
deventmax=distmu+nsigma*distsigma
deventrange=deventmax-deventmin

dtype=[('ID', 'i8'), ('RA', 'f8'), ('DEC', 'f8'), ('ZPHOTO', 'f8'), 
('ZPHOTO_ERR', 'f8'), ('MAG_I', 'f8'), ('DISTMU', 'f8'), ('DISTSIGMA', 'f8'),
('PROB','f8'), ('CATALOG','S10')]

galaxies=np.array([],dtype=dtype)

for i in idx_sort_cut:
    pixnum='{:05d}'.format(i)
    filename=config['catalog']['name']+pixnum+config['catalog']['extension']
    f=os.path.join(d,config['catalog']['version'],filename)
    g=fits.open(f)[1].data
    dgalmin=cosmo.luminosity_distance(g['ZPHOT']-nsigma*g['ZPHOT_ERR'])
    dgalmax=cosmo.luminosity_distance(g['ZPHOT']+nsigma*g['ZPHOT_ERR'])
    dgalrange=dgalmax-dgalmin
    s=np.maximum(dgalmax,deventmax[i])-np.minimum(dgalmin,deventmin[i])
    overlap=g[np.where(s<dgalrange+deventrange[i])]
    overlap[6]=distmu[i]
    overlap[7]=distsigma[i]
    overlap[8]=prob[i]
    overlap[9]=pixnum
    galaxies=np.append(galaxies,overlap)

# print number of selected galaxies

# save galaxy cat file

# find templates

# find search images

# make a plot showing the sky map, selected galaxies and templates



#raw_images = os.path.abspath(config['path']['raw_images'])
#reduced_images = os.path.abspath(config['path']['reduced_images'])

