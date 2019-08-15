import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import time
import os
import glob
import fitsio
import healpy as hp
from math import fabs
from astropy import wcs
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy import units as u

#======================================================
# work the list of galaxies
# 
# given a ra,dec of a search image, 
# return the list of gals and thumbs
#
def get_list_of_gals(ra,dec, gal_cat, 
        id_col_num, ra_col_num, dec_col_num, verbose=True) :
    usecols=(id_col_num, ra_col_num,dec_col_num)
    id, dc_ra, dc_dec = np.genfromtxt(
        gal_cat,unpack=True,usecols=usecols,delimiter=",", skip_header=1)
    
    distance = np.sqrt( 
        ((ra-dc_ra)*np.cos(dc_dec*2*np.pi/360.))**2 + (dec-dc_dec)**2)
    #distance = support.gc_separation(ra, dec, dc_ra, dc_dec)
    ix = distance < 1.2
    gals_id, gals_ra, gals_dec = id[ix], dc_ra[ix], dc_dec[ix]
#    for i in range(gals_ra.size):
#        if verbose: print
#            print "{} {:10d} {:10.5f} {:10.5f} \t {}".format(
#                i, gals_id[i], gals_ra[i], gals_dec[i])
    return gals_id, gals_ra, gals_dec


#======================================================
#
# the main routine, that shows both the search and template images
#
# search_dir is places like: /data/des51.b/data/DTS/src/20190510/
def work_image(search_dir, expid, image_ra,image_dec, filter, gal_list,
               id_col_num=0, ra_col_num=1, dec_col_num=2,
               do_num = -1, do_ccdnum=-1): #Stampsize is half width
    infilename = search_dir + "DECam_00{}.fits.fz".format(expid)
    do_bliss = True; 
    
    print("Opening file: " + infilename)
    hdulist = fits.open(infilename)
    
    gid,gra,gdec = get_list_of_gals(image_ra, image_dec, gal_list, 
        id_col_num, ra_col_num, dec_col_num)
    match_cutout = 15
    go_big = 25
    #match_cutout = 15
    #if do_cutout: stampSize = match_cutout/0.27
    if do_bliss:   stampSize = go_big/0.27

    count = 0
    for ccd in range(1,63):
        if (do_num != -1) and (do_num != ccd) : continue
        header = fits.open(infilename)[ccd].header #pyfits.getheader(infilename, 0)
        ccdnum = header["ccdnum"]
        if (do_ccdnum != -1) and (do_ccdnum != ccdnum) : continue

        print ccd,
        #print header

        have_read_data = False
        #Find pixel at the candidate ra dec
        w = WCS(header=header) 
        for i in range(0,gra.size) :
            # is the object on this ccd?
            corners = w.calc_footprint()
            ccd_ra_min = corners[:,0].min()
            ccd_ra_max = corners[:,0].max()
            ccd_dec_min = corners[:,1].min()
            ccd_dec_max = corners[:,1].max()
            #print ccd_ra_min, gra[i],ccd_ra_max
            #print ccd_dec_min, gdec[i], ccd_dec_max
            if ((gra[i] < ccd_ra_min) or (gra[i] > ccd_ra_max)) :  continue
            if ((gdec[i] < ccd_dec_min) or (gdec[i] > ccd_dec_max)) :  continue
            
            if not have_read_data:
                have_read_data = True
                fixed_data = get_biased_data(infilename, ccd)
                try:
                    ylen, xlen = fixed_data.shape[0], fixed_data.shape[1]
                except:
                    print "no data"
                    continue    
    
            px, py = w.all_world2pix(gra[i], gdec[i], 1)
            #objcoord = [px, py]
            #if np.isnan(px) or np.isnan(py):
            #    print "isnan"
            #    continue
            #if (px < 0 or px > 2048): continue
            #if (py < 0 or py > 4096): continue 

            #cutting image:
            X = int(px)  # cols before rot90 and flip
            Y = int(py)  # rows  ditto
            if ccdnum < 30 :
                Y = int(py+100)
            if ccdnum < 25 :
                X = int(px+50)
            if ccdnum > 54 :
                X = int(px+50)
                
            sizy = np.min(np.array([stampSize, fabs(ylen-stampSize)]))
            sizx = np.min(np.array([stampSize, fabs(xlen-stampSize)]))  
            siz = int(np.min(np.array([sizy,sizx])))
            row_min, row_max = np.int(Y-stampSize),np.int(Y+stampSize)+1
            col_min, col_max = np.int(X-stampSize),np.int(X+stampSize)+1
            if col_min < 0: col_min = 0
            if row_min < 0: row_min = 0
            if row_max > ylen: row_max = ylen
            if col_max > ylen: col_max = xlen
            print "cutout: rows,cols {}:{}".format(row_min, row_max), 
            print " {}:{}".format(col_min,col_max+1)

      
            data = np.copy(fixed_data[row_min:row_max, col_min:col_max+1])
            data = np.rot90(data)
            data = np.fliplr(data)
            search = data


            fig,axs=plt.subplots(1,2,figsize=(15,5))
            axs=axs.ravel()
            fig.suptitle('RA = '+str(round(gra[i],6))+', DEC = '+str(round(gdec[i],6)))
    
                
            if do_bliss:
                template = find_BLISS_image(gra[i],gdec[i],stampSize, filter)
                label = "Bliss single exposure"
            search_label="DECam_00{}.fits.fz ccd={}".format(expid,ccdnum)
                
            print "trying search + {}\n".format(label)
            try:
                axs[0].imshow(np.log10(template-np.amin(template)+1),origin='lower',cmap='gray')
                axs[0].set_title('template for gal {} {}'.format(i, label))
                axs[1].imshow(np.log10(search-np.amin(search)),origin='lower',cmap='gray')
                axs[1].set_title('search {}'.format( search_label))
                plt.show()
                count += 1
            except:
                print "exception, search image has shape",np.shape(search)
                pass
    print "\n Did {} galaxies!\n".format(count)


#======================================================
def find_BLISS_image(cand_ra,cand_dec,stampSize, filter):
    
    bands = ['g','r','i','z']
    cat_path = '/data/des81.a/data/luidhy/BLISS_allsky_try1/hpx/' 
    exp_path1 = '/data/des50.b/data/BLISS/'
    exp_path2 = '/data/des60.b/data/BLISS/'
    exp_path3 = '/data/des61.b/data/BLISS/'
    cat_nside = 32

    hpix = hp.ang2pix(cat_nside, cand_ra, cand_dec,lonlat=True)

    for band in bands:
        if band != filter : continue
        if hpix<10000:
            cat_file = cat_path+band+"/hpx_"+band+"_0"+str(hpix)+".fits"
        else:
            cat_file = cat_path+band+"/hpx_"+band+"_"+str(hpix)+".fits"
        if not os.path.isfile(cat_file):
            print "There is no source catalog in ", band
            exp_file = 0
            return ""
        else:
            h=fits.open(cat_file)[1].data
            print "Exposure available in ",band
            c1 = SkyCoord(cand_ra*u.deg, cand_dec*u.deg, frame='fk5')
            cat = SkyCoord(h['RA']*u.deg, h['DEC']*u.deg, frame='fk5')
            idx, d2d, d3d = c1.match_to_catalog_sky(cat)
            print "Closest object is at distance", d2d
            print "(RA,DEC)=",h['RA'][idx], h['DEC'][idx]
            print "Filter, mag, magerr", band, h['MAG_AUTO'][idx], h['MAGERR_AUTO'][idx]
            expnum = h['EXPNUM'][idx]
            ccdnum = h['CCDNUM'][idx]
            exp_fold = str(expnum)[:-2]+"00/"
            if (ccdnum<10): 
                ccdnum_str = "0"+str(ccdnum)
            else: 
                ccdnum_str = str(ccdnum)

            #Now open image

            exp_file1 = exp_path1+exp_fold+str(expnum)+"/D00"+str(expnum)+"_"+band+"_"+ccdnum_str+"_r1p1_immask.fits.fz"
            exp_file2 = exp_path2+exp_fold+str(expnum)+"/D00"+str(expnum)+"_"+band+"_"+ccdnum_str+"_r1p1_immask.fits.fz"
            exp_file3 = exp_path3+exp_fold+str(expnum)+"/D00"+str(expnum)+"_"+band+"_"+ccdnum_str+"_r1p1_immask.fits.fz"

            if os.path.isfile(exp_file1):
                exp_file = exp_file1
            elif os.path.isfile(exp_file2):
                exp_file = exp_file2
            elif os.path.isfile(exp_file3):
                exp_file = exp_file3
            else:
                print "Error: no exposure found in any path"
                exp_file = 0      

            #If image was found, make a cutout
            if (exp_file!=0):
                outfile = './cutouts/'+str(int(cand_ra))+str(int(cand_dec))+band+'.fits'
                img = cutout(exp_file,cand_ra,cand_dec,stampSize) #cand_ra,cand_dec
                img = np.rot90(img)
                img = np.fliplr(img)
      
            return img


def cutout(infilename,ra,dec,stampSize): #Stampsize is half width
    
    ptsInside2Rp = []
    print("Opening file: " + infilename)
    hdulist = fits.open(infilename)
    header = fits.open(infilename)[1].header #pyfits.getheader(infilename, 0)
    data = fits.open(infilename)[1].data
    ylen, xlen = data.shape[0], data.shape[1]
    
    #Find pixel at the candidate ra dec
    w = WCS(header=header) 
    px, py = w.all_world2pix(ra, dec, 1)
    objcoord = [px, py]
   

    X = int(px)
    Y = int(py)
    
    #cutting image:
    
    sizy = np.min(np.array([stampSize, fabs(ylen-stampSize)]))
    sizx = np.min(np.array([stampSize, fabs(xlen-stampSize)]))  
    siz = int(np.min(np.array([sizy,sizx])))
    print "cutout:",siz,px,py

    data = data[Y-siz:Y+siz+1,X-siz:X+siz+1]  

    return data


#======================================================
#
# deal with bias in decam header and in decam images
#
# this notation is [col_min, col_max, row_min, row_max]
# http://iraf.noao.edu/projects/ccdmosaic/imagedef/imagedef.html
def interp_header_keyword(header):
    section = header[1:-1]
    rows=section.split(",")[1]
    cols=section.split(",")[0]
    rmin= np.int(rows.split(":")[0])
    rmax= np.int(rows.split(":")[1])
    cmin= np.int(cols.split(":")[0])
    cmax= np.int(cols.split(":")[1])
    return rmin, rmax, cmin, cmax

def get_biased_data(infilename, num):
    
        data = fits.open(infilename)[num].data            
        header = fits.open(infilename)[num].header

        biassec = header["biasseca"]
        rmin, rmax, cmin, cmax = interp_header_keyword(biassec)
        bias = np.median(data[rmin:rmax, cmin:cmax])
        #print "bias a", bias
        trimsec =  header["trimseca"]
        rmin, rmax, cmin, cmax = interp_header_keyword(trimsec)
        data[rmin:rmax,cmin:cmax] = \
            data[rmin:rmax,cmin:cmax] - bias
            
        biassec = header["biassecb"]
        rmin, rmax, cmin, cmax = interp_header_keyword(biassec)
        bias = np.median(data[rmin:rmax, cmin:cmax])
        #print "bias b", bias
        trimsec =  header["trimsecb"]
        rmin, rmax, cmin, cmax = interp_header_keyword(trimsec)
        data[rmin:rmax,cmin:cmax] = \
            data[rmin:rmax,cmin:cmax] - bias
        
        
        #x="trimsec";print x,header[x]
        #x="datasec";print x,header[x]
        #x="trimseca";print x,header[x]
        #x="trimsecb";print x,header[x]
        #x="dataseca";print x,header[x]
        #x="datasecb";print x,header[x]
        #x="biasseca";print x,header[x]
        #x="biassecb";print x,header[x]
        #x="ccdseca";print x,header[x]
        #x="ccdsecb";print x,header[x]

        return data

#======================================================
#
# these are not used.
#
def find_search_image(image_name, cand_ra,cand_dec,stampSize, 
        dir="/data/des51.b/data/DTS/src/20190510/"):
    exp_file = dir+image_name
    img = raw_cutout(exp_file,cand_ra,cand_dec,stampSize) #cand_ra,cand_dec
    try:
        if img == -1: return -1
    except:
        img = np.rot90(img)
        img = np.fliplr(img)
      
    return img


#
#  raw decam image cutout
# 

def raw_cutout(infilename, ra,dec,stampSize, verbose=True): #Stampsize is half width
    
    ptsInside2Rp = []
    print("Opening file: " + infilename)
    hdulist = fits.open(infilename)
    
    found = False
    for ccd in range(1,63):
        if verbose: print ccd,
        header = fits.open(infilename)[ccd].header #pyfits.getheader(infilename, 0)
        #print header
        data = fits.open(infilename)[ccd].data
        try:
            ylen, xlen = data.shape[0], data.shape[1]
        except:
            if verbose: print "no data"
            continue
    
        #Find pixel at the candidate ra dec
        w = WCS(header=header) 
       
        px, py = w.all_world2pix(ra, dec, 1)
        objcoord = [px, py]
        if np.isnan(px) or np.isnan(py):
            if verbose: print "isnan"
            continue
        corners = w.calc_footprint()
        if (px < -1024 or px > 1024): continue
        if (py < -2048 or py > 2048): continue 
        print corners
        print ra,dec
        found = True
        if found: 
            print objcoord
            break
        
    if verbose: 
        if not found: 
            print "Not Found"
            return -1

    X = int(px)
    Y = int(py)
    
    #cutting image:
    
    sizy = np.min(np.array([stampSize, fabs(ylen-stampSize)]))
    sizx = np.min(np.array([stampSize, fabs(xlen-stampSize)]))  
    siz = int(np.min(np.array([sizy,sizx])))
    if verbose: print "cutout:",siz,px,py

    data = data[Y-siz:Y+siz+1,X-siz:X+siz+1]  

    return data


#======================================================
#
#
#

#======================================================
#
#
#


