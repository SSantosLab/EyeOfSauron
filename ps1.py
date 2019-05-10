import os
from splinter import Browser
from selenium import webdriver

# open headless browser
# laptop
# driver=webdriver.Firefox(executable_path='/Users/jmetzger/anaconda3/bin/geckodriver')
# des machines
driver=webdriver.Firefox(executable_path='/data/des30.a/data/annis/dae-haven/py-lib/lib/python2.7/site-packages/geckodriver/geckodriver')
browser=Browser(headless=True)
# the geckodriver executable needs to be in the environmental variable $PATH
# export PATH=$PATH:/data/des30.a/data/annis/dae-haven/py-lib/lib/python2.7/site-packages/geckodriver/

# to decimal degrees
def RA_convert(RA):
    RA=RA.split(':')
    RA=[float(x) for x in RA]
    return (RA[0] + RA[1]/60. + RA[2]/3600.)*(360/24.)

# to decimal degrees
def DEC_convert(DEC):
    DEC=DEC.split(':')
    DEC=[float(x) for x in DEC]
    if DEC[0]!=0: return (abs(DEC[0]) + DEC[1]/60. + DEC[2]/3600.)*abs(DEC[0])/DEC[0]
    else: return DEC[1]/60. + DEC[2]/3600.


#
#   Get Panstarrs template image
#       save template images to new folder
#
def get_template_image(RA,DEC,size,browser,path):
    #open PS1 query
    url='http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos='+str(RA)+'%2C'+str(DEC)+\
        '&filter=i&filetypes=stack&auxiliary=data&size='+str(size)+'&output_size=0&verbose=0&autoscale=99.500000&catlist='
    browser.visit(url)
    fitsfile=browser.find_link_by_partial_text('FITS-cutout')[0]['href']
    
    newfile=fits.open(fitsfile)
    
    # write file
    new_filename=path+'RA'+str(RA)+'_DEC'+str(DEC)+'.fits'
    try: newfile.writeto(new_filename)
    except OSError:
        try:
            os.remove(new_filename)
            newfile.writeto(new_filename)
        except OSError: print('file saving error')
    return newfile[0].data
