import numpy as np

#
# great circle separations
#
def gc_separation(ra1, dec1, ra2, dec2) :
    delDec = dec1-dec2
    delRa = ra1-ra2
    dhav = haversine(delDec)
    rhav = haversine(delRa)
    hav = dhav + np.cos(dec1)*np.cos(dec2)*rhav
    gc_distance = ahaversine(hav)
    return gc_distance

def haversine(theta) :
    hav = np.sin(theta/2.)**2
    return hav
def ahaversine(x) :
    ahav = 2*np.arcsin(np.sqrt(x))
    return ahav
