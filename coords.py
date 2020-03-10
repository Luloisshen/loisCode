import numpy as np
import pandas as pd

from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy import stats
from astropy import units as u
from astropy import cosmology
cosmo = cosmology.FlatLambdaCDM(H0=70, Om0=0.27)


def exportRegions(ra, dec, size=5, outputfile="regions.reg", addID=True, ID=[]):
    """
        export region file
        add id if addID = True
    """
    if len(ra) != len(dec):
        return "RA and Dec not match"

    f = open(outputfile, "w")
    f.write("# Region file format: DS9 version 4.0\n")
    f.write("global color=green font=\"helvetica 10 normal roman\" edit=1 move=1 delete=1 highlite=1 include=1 wcs=wcs\n")
    f.write("fk5")

    for i in range(len(ra)):
        f.write("\ncircle(" + str(ra[i]) + "," + str(dec[i]) + "," + str(size) + "\")")
        if addID:
            #text format circle(100,100,20) # text={This message has both a " and ' in it} textangle=30
            if len(ID)==0:
                ID = np.arange(len(ra))
            f.write("#text={%s}"%ID[i])
    f.close()
