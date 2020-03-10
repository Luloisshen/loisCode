import numpy as np
import matplotlib.pyplot  as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

import aplpy

from astropy.coordinates import SkyCoord  # High-level coordinates
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS

#example to run two panel figure 
#fig = plt.figure(figsize=(12, 6))
#contourR(radio_file, data.loc[i, '#SAD_id'], data.loc[i, 'RA_radio'], data.loc[i, 'Dec_radio'], frame_id=(1,2,1), cutout_size=10.)
#contourO(radio_file, opt_file, data.loc[i, '#SAD_id'], data.loc[i, 'RA_radio'], data.loc[i, 'Dec_radio'], frame_id=(1,2,2), cutout_size=10., opt_vmin=-0.05, opt_vmax=0.2)
#plt.subplots_adjust(wspace=0, hspace=0)
#fig.tight_layout()

def contourR(fig, fits_radio, id, ra, dec, frame_id=(1,1,1), frame_name='radio image', cutout_size=25., radio_rms=0.00001, contour_scale=[4, 16, 64]):
    """
    plot radio contour on radio image

    input:
        fits_radio: radio image in fits
        id, ra, dec: info of target
        frame_id: use for multipanel plots
        cutout_size: the cutout size in arcmin, default is 25arcmin
        radio_rms: the rms of the radio image, in unit of the image, default is 0.00001, 
        contour_scale: scale of the contour, default is [4, 16, 64]

    How to Use:
        contourR(fits_radio, Rid, Rra, Rdec, frame_id=(1,1,1), radio_rms=rrms)
        fig.canvas.draw()
        fig.savefig('radio_contour_on_radio.pdf')

    """
    fradio = aplpy.FITSFigure(fits_radio, figure=fig, subplot=frame_id)
    fradio.show_colorscale()
    fradio.recenter(ra, dec, radius=cutout_size/3600)
    fradio.show_contour(levels=np.array(contour_scale)*radio_rms, colors='red')

    fradio.axis_labels.hide_y() #no label and tick show
    fradio.tick_labels.hide_y()
    fradio.axis_labels.hide_x()
    fradio.tick_labels.hide_x()

    fradio.add_label(0.05, 0.95, '%s\nid = %s'%(frame_name, idx), relative=True, horizontalalignment='left', verticalalignment='top') #add id in top right og the panel
    fradio.add_scalebar(5./3600, "5''", color='white', corner='top right') #add a 5 arcmin scale bar

def contourO(fig, fits_radio, fits_opt, idx, ra, dec, frame_id=(1,1,1), frame_name='optical image', cutout_size=25., radio_rms=0.00001, contour_scale=[4, 16, 64], opt_vmin=None, opt_vmax=None):
    """
    plot radio contour on optical image

    input:
        fits_radio, fits_opt: radio and optical image in fits
        id, ra, dec: info of target
        frame_id: use for multipanel plots
        cutout_size: the cutout size in arcmin, default is 25 arcmin
        radio_rms: the rms of the radio image, in unit of the image, default is 0.00001, for me it is Jy
        contour_scale: scale of the contour, add as you like, default is [4, 16, 64]
        opt_vmin=None, opt_vmax=None: setup in case aplpy fail to output great optical image

    How to Use:
        contourO(fits_radio, fits_opt, Rid, Rra, Rdec, frame_id=(1,1,1), radio_rms=rrms)
        fig.canvas.draw()
        fig.savefig('radio_contour_on_optical.pdf')

    """
    
    fopt = aplpy.FITSFigure(fits_opt, figure=fig, subplot=frame_id)
    fopt.show_grayscale(vmin=opt_vmin, vmax=opt_vmax, invert=True)
    fopt.recenter(ra, dec, radius=cutout_size/3600)
    fopt.show_contour(fits_radio, levels=np.array(contour_scale)*radio_rms, colors='red')
    fopt.axis_labels.hide_y()
    fopt.tick_labels.hide_y()
    fopt.axis_labels.hide_x()
    fopt.tick_labels.hide_x()
    
    fopt.add_label(0.05, 0.95, '%s \nid = %s'%(frame_name, idx), relative=True, horizontalalignment='left', verticalalignment='top')
    
    fopt.add_scalebar(5./3600, "5''", color='white', corner='top right') #add a 5 arcmin scale bar

