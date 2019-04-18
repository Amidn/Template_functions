import pyfits
import healpy as H
import pylab as P
from optparse import OptionParser
from matplotlib.ticker import FormatStrFormatter
import numpy as np

import matplotlib.pyplot as plt #

degree = P.pi / 180.0
 
def SetupColorBar(fig, title=None, ticks=None):
    """ Create the color bar for a HEALPix figure
    """
    for ax in fig.get_axes():
        if type(ax) is H.projaxes.HpxMollweideAxes:
            cb = fig.colorbar(ax.get_images()[0], ax=ax,
                              orientation="horizontal",
                              shrink=0.8, aspect=35,
                              pad=0.05, fraction=0.1,
                              ticks=ticks,
                              format=FormatStrFormatter("%g"))
            for label in cb.ax.get_xticklabels():
                label.set_fontsize("large")
            #cb.set_label(title, size="x-large")
            cb.set_label(title, size="large")
            ax.annotate("-180$^\circ$", xy=(1.65, 0.625), size="large")
            ax.annotate("180$^\circ$", xy=(-1.9, 0.625), size="large")

if __name__ == "__main__":
    # Set up command line options
    usage = "usage: %prog [options] INPUT.fits"
    parser = OptionParser(usage)
    parser.add_option("-T", "--title", dest="title", default="",
                      help="Plot title")
    parser.add_option("-L", "--label", dest="label", default="",
                      help="Color bar label")
    parser.add_option("-c", "--color-map", dest="cmap", default="jet",
                      help="Color map for plot")
    parser.add_option("-m", "--min", dest="min", type=float,
                      help="Plot minimum value")
    parser.add_option("-M", "--max", dest="max", type=float,
                      help="Plot maximum value")
    parser.add_option("-e", "--energy-bin", dest="ebin", type=int, default=1,
                      help="Energy bin [1..35]")
    parser.add_option("-t", "--ticks", dest="ticks", default=None,
                      help="Ticks to use in plot")
    parser.add_option("-l", "--logz", action="store_true", dest="logz",
                      default=False, help="Plot z-axis on a log scale")
    parser.add_option("-b","--batchmode", action="store_true", dest="batchMode",
                      default=False, help="Execute without interaction")
    parser.add_option("-o", "--output", dest="output", default=None,
                      help="Output image file")
    parser.add_option("-p", "--process", dest="process", default='pi0',
                      help="Process (pi0, brems, IC, PROMPT) default pi0")
    options, args = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
 
    P.rcParams.update({ "font.family" : "serif",
                        "font.size"   :  14.0 })
 
    # Read DRAGON image from FITS
    hdu = pyfits.open(args[0])
    #gammaTable = hdu[0].data[0] #case with several cases: IC, brems, pi0
    for indexProcess in range(1,len(hdu)):
        current_header = hdu[indexProcess].header
        current_process_name = current_header['PROCESS']
        print indexProcess
        print "reading HDU with ", current_process_name
        print "looking for ", options.process
        if (current_process_name == options.process):
            print "FOUND", options.process
            gammaTable = hdu[indexProcess].data

    gammaHeader = hdu[0].header
    nE=gammaHeader['NE']
    print "nE", nE
    
    Eo=P.log10(gammaHeader['EMIN'])
    deltaE =P.log10(gammaHeader['EFACTOR'])

    print "Eo ", Eo
    print "deltaE ", deltaE
    # Image from one energy bin
    image = gammaTable[options.ebin]
    #nB, nL = image.shape
    nB=360
    nL=720
    #nB=gammaHeader['NAXIS1']
    #nL=gammaHeader['NAXIS2']
 
    # Create a skymap and get data
    nside = 256
    npix = H.nside2npix(nside)
 
    map = P.zeros(npix, dtype=float)

    #Fill the map pixels
    for i in range(0, npix):
        if options.logz:
            map[i] = P.log10(image[i]) 
        else:
            map[i] = image[i] #image is assumed already to be at a given Energy
 
       # print i, b, l, ib, il

    # Set up map z-axis limits
    dMin, dMax = (map[P.argmin(map)], map[P.argmax(map)])
    if options.min:
        dMin = options.min
    if options.max:
        dMax = options.max
 
    # Set up tick labels (if any)
    ticks = []
    if options.ticks == None:
        ticks = [dMin, dMax]
    else:
        ticks = [float(t.strip()) for t in options.ticks.split(",")]
 
    # Draw the color map and color bar
    energy = 10**(Eo+(options.ebin)*deltaE)
    unit = "GeV"
    if energy > 1e3:
        energy *= 1e-3
        unit = "TeV"
 
    H.mollview(map, fig=1, coord="G",
               min=dMin,
               max=dMax,
               cbar=False,
               notext=True,
               title=": ".join([options.title, "%3.1f %s" % (energy, unit)]),
               cmap=P.get_cmap(options.cmap))
    H.graticule()

############################

    fig = P.figure(1)
    SetupColorBar(fig, title=options.label, ticks=ticks)
 
    if options.output:
        fig.savefig(options.output)
 
    if not options.batchMode:
        P.show()

############################################################################
    flux = np.zeros([nE,nB,nL],dtype=np.float32)

    for j in range(0,nE):
        print j
        image = gammaTable[j]
        
        energy = 10**(Eo+(j)*deltaE)
        print "energy in GeV", energy
        
        for i in range(0, npix):
            map[i] = image[i] #image is assumed already to be at a given Energy

        #using the cartview projection to get a 2D array with the values
        map_rot= H.cartview(map, fig=1, coord="G",
            min=dMin,
            max=dMax,
            cbar=False,
            notext=True,          
            xsize=nL, ysize=nB, #bins in Gal. lon and lat 
            title=": ".join([options.title, "%3.1f %s" % (energy, unit)]),
            cmap=P.get_cmap(options.cmap),return_projected_map=True) 

        nk,nl = map_rot.shape
        #print nk, nl
        #save the values
        for k in range(0,nk):
            for l in range(0,nl):
                flux[j][k][l]=map_rot[k][l]*1e3 #flux is in GeV s-1 cm-2 sr-1, convert GeV to MeV
#####################################################################        
        
# Create the FITS file #            
    newhdu = pyfits.PrimaryHDU(flux)
    newhdr = newhdu.header

    newhdr["NAXIS"] = (3,"data axes")
    newhdr["NAXIS1"]= (nL,"length of data in axis 1 (longitude)")
    newhdr["NAXIS2"]= (nB,"length of data in axis 2 (latitude)")
    newhdr["NAXIS3"]= (nE,"length of data in axis 3 (energy)")
    
    newhdr["CRVAL1"]= (179.75,"Start of axis 1 (longitude, degree)")
    newhdr["CDELT1"]= (-0.5,"Increment of axis 1 (longitude, degree)")
    
    newhdr["CRVAL2"]= (-89.75,"Start of axis 2 (latitude,degree)")
    newhdr["CDELT2"]= (0.5,"Increment of axis 2 (latitude,degree)")
    
    newhdr["CRVAL3"]= (3,"Start of axis 3, log(energy/MeV)") #1000 MeV
    newhdr["CDELT3"]= (0.130333768495006136,"Increment of axis 3, log(Ei+1) - log(Ei)")
    
    hduList = pyfits.HDUList([newhdu])
    hduList.writeto(str(options.process)+'_DRAGON.fits',clobber=True)
