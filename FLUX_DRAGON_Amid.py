#!/bin/bash/python

import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import numpy as np
import healpy as hp
import pyfits
import os
import math
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import interp1d

# ------------------------  EDIT HERE --------------------------------------------------------------------------------------------------------
#windows_l = [0., 10., 20., 30., 40., 50., 60., 90., 120., 150., 180.] #Galactic longitudes of each window; latitude is from -5 to 5
#windows_l = [10., 11., 12., 13., 14., 15., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 50.] #Galactic longitudes of each window; latitude is from -5 to 5
#windows_l = [55.,56.,57.,58.,59., 60., 61., 62., 63., 64., 65.]
min_En = 0.3e3
max_En = 100e3
size_window= 0.5

windows_l = []
limit_lon_min= 0.
limit_lon_max= 85.5

step_window=1.

#TEMP
window_range=(limit_lon_max-limit_lon_min)/step_window

for i in xrange(0,int(window_range)+1):
  windows_l.append(i*step_window+limit_lon_min)

print 'windows_l', windows_l

#windows_b = [0., 5.] #+-5
windows_b = [ -2. , 2. ] #+-2
# --------------------------------------------------------------------------------------------------------------------------------------------

#filename = ['GammaModelPlus_FerriereSource_ringModel_res8_fullSky.fits']
filename = ['BaseModelPlus_FerriereSource_ringModel_res8_fullSky.fits'] #base

stringTitle = "Amid"

norm = [1.] #0.2*1.5
colors = ["red"]

nModels = len(norm)

#f = open('GammaModel.txt', 'w')
#fE = open('Espectrum_GammaModel_Evec.txt', 'w')
#fP = open('Espectrum_GammaModel-pi0.txt', 'w')
#fT = open('Espectrum_GammaModel-tot.txt', 'w')

f = open('BaseModel.txt', 'w')
#fE = open('Espectrum_BaseModel_Evec.txt', 'w')
#fP = open('Espectrum_BaseModel-pi0.txt', 'w')
#fT = open('Espectrum_BaseModel-tot.txt', 'w')
# --------------------------------------------------------------------------------------------------------------------------------------------

#####################################################

def pixnumrectangle_symmetric(nside,lowb_North,highb_North,lowb_South,highb_South,lowl,highl): #returns a mask with the pixels inside the ROI
  
  # moving to radiants
  lowb_North = lowb_North/180.*np.pi
  lowb_South = lowb_South/180.*np.pi
  lowl= lowl/180.*np.pi
  highb_North = highb_North/180.*np.pi
  highb_South = highb_South/180.*np.pi
  highl= highl/180.*np.pi

  npix=12*nside**2
  #print "npix", npix

  listpix = np.arange(npix)
  theta,phi = hp.pixelfunc.pix2ang(nside,listpix)
  b = np.pi/2.-theta
  l = phi
  mask = []
  
  for i in np.arange(npix):
    
    if(l[i]>np.pi):
      l[i]-=2.*np.pi
  
    if((((b[i] >= lowb_North) and (b[i] <= highb_North) ) or ( (b[i] >= lowb_South) and (b[i] <= highb_South))) and ((l[i] >= lowl) and (l[i] <= highl) or (l[i] >= -highl) and (l[i] <= -lowl))):
      mask.append(1)
    else:
      mask.append(0)

  #print mask
  return mask

#####################################################

def find_map_process(filename, process, npix, dimE):
  print "*** process: ", process, "; reading map..."
  map_ = np.zeros((npix,dimE))
  #print map_
  print map_.shape
  print "hdulist length = ", len(hdulist)
  for indexProcess in range(1,len(hdulist)):
     current_header = hdulist[indexProcess].header
     current_process_name = current_header['PROCESS']
     print indexProcess
     print "reading HDU with ", current_process_name
     print "looking for ", process
     if (current_process_name == process):
       print "FOUND", process
       current_hdu = hdulist[indexProcess].data
#      map_ = current_hdu	
       for i in range(0,npix):
         for ie in range(0,dimE):
           map_[i][ie] = current_hdu[ie][i]
  print "done"
  return map_

############### PACO to find nearest ################
#def find_nearest(array,value):
#    idx = (np.abs(array-value)).argmin()
#    return array[idx]
#####################################################
def find_nearest(array, value):
    ''' Find nearest value in an array '''
    idx = (np.abs(array-value)).argmin()
    return idx
#####################################################  

radius=math.radians(2.)#TEMP
solid_angle= 2*np.pi*(1-math.cos(radius))

for ind in xrange(len(windows_l)): #paco -1 removed

    ###begin preamble

    rc('text', usetex=True)
    rc('font', family='serif')
    rc('font', serif='Helvetica Neue')
    rc('xtick', labelsize=18)
    rc('ytick', labelsize=18)
    rcParams['legend.numpoints'] = 1
    rcParams['lines.linewidth'] = 3
    rcParams['figure.autolayout'] = True

    fig = plt.figure(figsize=(10., 7.))
    ax = fig.add_subplot(1, 1, 1)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.5)

    ax.minorticks_on()
    ax.tick_params('both', length=15, width=1.5, which='major', pad=6)
    ax.tick_params('both', length=10, width=1.3, which='minor', pad=6)

    plt.xticks(size=15)
    plt.yticks(size=15)

    ###end preamble

    bmin_south = windows_b[0]
    bmax_south = windows_b[1]
    bmin_north = windows_b[0]
    bmax_north = windows_b[1]

#    bmin_south = -windows_b[1]
#   bmax_south = -windows_b[0]
#   bmin_north = windows_b[0]
#   bmax_north = windows_b[1]

    print
    print
    # --------------------------------------------------------------------------------------------------------------------------------------------
    print "Reading models..."
    # --------------------------------------------------------------------------------------------------------------------------------------------
    print
    print

    for iModel in xrange(nModels):

        hdulist = pyfits.open(filename[iModel])
        n_ext  = len(hdulist)

        first_hdu = hdulist[0]
        first_hdu_header = first_hdu.header
        nside = first_hdu_header['nside']
        npix = 12*nside**2
        dimE = first_hdu_header['nE']
        Emin = first_hdu_header['Emin']
        Emax = first_hdu_header['Emax']
        E_factor = first_hdu_header['Efactor']

        E_vec = np.zeros(dimE)
        for i in xrange(0,dimE):
          E_vec[i] = Emin*pow(E_factor,i)

        print "npix model = ", npix
        print "dimE = ", dimE
        print "npix*dimE =", npix*dimE

        map_brems = np.zeros((npix,dimE))
        map_brems[:,:] = find_map_process(hdulist, 'brems', npix, dimE)

        map_IC = np.zeros((npix,dimE))
        map_IC[:,:] = find_map_process(hdulist, 'IC', npix, dimE)

        print map_IC

        map_pi0 = np.zeros((npix,dimE))
        map_pi0[:,:] = find_map_process(hdulist, 'pi0', npix, dimE)

        hdulist.close

        # --------------------------------------------------------------------------------------------------------------------------------------------
        print "Integrating model ", iModel, " over sky window..."
        # --------------------------------------------------------------------------------------------------------------------------------------------

        current_mask = pixnumrectangle_symmetric(nside, bmin_north, bmax_north, bmin_south, bmax_south, windows_l[ind]-size_window, windows_l[ind]+size_window)
        print "current mask -----------------: ", nside, bmin_north, bmax_north, bmin_south, bmax_south, windows_l[ind]-size_window, windows_l[ind]+size_window #current_mask






        current_spectrum_brems  = np.zeros(dimE)
        current_spectrum_IC = np.zeros(dimE)
        current_spectrum_pi0 = np.zeros(dimE)
        current_spectrum_tot = np.zeros(dimE)

        for ip in xrange(0,dimE-1):

          print "--> energy bin number = ", ip

          npix_brems = 0

          for ipix in np.arange(npix):

            if(current_mask[ipix]==1):
                npix_brems += 1
                current_spectrum_brems[ip] += (norm[iModel]*map_brems[ipix,ip]/(E_vec[ip]**2.))
            if(current_mask[ipix]==1):
                current_spectrum_IC[ip] += (map_IC[ipix,ip]/(E_vec[ip]**2.))
            if(current_mask[ipix]==1):
                current_spectrum_pi0[ip] += (norm[iModel]*map_pi0[ipix,ip]/(E_vec[ip]**2.))

          print "pixels: ", npix_brems
          if (npix_brems > 0):
            current_spectrum_brems[ip] /= npix_brems
            current_spectrum_IC[ip] /= npix_brems
            current_spectrum_pi0[ip] /= npix_brems

        current_spectrum_tot = current_spectrum_brems/(E_vec**2.) + current_spectrum_IC/(E_vec**2.) + current_spectrum_pi0/(E_vec**2.) 
        print current_spectrum_tot
        print "longitude", windows_l[ind]
        print "E_vec"
        print E_vec
        print "IC spectrum"
        print current_spectrum_IC #*solid_angle
        print "pi0 spectrum"
        print current_spectrum_pi0 #*solid_angle
        print "brems spectrum"
        print current_spectrum_brems #*solid_angle
        print "tot spectrum"
        print (current_spectrum_brems+current_spectrum_IC+current_spectrum_pi0) #*solid_angle

        total = (current_spectrum_brems+current_spectrum_IC+current_spectrum_pi0) #*solid_angle
        #select the total value at ~7TeV
# value = 7.e3 #7TeV
#     nearest = find_nearest(E_vec, value)
#       print "nearest", nearest
#       print "value nearest", total[nearest]

#       value_interp=np.interp(value, E_vec, total)
#       value_IC=np.interp(value, E_vec, current_spectrum_IC)
#       value_pi0=np.interp(value, E_vec, current_spectrum_pi0)
#       value_brems=np.interp(value, E_vec, current_spectrum_brems)
#       print "value interpolated", value_interp
#       print "Pi0",value_pi0
#       print "IC", value_IC
#       print "brem",value_brems
#       print "total", value_pi0 + value_brems + value_IC
        
        E_start=0
        E_finish = 40
        for n in range (0, 40):
            if E_vec[n] < min_En:
            #print E_vec[n]
                E_start += 1
            if E_vec[n]>max_En:
               # print "---", E_vec[n]
                E_finish += -1
        print "_____________________", E_start, E_finish

        spec_pi0 = 0
        for E in range (E_start, E_finish):
            spec_pi0 = (E_vec[E + 1] - E_vec[E]) * ( current_spectrum_pi0[E] + current_spectrum_pi0[E + 1] ) / 2.0
            spec_pi0 += spec_pi0

        spec_tot = 0
        for E in range (E_start, E_finish):
            spec_tot = (E_vec[E + 1] - E_vec[E]) * ( total[E] + total[E + 1] ) / 2.0
            spec_tot += spec_tot
        
        # Spectrum_total =  [x + y + z for x, y, z  in zip(current_spectrum_IC + current_spectrum_pi0 + current_spectrum_brems)]
        
        f.write("%f\t"   % windows_l[ind])
        f.write("%.3e\t" % E_start)
        f.write("%.3e\t" % E_finish)
        f.write("%.3e\t" % spec_tot)
        f.write("%.3e\n" % spec_pi0)
        
        # fE.write("%f\t"   % windows_l[ind])
        #for j in range (len(E_vec)):
        #   fE.write("%.3e\t" %E_vec[j])
        #fE.write("\n" )

        #fP.write("%f\t"   % windows_l[ind])
        #for j in range (len(E_vec)):
        #   fP.write("%.3e\t" %current_spectrum_pi0[j])
        #fP.write("\n" )


        #fT.write("%f\t"   % windows_l[ind])
        #for j in range (len(E_vec)):
        #   TOTAL0 = current_spectrum_brems[j]+current_spectrum_IC[j]+current_spectrum_pi0[j]
        #   fT.write("%.3e\t" % TOTAL0)
        #fT.write("\n" )

          
        
        # units of current_spectrum --> gammasky units [GeV^2/( cm^2 s sr GeV)]

# print "Plotting"

        #Just 1 model   
#       spec_tot, = plt.plot(E_vec, (current_spectrum_brems+current_spectrum_IC+current_spectrum_pi0), ls='-' , color=colors[iModel],  lw=2.2)
#       spec_brems, = plt.plot(E_vec, (current_spectrum_brems), ls=':' , color='red',  lw=2)
#       spec_IC,  = plt.plot(E_vec, (current_spectrum_IC), ls=':' , color='blue',  lw=2)
#       spec_pi0, = plt.plot(E_vec, (current_spectrum_pi0), ls=':' , color='green',  lw=2)
         

f.close()
#fE.close()
#fP.close()
#fT.close()
exit ()



    #current_spectrum_tot = current_spectrum_brems/(E_vec**2.) + current_spectrum_IC/(E_vec**2.) + current_spectrum_pi0/(E_vec**2.) + current_spectrum_fermi_exgb
    
#   print "IC spectrum"
#       print current_spectrum_IC
#       print "pi0 spectrum"
#    print current_spectrum_pi0
        
        # units of current_spectrum --> gammasky units [GeV^2/( cm^2 s sr GeV)]
        
#  print "Plotting"

#spec_brems, = plt.plot(E_vec, current_spectrum_brems*E_vec**0.8, ls='-' , color="red", lw=1.9)
#spec_IC, = plt.plot(E_vec, current_spectrum_IC*E_vec**0.8, ls='-' , color="blue",   lw=1.9)
#spec_pi0, = plt.plot(E_vec, current_spectrum_pi0*E_vec**0.8, ls='-' , color="green",  lw=1.9)
