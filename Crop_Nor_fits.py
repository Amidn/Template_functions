import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
#import pyfits
from astropy.io import fits
import os
mpl.rc("font", family="serif", size=14)
from decimal import *

Model_source = 'pi0_DRAGON_BASE_v4_Amid.fits'
DIR = ''
Model  = os.path.join(DIR,Model_source)

min_En = 0.3e3
max_En = 100e3
region = 'R_TEST' # "A1" # 'R_TEST' # "E" ,  "F"

# ======================================
if region == "R_TEST":
    maxl=  120.  #
    minl=  0.0   #
    
    maxb=  30.   # 190.5
    minb= -30.   # 170.5

if region == "A1":
    maxl=  15.  #
    minl=  10.   #
    
    maxb=  5.   # 190.5
    minb= -5.   # 170.5

if region == "F":
    maxl=  73.  #
    minl=  64.   #
    
    maxb=  5.   # 190.5
    minb= -5.   # 170.5

if region == "E":
    maxl=  64.  # 488.5
    minl=  56.   # 472.5

    maxb=  5.   # 190.5
    minb= -5.   # 170.5

interv_l = round(maxl - minl)
interv_b = round(maxb - minb)



maxlpix = 360.5 - 2.0 * maxl
minlpix = 360.5 - 2.0 * minl

maxbpix =  180.5 + 2.0 * maxb
minbpix =  180.5 + 2.0 * minb


n_pix_l = abs(int(maxlpix - minlpix )) +1
n_pix_b = abs(int( maxbpix - minbpix + 1))
print "---------------------------------------"
print "n_pixel: ", n_pix_l, n_pix_b
#interv_l_npix =
#interv_b_npix =


hdulist = fits.open( Model )
print hdulist.info()
print hdulist[0].data.shape

print  Model, " is read"

hdu_num = 0
hdr_new = hdulist[hdu_num].header
#data = hdulist[hdu_num].data # I'll also take the data for comparison.
new_data = np.zeros((n_pix_b, n_pix_l))
print "new_data", np.shape(new_data)



Emin = 1.0
E_factor = 1.35
E_vec = np.zeros(40)
for i in xrange(0,40):
    E_vec[i] = Emin*pow(E_factor,i)
print "E_vec >>>>>>>>>>>>>>>>>", E_vec

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
for i in range (E_start, E_finish):
    print "Energy range >>>",  E_vec[i]
print "--------------------"




xx =  0
spec_tot=0
for x in range (int(maxlpix),int(minlpix +1 )):
    
    yy = 0
    # print  " X=", xx
    for y in range (int(minbpix), int(maxbpix +1  )):
        
        #  print "yy=", yy
        spect=0
        for E in range (E_start, E_finish):
            spec = (E_vec[E + 1] - E_vec[E]) * ( hdulist[0].data[ E , y, x ]/(E_vec[E]**2.)   +  hdulist[0].data[ E + 1 , y, x ]/(E_vec[E + 1]**2.)    ) / 2.0
            spect += spec
            spec_tot += spec
        #print "for ", x ," - ", y, "the energy is ",  Ebin
        new_data[ yy, xx ] = spect
        yy += 1
#print xx
    xx += 1

print "total spectrum>>>>>>>>>>>>>>>>>>>>>>>> ", spec_tot

xx =  0
for x in range (int(maxlpix),int(minlpix +1 )):
    yy = 0
    for y in range (int(minbpix), int(maxbpix +1)):
        new_data[ yy, xx ] =  new_data[ yy, xx ] / spec_tot
        yy += 1
    xx += 1


xx =0
Test=0
for x in range (int(maxlpix),int(minlpix +1 )):
    yy = 0
    for y in range (int(minbpix), int(maxbpix +1)):
        Test +=  new_data[ yy, xx ]
        yy += 1
    xx += 1
print "---------------------------------"

print "Total Spectrum = ", Test
getcontext().prec = 7
if Test == 1 :
    print "File is Normalized Successfully"
else:
    print "File is Not Normalized Successfully"

print "--------------------------------"


hdu_new = fits.PrimaryHDU(new_data)

newhdr = hdu_new.header




newhdr["NAXIS"] = (2,"data axes")
#newhdr["BITPIX"] = (-32,"array data type")
newhdr["NAXIS1"]=   (interv_l,"length of data in axis 1 (longitude)")
newhdr["NAXIS2"]=   (interv_b ,"length of data in axis 2 (latitude)")

newhdr["CTYPE1"] = 'GLON' # "RA" # 'pixel' #
newhdr["CTYPE2"] = 'GLAT'    # "DEC"# ' 'pixel'  #

newhdr['CRPIX1'] = 1
newhdr['CRPIX2'] = 1

newhdr["CRVAL1"]= (maxl,"Start of axis 1 (longitude, degree)")
newhdr["CDELT1"]= (-0.5,"Increment of axis 1 (longitude, degree)")

newhdr["CRVAL2"]= (minb,"Start of axis 2 (latitude,degree)")
newhdr["CDELT2"]= (0.5,"Increment of axis 2 (latitude,degree)")

hdulist_new = fits.HDUList([hdu_new])

hdulist_new.writeto(str(region)+'_2D_' + str(Model_source) ,clobber=True)
print "-----------------"
print hdulist_new.info()
print "-----------------"



hdulist_new.close()


print "---------------------------------"

print "Total Spectrum = ", Test
getcontext().prec = 12
#print Decimal(Test) + Decimal('0.0')
Desi_Test = Decimal(Test) + Decimal('0.0')
if Desi_Test == 1 :
    print "File is Normalized Successfully"
else:
    print "File is Not Normalized Successfully"

print "--------------------------------"



print "Finished Successfully"

