#----------------------------------------------------------#
# Research Techniques in Astronomy        		           #
# Identifying solar system objects						   #
#														   #
#----------------------------------------------------------#

import numpy as np
from astropy.io import fits
import matplotlib.pylab as pl
#import os
from scipy.ndimage import shift

# 2013091
print(' ')
print('This program is designed to identify Neptune and its moons from a')
print('series of images located in a folder called "Questdata"')
print(' ')

# Reading in the Science images and headers
sci_files = [234901, 1020246, 1040543]
sci_head = []
master_sci_data = []

print(' ')
print('The exposure times for the science observations are:')
for expId in sci_files:
	sci_data = np.array(fits.getdata("QUESTdata/science/2013091%07ds.C22.fits" % expId),
					  dtype = np.float64)
	# assign the header to a temporary variable
	header = fits.getheader("QUESTdata/science/2013091%07ds.C22.fits" % expId)
	# extract exposure time from the header
	print(header['EXPTIME'])
	# add the exposure time into the empty list created above
	sci_head.append(header)
	#same for images
	master_sci_data.append(sci_data)

master_sci_data = np.array(master_sci_data)


# Same process for the flats folder
flat_files = ["0231708e", "0231759e", "0231853e", "0231944e", "0232035e", "0232127e", "0232217e", "0232309e", "0232401e", "1095457m", "1095548m", "1095641m", "1095732m", "1095822m", "1095915m"]
flat_head = []
master_flat_data = []

print(' ')
print('The exposure times for the flat observations are:')
for expId in flat_files:
	flat_data = np.array(fits.getdata("QUESTdata/flats/2013091%s.C22.fits" % expId),
					  dtype = np.float64)
	header = fits.getheader("QUESTdata/flats/2013091%s.C22.fits" % expId)
	print(header['EXPTIME'])
	flat_head.append(header)
	master_flat_data.append(flat_data)

master_flat_data = np.array(master_flat_data)

# Same process for the darks folder
dark_files = [10, 180]
dark_head = []
master_dark_data = []

print(' ')
print('The exposure times for the dark observations are:')
for expId in dark_files:
	dark_data = np.array(fits.getdata("QUESTdata/darks/dark_%d.C22.fits" % expId),
					  dtype = np.float64)
	header = fits.getheader("QUESTdata/darks/dark_%d.C22.fits" % expId)
	print(header['EXPTIME'])
	dark_head.append(header)
	master_dark_data.append(dark_data)

master_dark_data = np.array(master_dark_data)

print(' ')
print('By comparing the exposure times we can deduce that the first dark field')
print('observation belongs to the flat observation, and the second belongs')
print('to the science observation.')
print(' ')


#-----------------------------------Bias----------------------------------------

# Create empty list to put the Bias overscan region into.
Master_Sci_Bias = []
# Measure the median bias from each row in the overscan region 
for i in range(3):
	# Create an average for the 597-640 section of each row.
	SciBiasPerRow = np.median(master_sci_data[i,:,597:640], axis = 1)
	Master_Sci_Bias.append(SciBiasPerRow) 
# Turn lists into array
Master_Sci_Bias = np.array(Master_Sci_Bias)
# Reshape the science images, clipping out the bias overscan region
master_sci_data = master_sci_data[:,:,0:597]
# Subtract the median bias of each row, from every row of the science images
for i in range(3):
	for j in range(597):
		master_sci_data[i,:,j] =master_sci_data[i,:,j] - Master_Sci_Bias[i,:]


# Same process for the flat images

Master_Flat_Bias = []
for i in range(15):
	FlatBiasPerRow = np.median(master_flat_data[i,:,597:640], axis =1)
	Master_Flat_Bias.append(FlatBiasPerRow)

Master_Flat_Bias = np.array(Master_Flat_Bias)
master_flat_data = master_flat_data[:,:,0:597]

for i in range(15):
	for j in range(597):
		master_flat_data[i,:,j] -= Master_Flat_Bias[i,:]



#-----------------------------Dark Subtraction----------------------------------

master_dark_data = master_dark_data[:,:,0:597]

# Subtract the already bias corrected dark data from the science and flat images

for i in range(3):
	master_sci_data[i,:,:] -= master_dark_data[1,:,:]

for i in range(15):
	master_flat_data[i,:,:] -= master_dark_data[0,:,:]



#-------------------------Creating a flat field --------------------------------

# Median combine all 15 images
mFlat = np.median(master_flat_data, axis = 0)
# Normalise the master flat field images, so that the media pixel value = 1.
nFlat = mFlat/ np.median(mFlat.flatten( ))

#-------------------Plotting Pixel Values as a Histogram------------------------

# replaced normed=True to density=True as it was depreciated
pl.hist(nFlat.flatten(), bins=50, range=(0.7,1.3),density=True,
		histtype='stepfilled',facecolor='m')
pl.title(r"$Normalised \ Pixel \ values \ in \ Flat \ Image$")
pl.xlabel(r"$Pixel \ Value$")
pl.ylabel(r"$Frequency \ Density$")
pl.show()


# clobber=True has been deprciated and changed to overwrite=True
fits.writeto("QuestMasterFlat.fits", nFlat, output_verify = 'ignore', overwrite=True)



#-------------------Flat field correction of science images---------------------
master_sci_data = master_sci_data/nFlat

#-------------------Subtracting the sky from reduced image----------------------

print('Subtracting the sky')
# Determine the minimum and maximum values of pixel intensity in the images
minpix = min(master_sci_data.flatten())
maxpix = max(master_sci_data.flatten())
# Use this as the range for histogram
rng = int(maxpix-minpix)

sciSky = []
for i in range(3):
	hf = np.histogram(master_sci_data[i].flatten(),bins=rng,range=(minpix,maxpix))
	skyval = hf[1][ hf[0].argmax() ]
	skySub = master_sci_data[i] - skyval
	sciSky.append(skySub)
print('Scisky Complete')
print('')

#----------------------------Shifting the images--------------------------------

# Overwrite the images to the shifted images 
sciSky = np.array(sciSky)

# Shifting frames 1 and 2 so that the background stars allign with the background
# stars of frame 0.
sciSky[1,:,:] = shift(sciSky[1,:,:], (17,-4))
sciSky[2,:,:] = shift(sciSky[2,:,:], (39,-7))

# Save the reduced images to a file
for i in range(3):
	hdu=fits.PrimaryHDU(master_sci_data[i])
	expid = i + 1
	fits.writeto("QUESTdata/unshiftedimage%d.fits" %expid ,master_sci_data[i],sci_head[i], output_verify = 'silentfix', overwrite= True)
	fits.writeto("QUESTdata/FinalImage%d.fits" %expid, sciSky[i], sci_head[i], output_verify = "silentfix", overwrite= True)
print('Reduced Images Saved to QUESTdata folder and ready for viewing')
print('')
