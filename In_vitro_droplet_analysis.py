### Script written by Simon Mergenthal
### Humboldt-Universität zu Berlin
### Mechanobiology (Prof. Klotzsch)
### 13.05.2023

### Version 2.0 for data with different amount of parameters
### Number of parameters has to be specified before running this script
### last updated: 05.02.2024

### Script to analyse droplets in solution in ImageJ

### input: single image .tif or .obf - files from Abberior STED
### images have to be all in a single folder, having all parameters in the image title
### parameters: protein construct, protein concentration, salt concentration and temperature
### some data only have part of parameters present (e.g. protein construct and temperature)

### import functions
from ij import IJ
from ij.io import DirectoryChooser
from java.io import File
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
import io
import time
from ij import process

### method to change signed 16-bit image to unsigned 16-bit image
def signed_16bit_to_unsigned (imp):
	stack = imp.getStack()
	if (stack.isVirtual()):
		IJ.error("Non-virtual stack required")
	cal = imp.getCalibration()
	if (not cal.isSigned16Bit()):
		IJ.error("Signed 16-bit image required")
	cal.disableDensityCalibration()
	ip = imp.getProcessor()
	min = ip.getMin()
	max = ip.getMax()
	stats = process.StackStatistics(imp)
	minv = stats.min
	for i in xrange(1,stack.getSize()+1):
		ip = stack.getProcessor(i)
		ip.add(-minv)
	imp.setStack(stack)
	ip = imp.getProcessor()
	ip.setMinAndMax(min-minv, max-minv)
	imp.updateImage()

### changing the amount of parameters given in the image title
# e.g. protein construct + protein concentration + salt concentration + temperature = 4
# e.g. protein construct + protein concentration = 2
# if you use other values than 2 or 4, please change section 'starting values' accordingly

parameters = 2

# to measure the execution time
st = time.time()

# for test purpose only (maybe also good to have later on)
IJ.run("Close All", "")
IJ.run("Clear Results", "")

dc = DirectoryChooser("File Directory")
path = dc.getDirectory()

print("")
print("Analysing folder: " + path)

filedic = File(path) # getting file object
imagelist = filedic.list() # getting all the file names in a list

# we only want to analyze the images
imagelist_raw = []

for image in imagelist:
	if ".obf" in image:
		imagelist_raw.append(image)
	elif ".tif" in image:
		imagelist_raw.append(image)
		
# creating a list to collect information about the number of images from different conditions
list_amount_images = []
counter_images_conditions = 0
counter_particles_conditions = 0
counter_particles_total = 0		

# starting value
if len(imagelist_raw) > 0:
	split_name = imagelist_raw[0].split("_")
	if parameters == 4:
		# protein construct + protein concentration + salt concentration + temperature
		old_name = split_name[0] + "_" + split_name[1] + "_" + split_name[3] + "_" + split_name[4]
	elif parameters == 2:
		# protein construct + parameter (protein concentration or temperature)
		old_name = split_name[0] + "_" + split_name[1]
	else:
		print("Please give another value of number of parameters or change section regarding starting value accordingly!")
		print("I will now only account for 2 parameters!")
else:
	print("There are no images to analyze in this folder!")

# making a reference for the RoiManager to work later with it
rm = RoiManager().getInstance()

# changing the displayed measurements in the results table
IJ.run("Set Measurements...", "area mean min shape feret's integrated median display redirect=None decimal=3")
rt = ResultsTable.getResultsTable()

results = []

for i in xrange(len(imagelist_raw)):
	# checking the name/condition of the image for analysis
	split_name = imagelist_raw[i].split("_")
	if parameters == 4:
		str_img = split_name[0] + "_" + split_name[1] + "_" + split_name[3] + "_" + split_name[4]
	elif parameters == 2:
		str_img = split_name[0] + "_" + split_name[1]
	# if the current image has the same condition than the last one, the counter goes up by 1
	if str_img == old_name:
		counter_images_conditions += 1
	# else the current values in old_name and the counter will be appended to the list with all condition
	# and the counter of the new condition will be reset to 1, whil the counter for cells with no particles is reset to 0
	else:
		list_amount_images.append(old_name)
		list_amount_images.append(str(counter_images_conditions))
		list_amount_images.append(str(counter_particles_conditions))
		old_name = str_img
		counter_images_conditions = 1
		counter_particles_conditions = 0

	# making sure the RoiManager is empty
	rm.reset()
	
	# open images from folder
	img = IJ.openImage(path+imagelist_raw[i])
	signed_16bit_to_unsigned(img)
	# saving the name of the image to use it as label in the results table
	label = imagelist_raw[i]
	
	# Background subtraction with Gaussian filter beforehand
	IJ.run(img, "Gaussian Blur...", "sigma=1") # Gaussian blur to get rid of many small particles with low intensity profile
	IJ.run(img, "Subtract Background...", "rolling=50 sliding disable") # subtract background with a rolling ball of 50 pixels
																 		# disable smoothing because already Gaussian blur before

	# duplicating the image to use auto threshold
	img_dupl = img.duplicate()
	
	# using Auto Threshold to get an binary image, which is needed to analyze particles
	IJ.run(img_dupl, "Auto Threshold", "method=Otsu white")
	IJ.run(img_dupl, "Invert", "")
		
	# using Analyze Particles Method
	# I get ROIs for every particle found in the image and one line in the results table
	# the old ROI from the cell will be deleted in the RoiManager
	IJ.run(img_dupl, "Analyze Particles...", "size=0.05-100 circularity=0.00-1.00 display exclude clear add")

	if rt.size() > 0: # otherwise an image, where no particle is found will cause the problem, that older ROIS will be measured
		
		counter_particles_conditions += rt.size()
		counter_particles_total += rt.size()

		# clearing the results to measure the ROIs in the original image with intensities
		IJ.run("Clear Results", "")
			
		# loop through all ROIs of one Image
		# measuring the intensity for each ROI in the original image
		# the results will be displayed in the "Results" table
		for index in xrange(rm.getCount()):
			rm.select(img, index)
			IJ.run(img, "Measure", "")

		# parsing the measured data on the intensity image in another array
		for row in xrange(rt.size()): # row is just the index for the row
			rt.setLabel(label, row)
			results.append(rt.getRowAsString(row))
	
	IJ.run("Close All", "")
	IJ.run("Clear Results", "")

# to get the last condition into the list
if len(imagelist_raw) > 0:
	list_amount_images.append(old_name)
	list_amount_images.append(str(counter_images_conditions))
	list_amount_images.append(str(counter_particles_conditions))

# saving all conditions in a file
file_conditions = open(path + "Conditions_of_used_protein.csv", "w")
file_conditions.write("Condition" + "\t" + "Analysed_Images" + "\t" + "Particles_Found" + "\n")
j = 0
while j < len(list_amount_images):
	file_conditions.write(list_amount_images[j].encode("utf-8") + "\t" + list_amount_images[j+1] + "\t" + list_amount_images[j+2] + "\n")
	j += 3
file_conditions.close()

# getting the first row for the output with all labels
headings_str = "Particle\tLabel\tArea\tMean\tMin\tMax\tCirc.\tFeret\tIntDen\tMedian\tRawIntDen\tFeretX\tFeretY\tFeretAngle\tMinFeret\tAR\tRound\tSolidity"

file_particles = io.open(path + "Results_particles.csv", "w", encoding="utf-8")
file_particles.write(headings_str + "\n")
for line in results:
	file_particles.write(line + "\n")
file_particles.close()

# calculating execution time
et = time.time()
elapsed_time = et - st

print("")
print("Task finished successfully")
print("Total number of analyzed images: " + str(len(imagelist_raw)))
print("Total number of particles detected: " + str(counter_particles_total))
print("Execution Time: " + str(time.strftime("%H:%M:%S", time.gmtime(elapsed_time))))
print("")