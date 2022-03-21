# Import the modules of interest
import pandas as pd
import numpy as np
import subprocess
import tqdm
import time
import datetime
import ee
import os
from functools import partial
from pathlib import Path
from contextlib import contextmanager
from scipy.spatial import ConvexHull
from itertools import combinations
from itertools import repeat
from sklearn.decomposition import PCA
from pathlib import Path

ee.Initialize()

# Configuration
####################################################################################
# Input the name of the username that serves as the home folder for asset storage
usernameFolderString = 'johanvandenhoogen'

# Input the normal wait time (in seconds) for "wait and break" cells
normalWaitTime = 5

# Input a longer wait time (in seconds) for "wait and break" cells
longWaitTime = 10

# Specify the column names where the latitude and longitude information is stored
latString = 'Pixel_Lat'
longString = 'Pixel_Long'

####################################################################################
# Image export settings
# Set pyramidingPolicy for exporting purposes
pyramidingPolicy = 'mean'

exportingGeometry = ee.Geometry.Polygon([[[-180, 88], [180, 88], [180, -88], [-180, -88]]], None, False);


####################################################################################
# Covariate input settings
# Input a list of the covariates being used
covariateList = [
'CGIAR_Aridity_Index',
'CGIAR_PET',
'CHELSA_BIO_Annual_Mean_Temperature',
'CHELSA_BIO_Annual_Precipitation',
'CHELSA_BIO_Isothermality',
'CHELSA_BIO_Max_Temperature_of_Warmest_Month',
'CHELSA_BIO_Mean_Diurnal_Range',
'CHELSA_BIO_Mean_Temperature_of_Coldest_Quarter',
'CHELSA_BIO_Mean_Temperature_of_Driest_Quarter',
'CHELSA_BIO_Mean_Temperature_of_Warmest_Quarter',
'CHELSA_BIO_Mean_Temperature_of_Wettest_Quarter',
'CHELSA_BIO_Min_Temperature_of_Coldest_Month',
'CHELSA_BIO_Precipitation_Seasonality',
'CHELSA_BIO_Precipitation_of_Coldest_Quarter',
'CHELSA_BIO_Precipitation_of_Driest_Month',
'CHELSA_BIO_Precipitation_of_Driest_Quarter',
'CHELSA_BIO_Precipitation_of_Warmest_Quarter',
'CHELSA_BIO_Precipitation_of_Wettest_Month',
'CHELSA_BIO_Precipitation_of_Wettest_Quarter',
'CHELSA_BIO_Temperature_Annual_Range',
'CHELSA_BIO_Temperature_Seasonality',
'CMS_SoilRespiration_mean',
'CSP_Global_Human_Modification',
# 'EarthEnvTexture_CoOfVar_EVI',
# 'EarthEnvTexture_Contrast_EVI',
# 'EarthEnvTexture_Correlation_EVI',
# 'EarthEnvTexture_Dissimilarity_EVI',
# 'EarthEnvTexture_Entropy_EVI',
# 'EarthEnvTexture_Evenness_EVI',
# 'EarthEnvTexture_Homogeneity_EVI',
# 'EarthEnvTexture_Maximum_EVI',
# 'EarthEnvTexture_Range_EVI',
'EarthEnvTexture_Shannon_Index',
'EarthEnvTexture_Simpson_Index',
# 'EarthEnvTexture_Std_EVI',
# 'EarthEnvTexture_Uniformity_EVI',
# 'EarthEnvTexture_Variance_EVI',
# 'EarthEnvTopoMed_1stOrderPartialDerivEW',
# 'EarthEnvTopoMed_1stOrderPartialDerivNS',
# 'EarthEnvTopoMed_2ndOrderPartialDerivEW',
# 'EarthEnvTopoMed_2ndOrderPartialDerivNS',
# 'EarthEnvTopoMed_AspectCosine',
# 'EarthEnvTopoMed_AspectSine',
# 'EarthEnvTopoMed_Eastness',
'EarthEnvTopoMed_Elevation',
# 'EarthEnvTopoMed_Northness',
# 'EarthEnvTopoMed_ProfileCurvature',
'EarthEnvTopoMed_Roughness',
# 'EarthEnvTopoMed_Slope',
# 'EarthEnvTopoMed_TangentialCurvature',
# 'EarthEnvTopoMed_TerrainRuggednessIndex',
'EarthEnvTopoMed_TopoPositionIndex',
'EarthEnvTopoMed_VectorRuggednessMeasure',
'EsaCci_AboveGroundBiomass',
'EsaCci_BurntAreasProbability',
'EsaCci_SnowProbability',
'FanEtAl_Depth_to_Water_Table_AnnualMean',
'GHS_Population_Density',
'GPWv4_Population_Density',
'GlobBiomass_AboveGroundBiomass',
'GlobBiomass_GrowingStockVolume',
'GoshEtAl_GlobalGDP2006',
'GranthamEtAl_ForestLandscapeIntegrityIndex',
'HansenEtAl_TreeCover_Year2000',
'HansenEtAl_TreeCover_Year2010',
'IPCC_Global_Biomass',
'MODIS_EVI',
# 'MODIS_FPAR',
# 'MODIS_GPP',
# 'MODIS_LAI',
'MODIS_NDVI',
# 'MODIS_NPP',
'SG_Absolute_depth_to_bedrock',
'SG_Bulk_density_015cm',
'SG_CEC_015cm',
'SG_Clay_Content_015cm',
'SG_Coarse_fragments_015cm',
'SG_Depth_to_bedrock',
'SG_H2O_Capacity_015cm',
'SG_Probability_of_occurrence_of_R_horizon',
'SG_SOC_Content_015cm',
'SG_SOC_Density_015cm',
'SG_SOC_Stock_005cm_to_015cm',
'SG_Sand_Content_015cm',
'SG_Saturated_H2O_Content_015cm',
'SG_Silt_Content_015cm',
'SG_Soil_H2O_Capacity_pF_20_015cm',
'SG_Soil_H2O_Capacity_pF_23_015cm',
'SG_Soil_H2O_Capacity_pF_25_015cm',
'SG_Soil_pH_H2O_015cm',
'SG_Soil_pH_KCl_015cm',
'SpawnEtAl_HarmonizedAGBiomass',
'SpawnEtAl_HarmonizedBGBiomass',
'Yasso15_MeanAnnualDecompCoeff_Humus',
'Yasso15_MeanAnnualDecompCoeff_Lignin',
'Yasso15_MeanAnnualDecompCoeff_SugarsCelluloseWaxes',
	 ]

# Desert mask
# desertMask = compositeToClassify.select('WWF_Biome').neq(13)

# Load the composite on which to perform the mapping, and subselect the bands of interest
# compositeToClassify = ee.Image("projects/crowtherlab/Composite/CrowtherLab_Composite_30ArcSec").select(covariateList).updateMask(desertMask)
compositeToClassify = ee.Image("projects/crowtherlab/Composite/CrowtherLab_Composite_30ArcSec").select(covariateList)

sampled_data = pd.read_csv('data/20211008_global_fungi_data_sampled.csv')[covariateList]

# ##################################################################################################################################################################
# # Univariate int-ext analysis
# ##################################################################################################################################################################
# # Create a feature collection with only the values from the image bands
# fcForMinMax = fcOI.select(covariateList)
#
# # Make a FC with the band names
# fcWithBandNames = ee.FeatureCollection(ee.List(covariateList).map(lambda bandName: ee.Feature(ee.Geometry.Point([0,0])).set('BandName',bandName)))
#
# def calcMinMax(f):
#   bandBeingComputed = f.get('BandName')
#   maxValueToSet = fcForMinMax.reduceColumns(ee.Reducer.minMax(),[bandBeingComputed])
#   return f.set('MinValue',maxValueToSet.get('min')).set('MaxValue',maxValueToSet.get('max'))
#
# # Map function
# fcWithMinMaxValues = ee.FeatureCollection(fcWithBandNames).map(calcMinMax)
#
# # Make two images from these values (a min and a max image)
# maxValuesWNulls = fcWithMinMaxValues.toList(1000).map(lambda f: ee.Feature(f).get('MaxValue'))
# maxDict = ee.Dictionary.fromLists(covariateList,maxValuesWNulls)
# minValuesWNulls = fcWithMinMaxValues.toList(1000).map(lambda f: ee.Feature(f).get('MinValue'))
# minDict = ee.Dictionary.fromLists(covariateList,minValuesWNulls)
# minImage = minDict.toImage()
# maxImage = maxDict.toImage()
#
# totalBandsBinary = compositeToClassify.gte(minImage.select(covariateList)).lt(maxImage.select(covariateList))
# univariate_int_ext_image = totalBandsBinary.reduce('sum').divide(compositeToClassify.bandNames().length()).rename('univariate_pct_int_ext')

##################################################################################################################################################################
# Multivariate (PCA) int-ext analysis
##################################################################################################################################################################

# Input the proportion of variance that you would like to cover
propOfVariance = 90

# PCA interpolation/extrapolation helper function
def assessExtrapolation(fcOfInterest, propOfVariance):
	# Compute the mean and standard deviation of each band, then standardize the point data
	meanVector = fcOfInterest.mean()
	stdVector = fcOfInterest.std()
	standardizedData = (fcOfInterest-meanVector)/stdVector

	# Then standardize the composite from which the points were sampled
	meanList = meanVector.tolist()
	stdList = stdVector.tolist()
	bandNames = list(meanVector.index)
	meanImage = ee.Image(meanList).rename(bandNames)
	stdImage = ee.Image(stdList).rename(bandNames)
	standardizedImage = compositeToClassify.subtract(meanImage).divide(stdImage)

	# Run a PCA on the point samples
	pcaOutput = PCA()
	pcaOutput.fit(standardizedData)

	# Save the cumulative variance represented by each PC
	cumulativeVariance = np.cumsum(np.round(pcaOutput.explained_variance_ratio_, decimals=4)*100)

	# Make a list of PC names for future organizational purposes
	pcNames = ['PC'+str(x) for x in range(1,fcOfInterest.shape[1]+1)]

	# Get the PC loadings as a data frame
	loadingsDF = pd.DataFrame(pcaOutput.components_,columns=[str(x)+'_Loads' for x in bandNames],index=pcNames)

	# Get the original data transformed into PC space
	transformedData = pd.DataFrame(pcaOutput.fit_transform(standardizedData,standardizedData),columns=pcNames)

	# Make principal components images, multiplying the standardized image by each of the eigenvectors
	# Collect each one of the images in a single image collection

	# First step: make an image collection wherein each image is a PC loadings image
	listOfLoadings = ee.List(loadingsDF.values.tolist())
	eePCNames = ee.List(pcNames)
	zippedList = eePCNames.zip(listOfLoadings)
	def makeLoadingsImage(zippedValue):
		return ee.Image.constant(ee.List(zippedValue).get(1)).rename(bandNames).set('PC',ee.List(zippedValue).get(0))
	loadingsImageCollection = ee.ImageCollection(zippedList.map(makeLoadingsImage))

	# Second step: multiply each of the loadings image by the standardized image and reduce it using a "sum"
	# to finalize the matrix multiplication
	def finalizePCImages(loadingsImage):
		PCName = ee.String(ee.Image(loadingsImage).get('PC'))
		return ee.Image(loadingsImage).multiply(standardizedImage).reduce('sum').rename([PCName]).set('PC',PCName)
	principalComponentsImages = loadingsImageCollection.map(finalizePCImages)

	# Choose how many principal components are of interest in this analysis based on amount of
	# variance explained
	numberOfComponents = sum(i < propOfVariance for i in cumulativeVariance)+1
	print('Number of Principal Components being used:',numberOfComponents)

	# Compute the combinations of the principal components being used to compute the 2-D convex hulls
	tupleCombinations = list(combinations(list(pcNames[0:numberOfComponents]),2))
	print('Number of Combinations being used:',len(tupleCombinations))

	# Generate convex hulls for an example of the principal components of interest
	cHullCoordsList = list()
	for c in tupleCombinations:
		firstPC = c[0]
		secondPC = c[1]
		outputCHull = ConvexHull(transformedData[[firstPC,secondPC]])
		listOfCoordinates = transformedData.loc[outputCHull.vertices][[firstPC,secondPC]].values.tolist()
		flattenedList = [val for sublist in listOfCoordinates for val in sublist]
		cHullCoordsList.append(flattenedList)

	# Reformat the image collection to an image with band names that can be selected programmatically
	pcImage = principalComponentsImages.toBands().rename(pcNames)

	# Generate an image collection with each PC selected with it's matching PC
	listOfPCs = ee.List(tupleCombinations)
	listOfCHullCoords = ee.List(cHullCoordsList)
	zippedListPCsAndCHulls = listOfPCs.zip(listOfCHullCoords)

	def makeToClassifyImages(zippedListPCsAndCHulls):
		imageToClassify = pcImage.select(ee.List(zippedListPCsAndCHulls).get(0)).set('CHullCoords',ee.List(zippedListPCsAndCHulls).get(1))
		classifiedImage = imageToClassify.rename('u','v').classify(ee.Classifier.spectralRegion([imageToClassify.get('CHullCoords')]))
		return classifiedImage

	classifedImages = ee.ImageCollection(zippedListPCsAndCHulls.map(makeToClassifyImages))
	finalImageToExport = classifedImages.sum().divide(ee.Image.constant(len(tupleCombinations)))

	return finalImageToExport

# PCA interpolation-extrapolation image
PCA_int_ext = assessExtrapolation(sampled_data.dropna(how='any'), propOfVariance).rename('PCA_pct_int_ext')

exportTask = ee.batch.Export.image.toAsset(
	image = PCA_int_ext.toFloat(),
	description = 'AMF_PCA_int_ext',
	assetId = 'users/johanvandenhoogen/2021_nat_micro/global_fungi_PCA_int_ext' ,
	crs = 'EPSG:4326',
	crsTransform = '[0.008333333333333333,0,-180,0,-0.008333333333333333,90]',
	region = exportingGeometry,
	maxPixels = int(1e13),
	pyramidingPolicy = {".default": pyramidingPolicy}
);
exportTask.start()
print('Image export task started, moving on')
