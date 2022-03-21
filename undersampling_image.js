// Mask
var mask = ee.Image("path").select('CGIAR_PET').add(1)

// Dataset, transform to image
var global_fungi = ee.FeatureCollection('path')
var global_fungi_img = ee.Image.constant(1).updateMask(global_fungi.reduceToImage(['database_ID'], 'first'))

// Distance cost image
var dist = ee.Image.constant(1).cumulativeCost({
  source: global_fungi_img,
  maxDistance: 2500000,
  geodeticDistance: true
}).updateMask(mask)

// Color palette
var magma = ["000004", "1D1147", "51127C", "822681", "B63679", "E65164", "FB8861", "FEC287", "FCFDBF"]

// Exporting geometry
var exportingGeometry = ee.Geometry.Polygon([[[-180, 88], [180, 88], [180, -88], [-180, -88]]], null, false);

// Export to assets
Export.image.toAsset({image: dist,
                      description: 'distance_to_global_fungi_DB',
                      assetId:'pat',
                      crs: 'EPSG:4326',
                    	crsTransform: [0.08333333333333333,0,-180,0,-0.08333333333333333,90],
                    	region: exportingGeometry,
                    	maxPixels: 1e13,
})

// Load image after export is completed
dist = ee.Image('path')

// PCA interpolation extrapolation image
var int_ext = ee.Image('path2')

// Get min/max values
var minmaxVals_dist = dist.reduceRegion({reducer: ee.Reducer.minMax(),
                                       geometry: exportingGeometry,
                                       scale: 927.6624232772797,
                                       maxPixels: 1e13})
var minmaxVals_dist = minmaxVals_dist.rename(minmaxVals_dist.keys(), ['max','min'])

var minmaxVals_int_ext = int_ext.reduceRegion({reducer: ee.Reducer.minMax(),
                                       geometry: exportingGeometry,
                                       scale: 927.6624232772797,
                                       maxPixels: 1e13})
var minmaxVals_int_ext = minmaxVals_int_ext.rename(minmaxVals_int_ext.keys(), ['max','min'])

// Scale images
var dist_scaled = ee.Image.constant(1).subtract(dist.unitScale(minmaxVals_dist.get('min'), minmaxVals_dist.get('max')))
var int_ext_scaled = int_ext.unitScale(minmaxVals_int_ext.get('min'), minmaxVals_int_ext.get('max'))

// Combine images
var combined_img = ee.ImageCollection([dist_scaled.rename('b1'), dist_scaled.rename('b1'), int_ext_scaled.rename('b1')]).mean()

Map.addLayer(dist_scaled, {min:0, max:1, palette: magma}, 'dist_scaled', false)
Map.addLayer(int_ext_scaled, {min:0, max:1, palette: magma}, 'int_ext_scaled', false)
Map.addLayer(ee.Image.constant(1).subtract(combined_img), {min:0, max:.41, palette: magma}, 'combined_img', true)

Map.addLayer(global_fungi, {color:'red'}, 'Currently available data')

// Export to drive
Export.image.toDrive({image: combined_img,
                      description: 'distance_to_global_fungi_DB_exportToDrive',
                      crs: 'EPSG:4326',
                    	crsTransform: [0.08333333333333333,0,-180,0,-0.08333333333333333,90],
                    	region: exportingGeometry,
                    	maxPixels: 1e13,
})
