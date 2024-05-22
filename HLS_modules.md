```JavaScript
/* -------------Modules (or functions) necessary to run the HLS code------------------*/
//Function to select the first image when there is more than one image per date
//This procedure must be used for proper processing of S2 imagery
exports.FirstDate = function (imgs){
  //Simplify date to exclude time of day
  imgs = imgs.map(function(img){
    var d = ee.Date(img.get('system:time_start'));
    var day = d.get('day');
    var m = d.get('month');
    var y = d.get('year');
    var simpleDate = ee.Date.fromYMD(y,m,day);
    // return img.set('simpleTime',simpleDate.millis());
    return ee.Image(img.copyProperties(img, img.propertyNames())).set('simpleTime',simpleDate.millis());
  });
  //Find the unique days
  var days = ee.Dictionary(imgs.aggregate_histogram('simpleTime')).keys();
  imgs = days.map(function(d){
    d = ee.Number.parse(d);
    d = ee.Date(d);
    var t = imgs.filterDate(d,d.advance(1,'day'));
    var f = ee.Image(t.first());
    // t = t.mosaic();//The moisaic() changes the original projection
    // t = t.set('system:time_start',d.millis());
    // t = t.copyProperties(f);
    // t = t.copyProperties(f, f.propertyNames());
    // return ee.Image(t);
    return f;
    
    });
    imgs = ee.ImageCollection.fromImages(imgs);
      // This needs to happen AFTER the mosaicking step or else we still have edge artifacts
    // imgs = imgs.map(function(img){return img.updateMask(img.mask().reduce(ee.Reducer.min()))});
    return imgs;
}

// ######################################################################################################
// ######################################################################################################
//                                    ### STEP ATMOSPHERIC CORRECTION  ###
// ######################################################################################################
//Atmospheric correction via SIAC (https://github.com/MarcYin/SIAC_GEE)
var siac = require('users/marcyinfeng/utils:SIAC');

var Select_L7_bands = ['B1','B2','B3','B4','B5','B7'];
exports.SIAC_L7 = function(image){
  var boa = siac.get_l7_sur(image);
  return ee.Image(boa.copyProperties(image, image.propertyNames()));//.select(Select_L7_bands);
                      // .addBands(image.select('BQA')).select(Select_L7_bands);//'BQA' is used for cloud masking
};

var Select_L8_bands = ['B2', 'B3','B4', 'B5','B6','B7'];
exports.SIAC_L8 = function(image){
  var boa = siac.get_l8_sur(image);
  return ee.Image(boa.copyProperties(image, image.propertyNames()));//.select(Select_L8_bands);
                      // .addBands(image.select('BQA')).select(Select_L8_bands);//'BQA' is used for cloud masking;
};

var Select_S2_bands = ['B2','B3', 'B4','B8','B11', 'B12'];
exports.SIAC_S2 = function(image){
  var boa = siac.get_sur(image);
  return ee.Image(boa.copyProperties(image, image.propertyNames()));//.select(Select_S2_bands);
                          // .addBands(image.select('B9').toDouble());//'B9' is used for cloud masking
};




 
// ######################################################################################################
// ######################################################################################################
//                                ### STEP. CLOUD AND SHADOW MASKING  ###
// ######################################################################################################   

//### - LANDSAT 457 - (BOA original data) ###   
//<https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C01_T1_SR>
exports.cloudMaskL457_SRoriginal = function(image) {
  var qa = image.select('pixel_qa');
  var cloud = qa.bitwiseAnd(1 << 5)//Cloud
                  .and(qa.bitwiseAnd(1 << 7))//Cloud confidence
                  .or(qa.bitwiseAnd(1 << 3));//Cloud shadow
  // Remove edge pixels that don't occur in all bands
  var mask2 = image.mask().reduce(ee.Reducer.min());
  return image.updateMask(cloud.not()).updateMask(mask2);
};


//<https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T1_SR>
exports.cloudMaskL8_SRoriginal = function(image) {
  var cloudShadowBitMask = (1 << 3);//Cloud shadow
  var cloudsBitMask = (1 << 5);//Cloud
  var qa = image.select('pixel_qa');
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
};



//### - LANDSAT 7-8 (TOA original data) ###  

var extractQABits = function (qaBand, bitStart, bitEnd) {
      var numBits = bitEnd - bitStart + 1;
      //bitStart: bit Start Confidence (from the BQA product)
      //bitEnd: bit End Confidence
      var RADIX = 2;  // Radix for binary (base 2) data.
      var qaBits = qaBand.rightShift(bitStart).mod(Math.pow(RADIX, numBits));
  return qaBits.gte(2);//Pixels labeled as medium (2) or high (3) confidence are flagged (i.e., gte(2)).
};

exports.cloudMaskL78_TOAoriginal = function(image) {
  // var qaBand = image.select('BQA');//BQA is for TOA; 'pixel_qa' for BOA..Collection-1 - deprecated
  var qaBand = image.select('QA_PIXEL');//BQA is for TOA; 'pixel_qa' for BOA..Collection-2
   
  // -- Cloud confidence bits
  var Cloud_confidence = extractQABits(qaBand, 8, 9);
  //-- Shadow confidence bits
  var Shadow_confidence = extractQABits(qaBand, 10, 11); 
  //-- Snow/ice confidence bits
  var Snow_confidence = extractQABits(qaBand, 12, 13); 
  //-- Cirrus confidence bits
  var Cirrus_confidence = extractQABits(qaBand, 14, 15);  
  
  // -- Cloud + Shadow + Snow + Cirrus  combined mask
  var maskComposite = (Cloud_confidence.or(Shadow_confidence).or(Snow_confidence).or(Cirrus_confidence)).not();
  return image.updateMask(maskComposite);
};



// ###  SENTINEL 2  ###
// Code from <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_CLOUD_PROBABILITY>///

// funcHLS.edgeJoin(S2_col, S2clouds_col).map(funcHLS.cloudMaskS2)

function cloudProbabilityS2 (MAX_CLOUD_PROBABILITY) {
  var temp = function(img){
    var clouds = ee.Image(img.get('cloud_mask')).select('probability');
    var isNotCloud = clouds.lt(MAX_CLOUD_PROBABILITY);
    return img.updateMask(isNotCloud);
  };
  return temp;
}
// //For an explanation on the S2 edge issue see <https://hls.gsfc.nasa.gov/wp-content/uploads/2019/01/HLS.v1.4.UserGuide_draft_ver3.1.pdf>
function maskEdgesS2 (s2_img) {// Exclude bad data at scene edges.
    return s2_img.updateMask(
    s2_img.select('B8A').mask().updateMask(s2_img.select('B9').mask()));
    }


function CloudMask (S2_col, criteria, MAX_CLOUD_PROBABILITY){
 var S2clouds_col = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
                   .filter(criteria);
  print('S2clouds_col', S2clouds_col)                 
 var S2_col_maskedges = S2_col.map(maskEdgesS2);
// Join S2 SR with cloud probability dataset to add cloud mask.
print('Join S2 SR with cloud probability dataset to add cloud mask.')
  var s2SrWithCloudMask = ee.Join.saveFirst('cloud_mask').apply({
  primary: S2_col_maskedges,
  secondary: S2clouds_col,
  condition: ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'})
});
print('did work!')
  s2SrWithCloudMask = ee.ImageCollection(s2SrWithCloudMask);
  var tt = s2SrWithCloudMask.map(cloudProbabilityS2(MAX_CLOUD_PROBABILITY));
  return tt;
}

//To check whether the data is in the 0-1 or 16-bit scale. 
function getNumberScale (S2_col, bandName, ground_sensor){
          var t1 = S2_col.first().select(bandName).reduceRegion({reducer: ee.Reducer.mean(), 
          geometry: ground_sensor.buffer(30).bounds(),
          scale: 30});
          return ee.Number(t1.get(bandName));
          }


//---Find CLOUD SHADOWS with TDOM: Return Cloud + shadows masked data
//TDOM (Temporal Dark Outlier Mask) for finding dark outliers in time series. Masks pixels that are dark, and dark outliers

exports.Cloud_ClShadowS2 = function(args){
 
  //Cloud masked series of S2 data
  var S2_col_masked = CloudMask(args.Sentinel2Col, args.FilterCriteria, args.cloud_prob);
  
  var bandName_Initial = (S2_col_masked.first().bandNames()).getInfo();
  //Check whether the common band names are already defined
  //TDOM uses the NIR and SWIR bands for computation
  var bandName = (S2_col_masked.first().bandNames().get(0)).getInfo();
  var shadowSumBands = [];//NIR and SWIR1 bands
    if(bandName == 'B1'){
        shadowSumBands = ['B8','B11']; //B8=nir, B11=swir1
    }else{
       shadowSumBands = ['nir','swir1'];
    }

  var shadowSumBands_temp = ['S2_nir', 'S2_swir1'];

  
   // Get some pixel-wise stats for the time series
  var irMean;var irStdDev;
  if(args.PrecomputedTDOMstats === null || args.PrecomputedTDOMstats === undefined){
    print('Computing irMean & StdDev for TDOM');
    var Band_Stats = S2_col_masked.select(shadowSumBands, shadowSumBands_temp).map(function(img){return img.divide(10000)});
    irMean =  Band_Stats.mean();
    irStdDev = Band_Stats.reduce(ee.Reducer.stdDev());
  }else{
    print('Using pre-computed irMean & StdDev for TDOM');
    irMean = args.PrecomputedTDOMstats.select(['S2_nir_mean','S2_swir1_mean']);
    irStdDev = args.PrecomputedTDOMstats.select(['S2_nir_stdDev','S2_swir1_stdDev']);
  }
  
  
  //Mask out dark outliers,i.e., cloud shadows
    S2_col_masked = S2_col_masked.map(function(img){
    var preZ = img.select(shadowSumBands).divide(10000);
    var z = preZ.subtract(irMean).divide(irStdDev);
    var irSum = preZ.reduce(ee.Reducer.sum());
    var m = z.lt(args.zShadowThresh).reduce(ee.Reducer.sum()).eq(2).and(irSum.lt(args.irSumThresh))
             .focal_min(args.contractPixels).focal_max(args.dilatePixels)
             .not();
    return img.updateMask(img.mask().and(m));
     });

  return S2_col_masked.select(bandName_Initial);
};

// -------Remove clouds and cloud shadows--------
//https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_CLOUD_SCORE_PLUS_V1_S2_HARMONIZED#bands

exports.Cloud_ScorePlusS2 = function(args){
  
  // Cloud Score+ image collection. Note Cloud Score+ is produced from Sentinel-2
// Level 1C data and can be applied to either L1C or L2A collections.
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

// Use 'cs' or 'cs_cdf', depending on your use case; see docs for guidance.
var QA_BAND = 'cs_cdf';

// The threshold for masking; values between 0.50 and 0.65 generally work well.
// Higher values will remove thin clouds, haze & cirrus shadows.
var CLEAR_THRESHOLD = 0.60;

// Apply the cloud mask.
var s2_masked = args.Sentinel2Col
                .linkCollection(csPlus, [QA_BAND])
                .map(function(img) {
                return img.updateMask(img.select(QA_BAND).gte(CLEAR_THRESHOLD));
                });
    
return s2_masked;
  
  };


// ######################################################################################################
// ######################################################################################################
//                                ### STEP 04. BRDF and TOPO CORRECTION  ###
// ######################################################################################################
//Original Source: https://doi.org/10.3390/rs11070831
//Adapted by <https://github.com/ndminhhus/geeguide/blob/master/04.topo_correction.md>


////// ---- BRDF CORRECTION ----//////


    var PI = ee.Number(3.14159265359); var MAX_SATELLITE_ZENITH = 7.5;
    var MAX_DISTANCE = 1000000;        var UPPER_LEFT = 0;
    var LOWER_LEFT = 1;                var LOWER_RIGHT = 2;
    var UPPER_RIGHT = 3;


exports.applyBRDF = function(image){
    var date = image.date();
    var footprint = ee.List(image.geometry().bounds().bounds().coordinates().get(0));
    var angles =  getsunAngles(date, footprint);
    var sunAz = angles[0];
    var sunZen = angles[1];
  
    var viewAz = azimuth(footprint);
    var viewZen = zenith(footprint);
  
  
    var kval = _kvol(sunAz, sunZen, viewAz, viewZen);
    var kvol = kval[0];
    var kvol0 = kval[1];
    var result = _apply(image, kvol.multiply(PI), kvol0.multiply(PI));
  
    return result;
}
function getsunAngles(date, footprint){
  var jdp = date.getFraction('year');
  var seconds_in_hour = 3600;
  var  hourGMT = ee.Number(date.getRelative('second', 'day')).divide(seconds_in_hour);
    
  var latRad = ee.Image.pixelLonLat().select('latitude').multiply(PI.divide(180));
  var longDeg = ee.Image.pixelLonLat().select('longitude');
    
  // Julian day proportion in radians
  var jdpr = jdp.multiply(PI).multiply(2);
    
  var a = ee.List([0.000075, 0.001868, 0.032077, 0.014615, 0.040849]);
  var meanSolarTime = longDeg.divide(15.0).add(ee.Number(hourGMT));
  var localSolarDiff1 = value(a, 0)
          .add(value(a, 1).multiply(jdpr.cos())) 
          .subtract(value(a, 2).multiply(jdpr.sin())) 
          .subtract(value(a, 3).multiply(jdpr.multiply(2).cos())) 
          .subtract(value(a, 4).multiply(jdpr.multiply(2).sin()));

  var localSolarDiff2 = localSolarDiff1.multiply(12 * 60);
  
  var localSolarDiff = localSolarDiff2.divide(PI);
  var trueSolarTime = meanSolarTime 
          .add(localSolarDiff.divide(60)) 
          .subtract(12.0);
    
  // Hour as an angle;
  var ah = trueSolarTime.multiply(ee.Number(MAX_SATELLITE_ZENITH * 2).multiply(PI.divide(180))) ;   
  var b = ee.List([0.006918, 0.399912, 0.070257, 0.006758, 0.000907, 0.002697, 0.001480]);
  var delta = value(b, 0) 
        .subtract(value(b, 1).multiply(jdpr.cos())) 
        .add(value(b, 2).multiply(jdpr.sin())) 
        .subtract(value(b, 3).multiply(jdpr.multiply(2).cos())) 
        .add(value(b, 4).multiply(jdpr.multiply(2).sin())) 
        .subtract(value(b, 5).multiply(jdpr.multiply(3).cos())) 
        .add(value(b, 6).multiply(jdpr.multiply(3).sin()));

  var cosSunZen = latRad.sin().multiply(delta.sin()) 
        .add(latRad.cos().multiply(ah.cos()).multiply(delta.cos()));
  var sunZen = cosSunZen.acos();

  // sun azimuth from south, turning west
  var sinSunAzSW = ah.sin().multiply(delta.cos()).divide(sunZen.sin());
  sinSunAzSW = sinSunAzSW.clamp(-1.0, 1.0);
  
  var cosSunAzSW = (latRad.cos().multiply(-1).multiply(delta.sin())
                    .add(latRad.sin().multiply(delta.cos()).multiply(ah.cos()))) 
                    .divide(sunZen.sin());
  var sunAzSW = sinSunAzSW.asin();
  
  sunAzSW = where(cosSunAzSW.lte(0), sunAzSW.multiply(-1).add(PI), sunAzSW);
  sunAzSW = where(cosSunAzSW.gt(0).and(sinSunAzSW.lte(0)), sunAzSW.add(PI.multiply(2)), sunAzSW);
  
  var sunAz = sunAzSW.add(PI);
  // # Keep within [0, 2pi] range
    sunAz = where(sunAz.gt(PI.multiply(2)), sunAz.subtract(PI.multiply(2)), sunAz);
  
  var footprint_polygon = ee.Geometry.Polygon(footprint);
  sunAz = sunAz.clip(footprint_polygon);
  sunAz = sunAz.rename(['sunAz']);
  sunZen = sunZen.clip(footprint_polygon).rename(['sunZen']);
  
  return [sunAz, sunZen];
}
function azimuth(footprint){
    function x(point){return ee.Number(ee.List(point).get(0))}
    function  y(point){return ee.Number(ee.List(point).get(1))}
    
    var upperCenter = line_from_coords(footprint, UPPER_LEFT, UPPER_RIGHT).centroid().coordinates();
    var lowerCenter = line_from_coords(footprint, LOWER_LEFT, LOWER_RIGHT).centroid().coordinates();
    var slope = ((y(lowerCenter)).subtract(y(upperCenter))).divide((x(lowerCenter)).subtract(x(upperCenter)));
    var slopePerp = ee.Number(-1).divide(slope);
    var azimuthLeft = ee.Image(PI.divide(2).subtract((slopePerp).atan()));
    return azimuthLeft.rename(['viewAz']);
  }
function zenith(footprint){
    var leftLine = line_from_coords(footprint, UPPER_LEFT, LOWER_LEFT);
    var rightLine = line_from_coords(footprint, UPPER_RIGHT, LOWER_RIGHT);
    var leftDistance = ee.FeatureCollection(leftLine).distance(MAX_DISTANCE);
    var rightDistance = ee.FeatureCollection(rightLine).distance(MAX_DISTANCE);
    var viewZenith = rightDistance.multiply(ee.Number(MAX_SATELLITE_ZENITH * 2)) 
          .divide(rightDistance.add(leftDistance)) 
          .subtract(ee.Number(MAX_SATELLITE_ZENITH)) 
          .clip(ee.Geometry.Polygon(footprint)) 
          .rename(['viewZen']);
    return viewZenith.multiply(PI.divide(180));
}
function _apply(image, kvol, kvol0){
      var f_iso = 0;
      var f_geo = 0;
      var f_vol = 0;
      //Factors from Claverie et al (2018)<https://doi.org/10.1016/j.rse.2018.09.002>
			var blue = _correct_band(image, 'blue', kvol, kvol0, f_iso=0.0774, f_geo=0.0079, f_vol=0.0372);
			var green = _correct_band(image, 'green', kvol, kvol0, f_iso=0.1306, f_geo=0.0178, f_vol=0.0580);
			var red = _correct_band(image, 'red', kvol, kvol0, f_iso=0.1690, f_geo=0.0227, f_vol=0.0574);
			//var re1 = _correct_band(image, 're1', kvol, kvol0, f_iso=0.2085, f_geo=0.0256, f_vol=0.0845);//S2 only bands <https://doi.org/10.3390/rs9121325>
			//var re2 = _correct_band(image, 're2', kvol, kvol0, f_iso=0.2316, f_geo=0.0273, f_vol=0.1003);
			//var re3 = _correct_band(image, 're3', kvol, kvol0, f_iso=0.2599, f_geo=0.0294, f_vol=0.1197);
      var nir = _correct_band(image, 'nir', kvol, kvol0, f_iso=0.3093, f_geo=0.0330, f_vol=0.1535);
      //var re4 = _correct_band(image, 're4', kvol, kvol0, f_iso=0.2907, f_geo=0.0410, f_vol=0.1611);
      var swir1 = _correct_band(image, 'swir1', kvol, kvol0, f_iso=0.3430, f_geo=0.0453, f_vol=0.1154);   
      var swir2 = _correct_band(image, 'swir2', kvol, kvol0, f_iso=0.2658, f_geo=0.0387, f_vol=0.0639);
			return image.select([]).addBands([blue, green, red, nir, swir1, swir2]);
			//return image.select([]).addBands([blue, green, red, nir,re1,re2,re3,nir,re4,swir1, swir2]);//For all the S2 bands
}

function _correct_band(image, band_name, kvol, kvol0, f_iso, f_geo, f_vol){
	//"""fiso + fvol * kvol + fgeo * kgeo"""
	var iso = ee.Image(f_iso);
	var geo = ee.Image(f_geo);
	var vol = ee.Image(f_vol);
	var pred = vol.multiply(kvol).add(geo.multiply(kvol)).add(iso).rename(['pred']);
	var pred0 = vol.multiply(kvol0).add(geo.multiply(kvol0)).add(iso).rename(['pred0']);
	var cfac = pred0.divide(pred).rename(['cfac']);
	var corr = image.select(band_name).multiply(cfac).rename([band_name]);
	return corr;
}
function _kvol(sunAz, sunZen, viewAz, viewZen){
	//"""Calculate kvol kernel.
	//From Lucht et al. 2000
	//Phase angle = cos(solar zenith) cos(view zenith) + sin(solar zenith) sin(view zenith) cos(relative azimuth)"""
			
	var relative_azimuth = sunAz.subtract(viewAz).rename(['relAz']);
	var pa1 = viewZen.cos().multiply(sunZen.cos());
	var pa2 = viewZen.sin().multiply(sunZen.sin()).multiply(relative_azimuth.cos());
	var phase_angle1 = pa1.add(pa2);
	var phase_angle = phase_angle1.acos();
	var p1 = ee.Image(PI.divide(2)).subtract(phase_angle);
	var p2 = p1.multiply(phase_angle1);
	var p3 = p2.add(phase_angle.sin());
	var p4 = sunZen.cos().add(viewZen.cos());
	var p5 = ee.Image(PI.divide(4));

	var kvol = p3.divide(p4).subtract(p5).rename(['kvol']);

	var viewZen0 = ee.Image(0);
	var pa10 = viewZen0.cos().multiply(sunZen.cos());
	var pa20 = viewZen0.sin().multiply(sunZen.sin()).multiply(relative_azimuth.cos());
	var phase_angle10 = pa10.add(pa20);
	var phase_angle0 = phase_angle10.acos();
	var p10 = ee.Image(PI.divide(2)).subtract(phase_angle0);
	var p20 = p10.multiply(phase_angle10);
	var p30 = p20.add(phase_angle0.sin());
	var p40 = sunZen.cos().add(viewZen0.cos());
	var p50 = ee.Image(PI.divide(4));

	var kvol0 = p30.divide(p40).subtract(p50).rename(['kvol0']);

	return [kvol, kvol0]}
function line_from_coords(coordinates, fromIndex, toIndex){
    return ee.Geometry.LineString(ee.List([
      coordinates.get(fromIndex),
      coordinates.get(toIndex)]));
}
function where(condition, trueValue, falseValue){
  var trueMasked = trueValue.mask(condition);
  var falseMasked = falseValue.mask(invertMask(condition));
      return trueMasked.unmask(falseMasked);
}
function invertMask(mask){
    return mask.multiply(-1).add(1);
}
function value(list,index){
    return ee.Number(list.get(index));
}




////// ---- TOPOGRAPHIC CORRECTION ----//////
    // var scale = 10;
    var dem_original = ee.Image("USGS/SRTMGL1_003");
    // var dem = dem_original;
    // var degree2radian = 0.01745;
exports.illuminationCondition = function(img){
  //Noticed that by setting the dem projection to the same as for the images,
  //it avoids problems in displaying (Map.addLayer) the processed images (overlapping pixels).
  var selected_crs = img.select('red').projection();
  var dem =  dem_original.reproject({ crs: selected_crs  });
    
  // Extract image metadata about solar position
  var SZ_rad = ee.Image.constant(ee.Number(img.get('SOLAR_ZENITH_ANGLE'))).multiply(3.14159265359).divide(180).clip(img.geometry().buffer(10000)); 
  var SA_rad = ee.Image.constant(ee.Number(img.get('SOLAR_AZIMUTH_ANGLE')).multiply(3.14159265359).divide(180)).clip(img.geometry().buffer(10000)); 
  // Creat terrain layers
  var slp = ee.Terrain.slope(dem).clip(img.geometry().buffer(10000));
  var slp_rad = ee.Terrain.slope(dem).multiply(3.14159265359).divide(180).clip(img.geometry().buffer(10000));
  var asp_rad = ee.Terrain.aspect(dem).multiply(3.14159265359).divide(180).clip(img.geometry().buffer(10000));
  
  // Calculate the Illumination Condition (IC)
  // slope part of the illumination condition
  var cosZ = SZ_rad.cos();
  var cosS = slp_rad.cos();
  var slope_illumination = cosS.expression("cosZ * cosS", 
                                          {'cosZ': cosZ,
                                          'cosS': cosS.select('slope')});
  // aspect part of the illumination condition
  var sinZ = SZ_rad.sin(); 
  var sinS = slp_rad.sin();
  var cosAziDiff = (SA_rad.subtract(asp_rad)).cos();
  var aspect_illumination = sinZ.expression("sinZ * sinS * cosAziDiff", 
                                          {'sinZ': sinZ,
                                            'sinS': sinS,
                                            'cosAziDiff': cosAziDiff});
  // full illumination condition (IC)
  var ic = slope_illumination.add(aspect_illumination);

  // Add IC to original image
  var img_plus_ic = ee.Image(img.addBands(ic.rename('IC')).addBands(cosZ.rename('cosZ')).addBands(cosS.rename('cosS')).addBands(slp.rename('slope')));
  return img_plus_ic;
  };
  
exports.illuminationConditionS2 = function(img){
  //Noticed that by setting the dem projection to the same as for the images,
  //it avoids problems in displaying (Map.addLayer) the processed images (overlapping pixels).
  var selected_crs = img.select('red').projection();
  var dem =  dem_original.reproject({ crs: selected_crs  });
  
  // Extract image metadata about solar position
  var SZ_rad = ee.Image.constant(ee.Number(img.get('MEAN_SOLAR_ZENITH_ANGLE'))).multiply(3.14159265359).divide(180).clip(img.geometry().buffer(10000)); 
  var SA_rad = ee.Image.constant(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')).multiply(3.14159265359).divide(180)).clip(img.geometry().buffer(10000)); 
  // Creat terrain layers
  var slp = ee.Terrain.slope(dem).clip(img.geometry().buffer(10000));
  var slp_rad = ee.Terrain.slope(dem).multiply(3.14159265359).divide(180).clip(img.geometry().buffer(10000));
  var asp_rad = ee.Terrain.aspect(dem).multiply(3.14159265359).divide(180).clip(img.geometry().buffer(10000));
  
  // Calculate the Illumination Condition (IC)
  // slope part of the illumination condition
  var cosZ = SZ_rad.cos();
  var cosS = slp_rad.cos();
  var slope_illumination = cosS.expression("cosZ * cosS", 
                                          {'cosZ': cosZ,
                                          'cosS': cosS.select('slope')});
  // aspect part of the illumination condition
  var sinZ = SZ_rad.sin(); 
  var sinS = slp_rad.sin();
  var cosAziDiff = (SA_rad.subtract(asp_rad)).cos();
  var aspect_illumination = sinZ.expression("sinZ * sinS * cosAziDiff", 
                                          {'sinZ': sinZ,
                                            'sinS': sinS,
                                            'cosAziDiff': cosAziDiff});
  // full illumination condition (IC)
  var ic = slope_illumination.add(aspect_illumination);

  // Add IC to original image
  var img_plus_ic = ee.Image(img.addBands(ic.rename('IC')).addBands(cosZ.rename('cosZ')).addBands(cosS.rename('cosS')).addBands(slp.rename('slope')));
  return img_plus_ic
  };


exports.illuminationCorrection = function(img){
    var props = img.toDictionary();
    //var st = img.get('system:time_start');
    
    var img_plus_ic = img;
    var mask1 = img_plus_ic.select('nir').gt(-0.1);
    var mask2 = img_plus_ic.select('slope').gte(5)
                            .and(img_plus_ic.select('IC').gte(0))
                            .and(img_plus_ic.select('nir').gt(-0.1));
    var img_plus_ic_mask2 = ee.Image(img_plus_ic.updateMask(mask2));
    
    // Specify Bands to topographically correct  
    var bandList = ['blue','green','red','nir','swir1','swir2']; 
    var compositeBands = img.bandNames();
    var nonCorrectBands = img.select(compositeBands.removeAll(bandList));
    
    var geom = ee.Geometry(img.get('system:footprint')).bounds().buffer(10000);
    
    function apply_SCSccorr(band){
      var method = 'SCSc';
      var out = img_plus_ic_mask2.select('IC', band).reduceRegion({
      reducer: ee.Reducer.linearFit(), // Compute coefficients: a(slope), b(offset), c(b/a)
      geometry: ee.Geometry(img.geometry().buffer(-100)), // trim off the outer edges of the image for linear relationship 
      scale: 10,
      maxPixels: 1000000000,
      });  

  if (out === null || out === undefined ){
      return img_plus_ic_mask2.select(band);
      }
  
  else{
      var out_a = ee.Number(out.get('scale'));
      var out_b = ee.Number(out.get('offset'));
      var out_c = out_b.divide(out_a);
      // Apply the SCSc correction
      var SCSc_output = img_plus_ic_mask2.expression(
        "((image * (cosB * cosZ + cvalue)) / (ic + cvalue))", {
        'image': img_plus_ic_mask2.select(band),
        'ic': img_plus_ic_mask2.select('IC'),
        'cosB': img_plus_ic_mask2.select('cosS'),
        'cosZ': img_plus_ic_mask2.select('cosZ'),
        'cvalue': out_c
      });
      
      return SCSc_output;
    }
      
    }
    
    var img_SCSccorr = ee.Image(bandList.map(apply_SCSccorr)).addBands(img_plus_ic.select('IC'));
    var bandList_IC = ee.List([bandList, 'IC']).flatten();
    img_SCSccorr = img_SCSccorr.unmask(img_plus_ic.select(bandList_IC)).select(bandList);
    
    // return img_SCSccorr.addBands(nonCorrectBands)
    //   .setMulti(props)
    //   .set('system:time_start',st);
    return ee.Image(img_SCSccorr.copyProperties(img, img.propertyNames()));//Modified Elias  

  };
  
  
  // ######################################################################################################
// ######################################################################################################
//                                ### STEP. REPROJECT and RESAMPLE  ###
// ######################################################################################################

exports.reproj_S2 = function(output_Scale){
  var resamp = function(img){
   var new_imag = img.reproject({
    crs: img.select('red').projection(),
    scale: output_Scale
  });
  return new_imag.set({'SATELLITE': 'SENTINEL_2' });//Add 'SATELLITE' for plotting purposes
  };   
  return resamp;
};


exports.aggregate_reproj_S2 = function(output_Scale){
  var resamp = function(img){
   var new_imag = img
    .reduceResolution({// Force the next reprojection to aggregate instead of resampling.
      reducer: ee.Reducer.mean(),
      maxPixels: 65535
    })   
   .reproject({
    crs: img.select('red').projection(),
    scale: output_Scale
  });
  return new_imag.set({'SATELLITE': 'SENTINEL_2' });//Add 'SATELLITE' for plotting purposes
  };   
  return resamp;
};


exports.reproj_L78 = function(selected_S2_crs_L){
  var apply_resample = function(img){
  var new_imag = img.reproject({crs: selected_S2_crs_L});
  return new_imag.set({'SATELLITE': img.get('SPACECRAFT_ID') });//Add 'SATELLITE' for plotting purposes
  };
  return apply_resample; 
};

exports.aggregate_reproj_L78 = function(selected_S2_crs_L){
  var apply_resample = function(img){
  var new_imag = img
    .reduceResolution({// Force the next reprojection to aggregate instead of resampling.
      reducer: ee.Reducer.mean(),
      maxPixels: 65535
    })  
  .reproject({crs: selected_S2_crs_L});
  return new_imag.set({'SATELLITE': img.get('SPACECRAFT_ID') });//Add 'SATELLITE' for plotting purposes
  };
  return apply_resample; 
};






// ######################################################################################################
// ######################################################################################################
//                                ###  BAND ADJUSTMENT   ###
// ######################################################################################################


// Cross-sensor transformation coefficients for TOA reflectance (for reflectance in the scale 0-1)

//  Band names_L8 =  ['B1-Aeros', 'B2-Blue', 'B3-Green', 'B4-Red', 'B5-NIR', 'B6-SWIR1', 'B7-SWIR2', 'B8-Pan', 'B9-Cir', 'B10_T1', 'B11_T2,'BQA' > Chastain et al (2019)
                      //'QA_RADSAT', 'SAA', 'SZA', 'VAA', 'VZA']
//Value = 1 is to maintain all bands in the output image
var interceptsL8 =   [  0       , -0.0107,    0.0026,    -0.0015,   0.0033,   0.0065,     0.0046  ,    0     ,    0    ,   0    ,     0   ,  0,
                          0,           0,      0,    0,    0];
var slopesL8 =       [  1       , 1.0946,     1.0043,    1.0524,    0.8954,   1.0049,     1.0002  ,    1     ,    1    ,   1    ,     1   ,  1, 
                           1,          1,      1,    1,    1];

//Band names_L7 =  ['B1-Blue',  'B2-Green', 'B3-Red', 'B4-NIR', 'B5-SWIR1', 'B6_VCID_1', B6_VCID_2, 'B7-SWIR2', 'B8-Pan', BQA,
                    //'QA_RADSAT', 'SAA', 'SZA', 'VAA', 'VZA']
var interceptsL7 = [-0.0139,      0.0041,    -0.0024,  -0.0076,    0.0041,       0,         0,        0.0086,       0,      0,
                     0,            0,      0,    0,    0];
var slopesL7     = [1.1060,       0.9909,    1.0568,   1.0045,     1.0361,       1,         1,        1.0401,       1,      1,
                       1,          1,      1,    1,    1];


exports.band_adjsut_L8 =  function(image){
  var b_adj = image.multiply(slopesL8).add(interceptsL8);
  return ee.Image(b_adj.copyProperties(image, image.propertyNames()));
};

exports.band_adjsut_L7 =  function(image){
  var b_adj = image.multiply(slopesL7).add(interceptsL7);
  return ee.Image(b_adj.copyProperties(image, image.propertyNames()));
};


// ######################################################################################################
// ######################################################################################################
//                                ### STEP. SPATIAL CO-REGISTRATION  ###
// ######################################################################################################

//Select the corrected cloud-free image pairs
exports.selDate = function(img_col, Selected_Date){
var Sel_1 = img_col.filterDate(Selected_Date, Selected_Date.advance(1, 'day'))
                      .toList(100);
return ee.Image(Sel_1.get(0));   
};


exports.Co_reg = function(L_Band, S2_L_Band, iterpolation_method){//Nested function
  var Disp = function(img){
    // Determine the displacement by matching a user-specified single band.
    var displacement_L_S2 = L_Band.displacement({
      referenceImage: S2_L_Band,
      //When tried values maxOffset>= 30 and 'nearest_neighbor', I noticed problems with the coregistered Landsat 8
      //Specifically, blank strips with no data
      maxOffset: 15.0,
      patchWidth: 100.0 
    });
    // Use the computed displacement to register all bands
    var L8_reg = img.displace({
      displacement: displacement_L_S2,
      mode: iterpolation_method//The interpolation mode to use: 'nearest_neighbor', 'bilinear' or 'bicubic'.)
    });
    return L8_reg;
  };
  return Disp;
};


exports.Repro = function(sel_CRS, sel_scale){//Nested function
  var Reprojected_img = function(img){
      var Repr_img = img.reproject({crs: sel_CRS, scale: sel_scale});
    return Repr_img;
  };
  return Reprojected_img;
};

// // ######################################################################################################
// // ######################################################################################################
// //                                    ### STEP: EXPORT COLLECTIONS ###
// // ######################################################################################################

// Define function to calculate NDVI.
exports.calc_NDVI = function (img) {
  var orig = img;
  img = img.normalizedDifference(['nir', 'red']).rename('NDVI');
  return ee.Image(img.copyProperties(orig, orig.propertyNames()));
};


//--Export paired images collection to GEE Asset---

/* 
 * Author: Rodrigo E. Principe
 * License: Apache 2.0
 * email: fitoprincipe82@gmail.com
https://github.com/fitoprincipe/geetools-code-editor/blob/60f41ab8932876d3ce55b9b221bd1669c142e85a/batch
*/
var get_options = function(def, options) {
  // fill options dict with default values if not present
  if (options !== undefined) {
    var opt = options;
    for (var key in def) {
      var value = def[key];
      if (opt[key] === undefined) {opt[key] = value}
    }
  } else {var opt = def}
  return opt;
};


var ImCol_toAsset = function(collection, assetFolder, options, ImPairName) {
  var root = ee.data.getAssetRoots()[0]['id'];
  var folder = assetFolder;
  if (folder !== null && folder !== undefined) {
    var assetFolder = root+'/'+folder+'/'
  } else {
    var assetFolder = root+'/'
  }
  
  var defaults = {
      name: null,
      scale: 1000,
      maxPixels: 1e13,
      region: null
    };
    
  var opt = get_options(defaults, options);

  
  var n = collection.size().getInfo();
  
    
  var colList = collection.toList(n);
  
  for (var i = 0; i < n; i++) {
    var img = ee.Image(colList.get(i));
    print('img.id().getInfo()', img.id().getInfo());
    print('image_+i.toString()', 'image_'+i.toString());
    
    var id = img.id().getInfo() || 'image_'+i.toString();
    print('id line 757', id);
    var region = opt.region || img.geometry().bounds().getInfo().coordinates;
    // var region = opt.region || img.geometry().bounds().getInfo()["coordinates"];

    // var assetId = assetFolder+id;
    var assetId = assetFolder+ImPairName.toString()+'_'+id;
    
    Export.image.toAsset({
      image: img,
      description: id,
      assetId: assetId,
      region: region,
      scale: opt.scale,
      maxPixels: opt.maxPixels});
  }
};



var funcStat = require('users/geeberra/Modules:HLS_Stats');

exports.Export_toAsset_ImgPairs = function  (args){
//----- Join pairs -----
var joinDiff = 24*args.CArgs.nDayDiff;//24 = 24 hours
var oliBns = args.bandOut.map(function(bn){return ee.String(bn).cat('_').cat('L8')});
var ETMBns = args.bandOut.map(function(bn){return ee.String(bn).cat('_').cat('L7')});

var L8_col = args.Collection_L8.select(args.bandIn_L8, args.bandOut).select(args.bandOut, oliBns);
var L7_col = args.Collection_L7.select(args.bandIn_L7, args.bandOut).select(args.bandOut, ETMBns);
var S2_col = args.Collection_S2.select(args.bandIn_S2, args.bandOut);

  // print('L8_col', L8_col)
  //   print('L7_col', L7_col)
 
  if(args.CArgs.ImPairs === null || args.CArgs.ImPairs === undefined){
    //This will join only when the three sats have corresponding obeservations within the nDayDiff
  var joined = funcStat.spatioTemporalJoin(L8_col, L7_col, joinDiff, 'L7');
  joined = ee.ImageCollection(funcStat.spatioTemporalJoin(joined, S2_col, joinDiff, 'S2')).map(MskPix);
  print('joining the L7, L8, S2 together');


  }else if(args.CArgs.ImPairs === 'L8_S2'){ //Only L8 x S2 observations 
   var joined_L8_S2 = ee.ImageCollection(funcStat.spatioTemporalJoin(L8_col, S2_col, joinDiff, 'S2')).map(funcStat.MskPix); 
      print('joined_L8_S2'); 
    var NT = args.ProcessingStep.toString();  
    ImCol_toAsset(joined_L8_S2, args.FolderAsset,
                {scale: args.CArgs.C_Scale, region: args.CArgs.studyArea, type: 'float'},
                NT+'_L8_S2/joined_L8_S2');  
    
    
  }else if(args.CArgs.ImPairs === 'L8_L7'){//Only L8 x L7 observations 
        var joined_L8_L7 = ee.ImageCollection(funcStat.spatioTemporalJoin(L8_col, L7_col, joinDiff, 'L7')).map(funcStat.MskPix);
   
      print('joined_L8_L7', joined_L8_L7);
    var NT = args.ProcessingStep.toString();  
    ImCol_toAsset(joined_L8_L7, args.FolderAsset,
                {scale: args.CArgs.C_Scale, region: args.CArgs.studyArea, type: 'float'},
                NT+'_L8_L7/joined_L8_L7');  
    
  }else if(args.CArgs.ImPairs === 'L7_S2'){//Only L7 x S2 observations 
        var joined_L7_S2 = ee.ImageCollection(funcStat.spatioTemporalJoin(L7_col, S2_col, joinDiff, 'S2')).map(funcStat.MskPix);
      print('joined_L7_S2', joined_L7_S2);
    
    var NT = args.ProcessingStep.toString();  
    ImCol_toAsset(joined_L7_S2, args.FolderAsset,
                {scale: args.CArgs.C_Scale, region: args.CArgs.studyArea, type: 'float'},
                NT+'_L7_S2/joined_L7_S2');  
    
  }else{
  var joined_L8_L7 = ee.ImageCollection(funcStat.spatioTemporalJoin(L8_col, L7_col, joinDiff, 'L7')).map(funcStat.MskPix);
  var joined_L8_S2 = ee.ImageCollection(funcStat.spatioTemporalJoin(L8_col, S2_col, joinDiff, 'S2')).map(funcStat.MskPix);
    L7_col = L7_col.select(args.bandOut, ETMBns);
  var joined_L7_S2 = ee.ImageCollection(funcStat.spatioTemporalJoin(L7_col, S2_col, joinDiff, 'S2')).map(funcStat.MskPix);
  print('joining S2xL8, S2xL7, L7xL8');
  // print('joined_L8_L7', joined_L8_L7)
  // print('joined_L8_S2', joined_L8_S2)
  // print('joined_L7_S2', joined_L7_S2)

//ImCol_toAsset = function(collection, assetFolder, options, ImPairName)
var NT = args.ProcessingStep.toString();
ImCol_toAsset(joined_L8_L7, args.FolderAsset,
                {scale: args.CArgs.C_Scale, region: args.CArgs.studyArea, type: 'float'},
                NT+'_L8_L7/joined_L8_L7');
ImCol_toAsset(joined_L8_S2, args.FolderAsset,
                {scale: args.CArgs.C_Scale, region: args.CArgs.studyArea, type: 'float'},
                NT+'_L8_S2/joined_L8_S2');
ImCol_toAsset(joined_L7_S2, args.FolderAsset,
                {scale: args.CArgs.C_Scale, region: args.CArgs.studyArea, type: 'float'},
                NT+'_L7_S2/joined_L7_S2');            
  }
};


// // ######################################################################################################
// // ######################################################################################################
// //                                    ### EXPORT FUNCTIONS ###
// // ######################################################################################################
exports.CloudMask = CloudMask;


```
