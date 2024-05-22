/* -------------- Harmonization of Landsat and Sentinel-2 -------------------
*/

```JavaScript

// -- Import the MODULES (functions)
var funcHLS = require('users/geeberra/Modules:HLS_Module_v2');

// -- Selected point of interest and clear sky (near)same dates Landsat-Sentinel observations
// OPTION 1 //// ######### Agriculture in Carazinho - Brazil
      var PoI = ee.Geometry.Point(-52.905090, -28.228550);//Point of Interest. Ground NDVI sensors in Carazinho
      var PointName = 'AgricultureCarazinho';
      //If precomputed TDOM image is in 16-bit scale, convert to 0-1 scale (/10000) to match the scale of the TDOM parameters
      var Pre_TDOM = ee.Image("users/geeberra/TDOM_Precomp/S2_mean_Std_Caraz_20152021").divide(10000);//Precomputed TDOM stats
      // - L8-S2 pair 
      var Date_L8 = ee.Date('2019-08-06');   
      var Date_S2_L8 = ee.Date('2019-08-06');
      // - L7-S2 pair
      var Date_L7 = ee.Date('2019-11-18');   
      var Date_S2_L7 = ee.Date('2019-11-19');

// OPTION 2 //// ######### Deciduous woodland - UK
      // var PoI = ee.Geometry.Point([-1.6961030702748747, 55.22054878343876]);//Point of Interest. 
      // var PointName = 'DeciduousWoodland';
      // // - L8-S2 pair 
      // var Date_L8 = ee.Date('2020-04-20');   
      // var Date_S2_L8 = ee.Date('2020-04-19');
      // // - L7-S2 pair
      // var Date_L7 = ee.Date('2020-04-21');   
      // var Date_S2_L7 = ee.Date('2020-04-19');
      
// OPTION 3 //// ######### Deciduous woodland - Oxford
      // var PoI = ee.Geometry.Point([-1.3334646332482958,51.77926643399076]);    
      // var PointName = 'Decid_Woodland_Oxford';
      // // - L8-S2 pair 
      // var Date_L8 = ee.Date('2020-06-25');   
      // var Date_S2_L8 = ee.Date('2020-06-25');
      // // - L7-S2 pair
      // var Date_L7 = ee.Date('2019-09-19');   
      // var Date_S2_L7 = ee.Date('2019-09-19');      

// Map.centerObject(PoI, 5);

// -- Choose which bands to get statistics from
var bandNames_img1_L7 = ee.List(['B3', 'B4']);//Red and NIR: L7
var bandNames_img1_L8 = ee.List(['B4', 'B5']);//Red and NIR: L8
var bandNames_img2_S2 = ee.List(['B4', 'B8']);//Red and NIR: S2
var bandNames_common = ee.List(['red', 'nir']);//Harmonized names

// -- Visualization parameters in an object literal.
var Color_comp_01 = {bands:"B4,B3,B2", min: 0.0, max: 0.4, gamma: 1};
var Color_comp_01_L7 = {bands:"B3,B2,B1", min: 0.0, max: 0.4, gamma: 1};
var Color_comp =    {bands:"B4,B3,B2", min:200, max:2000, gamma: 1};
var viz =  {min:200,max:2000,bands:['red','green','blue']};
var viz2 =  {min:0,max:0.4,bands:['red','green','blue']};


//Identify the equivalent bands across the three sensors
    //For L7:     B1=blue, B2=green, B3=red, B4=nir, B5=swir1, B7=swir2  
    var bandIn_L7 = ['B1',   'B2',    'B3',   'B4',   'B5',    'B7'];
    
    //For L8:      B2=blue, B3=green, B4=red, B5=nir, B6=swir1, B7=swir2
    var bandIn_L8 = ['B2',   'B3',    'B4',    'B5',   'B6',   'B7'];
    
    //For S2:    B2=blue, B3=green, B4=red, B8=nir, B11=swir1, B12=swir2  
    var bandIn_S2 = [  'B2',  'B3',  'B4',    'B8',    'B11',    'B12'];
    
    //BandOut will have the same name for all the sensors
    var bandOut = ['blue','green','red','nir','swir1','swir2'];
    

// ######################################################################################################
// ######################################################################################################
//                                    ### STEP 01: GENERAL FILTERING  ###
// ######################################################################################################

// -- Region of Interest (ROI) and time span 
// var polygon = PoI.buffer(300).bounds();//buffer around the point
var start_date = '2016-01-01';
var end_date   = '2020-12-31';


var criteria = ee.Filter.and(
    ee.Filter.bounds(PoI), ee.Filter.date(start_date, end_date));
var cloud_perc = 60;//Max cloud percentile per scene.    

// -- Collections of Landsat 7, 8 and Sentinel 2
var L8_col = ee.ImageCollection("LANDSAT/LC08/C02/T1_TOA")   //--Collection-2              
            // ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA')  //--Collection 1 >> T1_TOA or T1_SR
                .filter(criteria)
                .filter(ee.Filter.lt('CLOUD_COVER', cloud_perc));

var L7_col = ee.ImageCollection("LANDSAT/LE07/C02/T1_TOA") //--Collection-2
              // ee.ImageCollection('LANDSAT/LE07/C01/T1_TOA')//--Collection 1 >> T1_TOA or T1_SR
                .filter(criteria)
                .filter(ee.Filter.lt('CLOUD_COVER', cloud_perc));

var S2_col = ee.ImageCollection("COPERNICUS/S2_HARMONIZED")
            // ee.ImageCollection("COPERNICUS/S2")//S2 is TOA; S2_SR is BOA
                .filter(criteria)
                .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloud_perc));
/*S2 may have mulitple observations in the same day on some dates. 'funcHLS.FirstDate' selects one scene per date,
but works only for a single path/row every time. It's quite slow. Think wheter to use it.
  */
// S2_col = funcHLS.FirstDate(S2_col);
var footprint_S2 = S2_col.first().geometry().bounds();
// var footprint_S2 = geometry2;
Map.addLayer(footprint_S2, [], 'footprint_S2')

// ######################################################################################################
// ######################################################################################################
//                                ### STEP 02.CLOUD and SHADOW MASKING  ###
// ######################################################################################################

// Tuning parameters for S2 cloud and cloud shadow detection and masking
var args = {'Sentinel2Col': S2_col,
    'FilterCriteria': criteria,//For cloud masking
    'cloud_prob': 50,//For S2: Cloud probability threshold ranging in between 0-100%
    'dilatePixels': 1, //Pixels to dilate around clouds
    'contractPixels': 2,//Pixels to reduce cloud mask and dark shadows by to reduce inclusion of single-pixel comission errors
    'zShadowThresh': -1,// zShadowThresh: Threshold for cloud shadow masking- lower number masks out less. Between -0.8 and -1.2 generally works well
    'irSumThresh': 0.35,// irSumThresh: Sum of IR bands to include as shadows within TDOM  (lower number masks out less)
           //If precomputed TDOM image is in 16-bit scale, convert to 0-1 scale (/10000) to match the scale of the TDOM parameters
    'PrecomputedTDOMstats': Pre_TDOM//Precomputed TDOM stats
    // 'PrecomputedTDOMstats': null//'null' computes TDOM stats on the fly (slower). The code assumes S2 colletions are in the 16-bit scale.
};
var L7_col_masked = L7_col.map(funcHLS.cloudMaskL78_TOAoriginal);
var L8_col_masked = L8_col.map(funcHLS.cloudMaskL78_TOAoriginal);
// var S2_col_masked = funcHLS.Cloud_ClShadowS2(args);
var S2_col_masked = funcHLS.Cloud_ScorePlusS2(args);

// -Visualization and Checking
// var List = 1;
// var orimg = ee.Image(S2_col.toList(100).get(List));
//   Map.addLayer(orimg, Color_comp, 'Original');
// var selimg = ee.Image(S2_col_masked.toList(100).get(List));
//   Map.addLayer(selimg, Color_comp, 'Masked');




// ######################################################################################################
// ######################################################################################################
//                                ### STEP 03. INTER-SENSOR BAND ADJUSTMENT   ###
// ######################################################################################################


//Adjust L7 and L8 to match the S2 TOA spectral bands
var L7_col_adj = L7_col_masked.map(funcHLS.band_adjsut_L7);
var L8_col_adj = L8_col_masked.map(funcHLS.band_adjsut_L8);
var S2_col_adj = S2_col_masked;//Rename the S2 collection just to keep consistency with the L7-8 naming

// -Visualization and Checking
// var orimg = ee.Image(L8_col.toList(100).get(9))
//   Map.addLayer(orimg, Color_comp_01, 'Original');
// var selimg = ee.Image(L8_col_adj.toList(100).get(9))
//   Map.addLayer(selimg, Color_comp_01, 'Adjusted');





// ######################################################################################################
// ######################################################################################################
//                                    ### STEP 04.ATMOSPHERIC CORRECTION  ###
// ######################################################################################################


//Retrieve bottom of atmosphere (BOA) reflectance via SIAC
//Give common names to all bands
var L7_col_boa = L7_col_adj.map(funcHLS.SIAC_L7).select(bandIn_L7, bandOut);
var L8_col_boa = L8_col_adj.map(funcHLS.SIAC_L8).select(bandIn_L8, bandOut); 
var S2_col_boa = S2_col_adj.map(funcHLS.SIAC_S2).select(bandIn_S2, bandOut);


// -Visualization and Checking
// var orimg = ee.Image(S2_col.toList(100).get(9))
//   Map.addLayer(orimg, Color_comp_01, 'Original');
// var selimg = ee.Image(S2_col_boa.toList(100).get(9))
//   Map.addLayer(selimg, viz2, 'Boa');




// ######################################################################################################
// ######################################################################################################
//                              ### STEP 05. BRDF  ###
// ######################################################################################################
/*Code from <https://github.com/ndminhhus/geeguide/blob/master/04.topo_correction.md>
In a first test, the topographic correction seemed very computing power demanding. 
Think whether is really necessary to apply it or not.
*/

//Apply BRDF
var imgL7_SR_BRDF = L7_col_boa.map(funcHLS.applyBRDF); 
var imgL8_SR_BRDF = L8_col_boa.map(funcHLS.applyBRDF); 
var imgS2_SR_BRDF = S2_col_boa.map(funcHLS.applyBRDF);

// -Visualization and Checking
// var orimg = ee.Image(S2_col_boa.toList(100).get(9))
//   Map.addLayer(orimg, viz2, 'Original');
// var selimg = ee.Image(imgL8_SR_BRDF.toList(100).get(9))
//   Map.addLayer(selimg, viz2, 'BRDF');






// ######################################################################################################
// ######################################################################################################
//                                ### STEP 06. SPATIAL CO-REGISTRATION  ###
// ######################################################################################################


// 1st) Select the cloud-free image pairs.
// 2nd) Choose which band to calculate the displacement from.
var Bd = 'nir';
var Cloud_free_L7_Bd = funcHLS.selDate(L7_col_masked.select(bandIn_L7, bandOut), Date_L7).select(Bd);
var Cloud_free_S2_L7_Bd = funcHLS.selDate(S2_col_masked.select(bandIn_S2, bandOut), Date_S2_L7).select(Bd);

var Cloud_free_L8_Bd = funcHLS.selDate(L8_col_masked.select(bandIn_L8, bandOut), Date_L8).select(Bd);
var Cloud_free_S2_L8_Bd = funcHLS.selDate(S2_col_masked.select(bandIn_S2, bandOut), Date_S2_L8).select(Bd);


var iterpolation_method = "bilinear";//The interpolation mode to use: 'nearest_neighbor', 'bilinear' or 'bicubic'.)
var L7_regist_col = imgL7_SR_BRDF.map(funcHLS.Co_reg(Cloud_free_L7_Bd, Cloud_free_S2_L7_Bd, iterpolation_method));
var L8_regist_col = imgL8_SR_BRDF.map(funcHLS.Co_reg(Cloud_free_L8_Bd, Cloud_free_S2_L8_Bd, iterpolation_method));
var S2_regist_col = imgS2_SR_BRDF;//Rename the S2 collection just to keep consistency with the L7-8 naming

//-- Visualization and Checking --
// var orimg = ee.Image(L8_col.toList(100).get(0))
//   Map.addLayer(orimg, Color_comp_01, 'Original', false);
// var selimg = ee.Image(L8_regist_col.toList(100).get(0))
//   Map.addLayer(selimg, viz2, 'Registered');





// ######################################################################################################
// ######################################################################################################
//                                ### STEP 07. REPROJECT and RESAMPLE  ###
// ######################################################################################################



// ########## OPTION 1 ##########
      //1st: reproject S2 data to 30 m scale
      // var output_Scale = 30;//scale in meters
      // var imgS2_30 = S2_regist_col.map(funcHLS.reproj_S2(output_Scale));
      // var selected_S2_crs_L = imgS2_30.first().select('red').projection();

      // //2nd: Reproject and rescale L7-L8 images according to the S2 30m grid.
      // //No resample() was set: the default nearest-neighbor is used  to compute pixels in the chosen projection 
      // var imgL7_30 = L7_regist_col.map(funcHLS.reproj_L78(selected_S2_crs_L));
      // var imgL8_30 = L8_regist_col.map(funcHLS.reproj_L78(selected_S2_crs_L));
      

// ########## OPTION 2 ##########
/*Note: After tests, as described in the code 'HLS_test_2'(https://code.earthengine.google.com/f9375eab80b8479c15c773522a058497), this option seems
to result in least uncertainty as plot analysis revealed ('ImagePairs_Correl_St06_07.ipynb', https://colab.research.google.com/drive/15hv46w82cKoD9-3hdTr4gsqCxHlfn0zt?usp=drive_link)
*/
        //1st: Select the Landsat crs to reproject and re-scale the S2 data to
        var selected_L_crs_S2 = L8_regist_col.first().select('red').projection();
        //2nd: Reproject and rescale S2 images according to the L8 30m grid.
        // var imgS2_30 = S2_regist_col.map(funcHLS.resample_L78_30(selected_L_crs_S2));
        var imgS2_30 = S2_regist_col.map(funcHLS.reproj_L78(selected_L_crs_S2));
        var imgL7_30 = L7_regist_col;
        var imgL8_30 = L8_regist_col;


//-- Visualization and Checking --
// var N = 0;
// var orimg = ee.Image(L8_col.toList(100).get(N))
//   Map.addLayer(orimg, Color_comp_01, 'Original', false);
// viz2 =  {min:0,max:0.4,bands:['nir']};  
// var selimg = ee.Image(imgL8_30.select('nir').toList(100).get(N))
//   Map.addLayer(selimg, viz2, 'Reprojected');





// // ######################################################################################################
// // ######################################################################################################
// //                                    ### STEP 08: PRINT & EXPORT COLLECTIONS ###
// // ######################################################################################################


var s_date = start_date;
var e_date = end_date;

var L7_NDVI = imgL7_30.filterDate(s_date, e_date).map(funcHLS.calc_NDVI)
              .map(function(img){return img.set({'SATELLITE': 'LANDSAT_7'})});//Add 'SATELLITE' for plotting purposes  
var L8_NDVI = imgL8_30.filterDate(s_date, e_date).map(funcHLS.calc_NDVI)
              .map(function(img){return img.set({'SATELLITE': 'LANDSAT_8'})});
var S2_NDVI = imgS2_30.filterDate(s_date, e_date).map(funcHLS.calc_NDVI)
              .map(function(img){return img.set({'SATELLITE': 'SENTINEL_2'})});


// //Merge the series
// //<https://gis.stackexchange.com/questions/354961/plotting-time-series-of-different-products-in-google-earth-engine>
var Merged_Series = L7_NDVI.merge(L8_NDVI).merge(S2_NDVI);
Merged_Series = ee.ImageCollection(Merged_Series);
print('Merged_Series', Merged_Series)

// Calculate median NDVI for pixels intersecting the AOI for
// each image in the collection. Add the value as an image property.
var VI = 'NDVI';
var DateID = 'DATE';
var Sat_Names = 'SATELLITE';

// var allObs = Merged_Series.map(function(img) {
//   var obs = img.reduceRegion({
//     geometry: PoI,
//     reducer: ee.Reducer.median(),
//     scale: 30
//   });
//   var date_t = (ee.Date(img.date())).format("YYYY-MM-dd");
//   return img.set(VI, obs.get(VI)).set(DateID, date_t);
// });
// print('allObs', allObs)

// // // Make a chart of all observations where color distinguishes sensor.
// var chartAllObs = ui.Chart.feature.groups(
//   allObs, 'system:time_start', VI, Sat_Names)
//   .setChartType('ScatterChart')
//   .setOptions({
//     title: 'All Observations',
//     colors: ['red', 'green', 'blue'],
//     hAxis: {title: 'Date'},
//     vAxis: {title: 'NDVI'},
//     pointSize: 6,
//     dataOpacity: 0.5
//   });
// print(chartAllObs);


// Export.table.toDrive({
//   collection: allObs,
//   description: 'Drive_HLS_at_AOI',
//   folder: 'NDVI_fromHLS',
//   fileNamePrefix: ee.String('HLS').cat('_').cat(VI).cat('_').cat(PointName).cat('_')
//                   .cat(start_date.toString()).cat('_').cat(end_date.toString()).getInfo(),
//   fileFormat: 'CSV',
//   selectors: ([Sat_Names, DateID, VI])//Which columns to export
// }); 


var chart_merged = ui.Chart.image.series(Merged_Series, PoI, ee.Reducer.mean(), 30).setChartType('ScatterChart');
print('Merged_Series', chart_merged);


// // //
Map.addLayer(PoI, {}, 'Point of Interest');

```
