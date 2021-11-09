#%%

import os
import geopandas as gpd
import ee, eemont, geemap
import numpy as np
import pandas as pd
import rasterio as rio
import rioxarray as rxr
import xarray as xr
from datetime import datetime, timedelta

Map = geemap.Map()

def run_full_model(geopandas_object):

    '''
    This model is set up to use the 'MODIS/006/MCD43A4'sensor. The sensor can be changed by editing the string 
    'MODIS/006/MCD43A4' and also the references to the spectral bands. For example if the sensor is changed to 
    Landsat 8 then Change 'Nadir_Reflectance_Band1' to SR_4
    
    '''

    gdf_all = geopandas_object.to_crs('EPSG:4326')

    model_list = []

    for i in gdf_all.index:

        gdf = gdf_all.loc[[i]]
      
        fire_coords = gdf.total_bounds
        bbox_ee = ee.Geometry.BBox(fire_coords[0], fire_coords[1], fire_coords[2], fire_coords[3])

        gdf['fih_date1'] = datetime.strptime(gdf['fih_date1'].values[0], '%Y-%m-%d')

        gdf['fih_date_pre'] = gdf['fih_date1'] + timedelta(days = -14)
        gdf['fih_date_pre_end'] = gdf['fih_date1'] + timedelta(days = -1)

        gdf['fih_date_post'] = gdf['fih_date1'] + timedelta(days = 14)
        gdf['fih_date_post_end'] = gdf['fih_date_post'] + timedelta(days = 14)

        gdf['fih_date_ndvi_120'] = gdf['fih_date1'] + timedelta(days = -120)
        gdf['fih_date_ndvi_end'] = gdf['fih_date1'] + timedelta(days = -106)

        pre = pd.to_datetime(gdf['fih_date_pre'].values[0]).strftime('%Y-%m-%d')
        pre_end = pd.to_datetime(gdf['fih_date_pre_end'].values[0]).strftime('%Y-%m-%d')
        post = pd.to_datetime(gdf['fih_date_post'].values[0]).strftime('%Y-%m-%d')
        post_end = pd.to_datetime(gdf['fih_date_post_end'].values[0]).strftime('%Y-%m-%d')
        ndvi_120 = pd.to_datetime(gdf['fih_date_ndvi_120'].values[0]).strftime('%Y-%m-%d')
        ndvi_end = pd.to_datetime(gdf['fih_date_ndvi_end'].values[0]).strftime('%Y-%m-%d')

        modis_pre = (ee.ImageCollection('MODIS/006/MCD43A4')
         
            .filterBounds(bbox_ee)
            .filterDate(pre, pre_end)
            .select(
                'Nadir_Reflectance_Band2',
                'Nadir_Reflectance_Band7' 
            )
            .first()
            )

        modis_post = (ee.ImageCollection('MODIS/006/MCD43A4')
                
            .filterBounds(bbox_ee)
            .filterDate(post, post_end)
            .select(
                'Nadir_Reflectance_Band2',
                'Nadir_Reflectance_Band7' 
            )
            .first()
            )

        modis_ndvi = (ee.ImageCollection('MODIS/006/MCD43A4')
                
            .filterBounds(bbox_ee)
            .filterDate(ndvi_120, ndvi_end)
            .select(
                'Nadir_Reflectance_Band1',
                'Nadir_Reflectance_Band2' 
            )
            .first()
            )

        modis_ndvi_end = (ee.ImageCollection('MODIS/006/MCD43A4')
                
            .filterBounds(bbox_ee)
            .filterDate(pre, pre_end)
            .select(
                'Nadir_Reflectance_Band1',
                'Nadir_Reflectance_Band2' 
            )
            .first()
        )

        pre_path = "imagery/pre_test.tif"

        try: 

            os.remove(pre_path)

        except: None

        geemap.ee_export_image(modis_pre, 
                               filename = pre_path, 
                               region = bbox_ee, 
                               scale = 500,
                               file_per_band = False)

        # Open and read in the MODIS file
        # Note that rxr is the alias for rioxarray
        raster_pre = rxr.open_rasterio(pre_path)

        pre_nbr = (raster_pre[0] - raster_pre[1]) / (raster_pre[0] + raster_pre[1])

        rxr.open_rasterio(pre_path).close()

        post_path = "imagery/post_test.tif"

        try: 

            os.remove(post_path)

        except: None

        geemap.ee_export_image(modis_post, 
                               filename = post_path, 
                               region = bbox_ee, 
                               scale = 500,
                               file_per_band = False)

        # Open and read in the MODIS file
        # Note that rxr is the alias for rioxarray
        raster_post = rxr.open_rasterio(post_path)

        post_nbr = (raster_post[0] - raster_post[1]) / (raster_post[0] + raster_post[1])

        rxr.open_rasterio(post_path).close()

        ndvi_path = "imagery/ndvi_test.tif"

        try: 

            os.remove(ndvi_path)

        except: None

        geemap.ee_export_image(modis_ndvi, 
                               filename = ndvi_path, 
                               region = bbox_ee, 
                               scale = 500,
                               file_per_band = False)

        # Open and read in the MODIS file
        # Note that rxr is the alias for rioxarray
        raster_ndvi = rxr.open_rasterio(ndvi_path)

        ndvi = (raster_ndvi[1] - raster_ndvi[0]) / (raster_ndvi[1] + raster_ndvi[0])

        rxr.open_rasterio(ndvi_path).close()

        post_ndvi_path = "imagery/post_ndvi_test.tif"

        try: 

            os.remove(post_ndvi_path)

        except: None

        geemap.ee_export_image(modis_ndvi_end, 
                               filename = post_ndvi_path, 
                               region = bbox_ee, 
                               scale = 500,
                               file_per_band = False)

        # Open and read in the MODIS file
        # Note that rxr is the alias for rioxarray
        raster_post_ndvi = rxr.open_rasterio(post_ndvi_path)

        post_ndvi = (raster_post_ndvi[1] - raster_post_ndvi[0]) / (raster_post_ndvi[1] + raster_post_ndvi[0])

        rxr.open_rasterio(post_ndvi_path).close()

        dnbr = (pre_nbr.values - post_nbr.values)
        dndvi = (post_ndvi.values - ndvi.values)

        dnbr_xr = xr.DataArray(dnbr, 
                                coords = [post_nbr.y.values, 
                                post_nbr.x.values], 
                                dims = ['y', 'x'])

        dndvi_xr = xr.DataArray(dndvi, 
                        coords = [post_ndvi.y.values, 
                                  post_ndvi.x.values], 
                        dims = ['y', 'x'])

        # Define dNBR classification bins
        dnbr_class_bins = [-np.inf, -.1, .1, .27, .66, np.inf]

        dnbr_class_xr = xr.apply_ufunc(np.digitize,
                                       dnbr_xr,
                                       dnbr_class_bins)

        dnbr_df = dnbr_xr.to_dataframe(name = 'dnbr').reset_index(drop = True)
        dnbr_class_df = dnbr_class_xr.to_dataframe(name = 'dnbr_class').reset_index(drop = True)
        
        pre_ndvi_df = ndvi.to_dataframe(name = 'ndvi_120').reset_index(drop = True)
        post_ndvi_df = post_ndvi.to_dataframe(name = 'ndvi_14').reset_index(drop = True)
        dndvi_df = dndvi_xr.to_dataframe(name = 'dndvi').reset_index(drop = True)

        reg_df = pd.concat([dnbr_df, 
                            dnbr_class_df, 
                            pre_ndvi_df,
                            post_ndvi_df,
                            dndvi_df], axis = 1)

        reg_df['fire_name'] = gdf['fih_number'].values[0]

        model_list.append(reg_df)

    return model_list
# %%
