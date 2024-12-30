#-------------------------------------------------------------------------------
# Name:        Threats to GDEs due to human physical impacts
# Purpose:     Determine whether GDEs are at risk of damage from people
#              Roads and housing density - redundant?
#              Surface water diversions from NDWR database - PODs?
#              
#
# Author:      sarah.byer
#
# Created:     October 2021

#-------------------------------------------------------------------------------

# Import ArcGIS modules and check out spatial analyst extension
import arcpy, os
import pandas as pd
import numpy as np
from arcpy import env
from arcpy.sa import *
arcpy.CheckOutExtension("spatial")

# Path to temporary geodatabase
path =  r"E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Maps\GDE_Threats.gdb"

# Environment settings
env.workspace = path
env.overwriteOutput = True
env.outputCoordinateSystem = arcpy.SpatialReference(26911) # Spatial reference NAD 1983 UTM Zone 11N. The code is '26911'

# GDE layers
springs = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Springs'
wet = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Wetlands'
lp = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Lakes_Playas'
rs = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Rivers_Streams'
phr2 = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Phreatophytes_Explode'

#-------------------------------------------------------------------------------
# Surface diversions == stressor
# 1 if within half mile of gde
# 0.1 if not

# Subset surface diversions from POD dataset
pod = arcpy.CopyFeatures_management(r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NDWR\POD_101719.shp', 'pods_surface')
[f.name for f in arcpy.ListFields(pod)]
arcpy.GetCount_management(pod)

surface = ['STO', 'STR', 'SPR', 'LAK', 'RES', 'OSW'] # storage, stream, spring, lake, rservoir, other surface water
with arcpy.da.UpdateCursor(pod, 'source') as cursor:
    for row in cursor:
        if (row[0] not in surface):
            print(row[0])
            cursor.deleteRow()
del cursor
arcpy.GetCount_management(pod)

# Same as water level trend data - buffer by a half-mile and overlay GDEs
pod_buff = arcpy.Buffer_analysis(pod, 'pods_surface_halfmile', '0.5 Mile', '', '', 'ALL')
pod_buff = r'pods_surface_halfmile'
[f.name for f in arcpy.ListFields(pod_buff)]

# Select and attribute GDEs, as above
# Springs
arcpy.AddField_management(springs, 'SurfaceDiv', 'DOUBLE')
arcpy.CalculateField_management(springs, 'SurfaceDiv', 0.1)
[f.name for f in arcpy.ListFields(springs)]
spr_select = arcpy.SelectLayerByLocation_management(springs, 'INTERSECT', pod_buff)
arcpy.CalculateField_management(spr_select, 'SurfaceDiv', 1)

# Phreatophytes
arcpy.AddField_management(phr2, 'SurfaceDiv', 'DOUBLE')
arcpy.CalculateField_management(phr2, 'SurfaceDiv', 0.1)
[f.name for f in arcpy.ListFields(phr2)]
phr_select = arcpy.SelectLayerByLocation_management(phr2, 'INTERSECT', pod_buff)
arcpy.CalculateField_management(phr_select, 'SurfaceDiv', 1)

# Wetlands
arcpy.AddField_management(wet, 'SurfaceDiv', 'DOUBLE')
arcpy.CalculateField_management(wet, 'SurfaceDiv', 0.1)
[f.name for f in arcpy.ListFields(wet)]
wet_select = arcpy.SelectLayerByLocation_management(wet, 'INTERSECT', pod_buff)
arcpy.CalculateField_management(wet_select, 'SurfaceDiv', 1)

# Rivers/Streams
arcpy.AddField_management(rs, 'SurfaceDiv', 'DOUBLE')
arcpy.CalculateField_management(rs, 'SurfaceDiv', 0.1)
[f.name for f in arcpy.ListFields(rs)]
rs_select = arcpy.SelectLayerByLocation_management(rs, 'INTERSECT', pod_buff)
arcpy.CalculateField_management(rs_select, 'SurfaceDiv', 1)

# Lakes/Playas
arcpy.AddField_management(lp, 'SurfaceDiv', 'DOUBLE')
arcpy.CalculateField_management(lp, 'SurfaceDiv', 0.1)
[f.name for f in arcpy.ListFields(lp)]
lp_select = arcpy.SelectLayerByLocation_management(lp, 'INTERSECT', pod_buff)
arcpy.CalculateField_management(lp_select, 'SurfaceDiv', 1)

#-------------------------------------------------------------------------------
# Road and housing density
#d = arcpy.Raster(r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\road_density_tiger2014.tif')
#dmax = d.maximum
#d_norm = d/dmax # divide by max value to normalize btw 0 and 1
#d_norm.save('road_density_norm')
#roads = arcpy.CopyRaster_management(d_norm, r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\road_density_norm.tif')

# Housing Density
hd_now = arcpy.Raster(r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\HousingDensity\density_theobald_mos_2010.tif')
hd_future = arcpy.Raster(r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\HousingDensity\density_theobald_2060.tif')

hd_now_norm = hd_now/hd_now.maximum
hd_future_norm = hd_future/hd_future.maximum

arcpy.CopyRaster_management(hd_now_norm, r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\HousingDensity\density_theobald_2010_norm.tif')
arcpy.CopyRaster_management(hd_future_norm, r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\HousingDensity\density_theobald_2060_norm.tif')

# Change in housing density
hd_change = hd_future - hd_now
arcpy.CopyRaster_management(hd_change, r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\HousingDensity\density_theobald_2010to2060.tif')
# Positive values == increase in housing density
# Make negative values (places that decreased in density) NULL
x = SetNull(hd_change, hd_change, "Value < 1")
x.save('density_change_null')
hd_change_norm = x/x.maximum
arcpy.CopyRaster_management(hd_change_norm, r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\HousingDensity\density_theobald_change_norm.tif')

#-------------------------------------------------------------------------------
# Threat - projected change in housing density
# From normalized increase in housing density raster (above)

hd = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\HousingDensity\density_theobald_change_norm.tif'

# Springs
arcpy.DeleteField_management(springs, 'HousingThr')
temp_spr = arcpy.sa.ExtractValuesToPoints(springs, hd, 'springs_house_dens_change')
[f.name for f in arcpy.ListFields(temp_spr)]
arcpy.AddField_management(temp_spr, 'HousingThr', 'DOUBLE')
arcpy.CalculateField_management(temp_spr, 'HousingThr', '!RASTERVALU!', 'PYTHON3')
arcpy.JoinField_management(springs, 'OBJECTID', temp_spr, 'OBJECTID', ['HousingThr'])
[f.name for f in arcpy.ListFields(springs)]

# Wetlands
arcpy.DeleteField_management(wet, 'HousingThr')
temp_wet = arcpy.sa.ZonalStatisticsAsTable(wet, 'OBJECTID', hd, 'wetlands_house_denschange_tbl', 'DATA', 'MEAN')
[f.name for f in arcpy.ListFields(temp_wet)]
arcpy.AddField_management(temp_wet, 'HousingThr', 'DOUBLE')
arcpy.CalculateField_management(temp_wet, 'HousingThr', '!MEAN!', 'PYTHON3')
arcpy.JoinField_management(wet, 'OBJECTID', temp_wet, 'OBJECTID_1', ['HousingThr'])
[f.name for f in arcpy.ListFields(wet)]

# PHreatophytes
arcpy.DeleteField_management(phr2, 'HousingThr')
temp_phr = arcpy.sa.ZonalStatisticsAsTable(phr2, 'OBJECTID', hd, 'phr_house_denschange_tbl', 'DATA', 'MEAN')
[f.name for f in arcpy.ListFields(temp_phr)]
arcpy.AddField_management(temp_phr, 'HousingThr', 'DOUBLE')
arcpy.CalculateField_management(temp_phr, 'HousingThr', '!MEAN!', 'PYTHON3')
arcpy.JoinField_management(phr2, 'OBJECTID', temp_phr, 'OBJECTID_1', ['HousingThr'])
[f.name for f in arcpy.ListFields(phr2)]

# Lakes/playas
arcpy.DeleteField_management(lp, 'HousingThr')
temp_lp = arcpy.sa.ZonalStatisticsAsTable(lp, 'OBJECTID', hd, 'lp_house_denschange_tbl', 'DATA', 'MEAN')
[f.name for f in arcpy.ListFields(temp_lp)]
arcpy.AddField_management(temp_lp, 'HousingThr', 'DOUBLE')
arcpy.CalculateField_management(temp_lp, 'HousingThr', '!MEAN!', 'PYTHON3')
arcpy.JoinField_management(lp, 'OBJECTID', temp_lp, 'OBJECTID_1', ['HousingThr'])
[f.name for f in arcpy.ListFields(lp)]

# Rivers/streams
arcpy.DeleteField_management(rs, 'HousingThr')
temp_rs = arcpy.sa.ZonalStatisticsAsTable(rs, 'OBJECTID', hd, 'rs_house_denschange_tbl', 'DATA', 'MEAN')
[f.name for f in arcpy.ListFields(temp_rs)]
arcpy.AddField_management(temp_rs, 'HousingThr', 'DOUBLE')
arcpy.CalculateField_management(temp_rs, 'HousingThr', '!MEAN!', 'PYTHON3')
arcpy.JoinField_management(rs, 'OBJECTID', temp_rs, 'OBJECTID_1', ['HousingThr'])
[f.name for f in arcpy.ListFields(rs)]

# Done with this part


#-------------------------------------------------------------------------------
# Stress - current housing and/or road density
# Housing and road density are very similar... just do housing density as stressor for now

hd = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\HousingDensity\density_theobald_2010_norm.tif'

# Springs
temp_spr = arcpy.sa.ExtractValuesToPoints(springs, hd, 'springs_house_dens')
[f.name for f in arcpy.ListFields(temp_spr)]
arcpy.AddField_management(temp_spr, 'HousingStr', 'DOUBLE')
arcpy.CalculateField_management(temp_spr, 'HousingStr', '!RASTERVALU!', 'PYTHON3')
arcpy.JoinField_management(springs, 'OBJECTID', temp_spr, 'OBJECTID', ['HousingStr'])
[f.name for f in arcpy.ListFields(springs)]

# Wetlands
temp_wet = arcpy.sa.ZonalStatisticsAsTable(wet, 'OBJECTID', hd, 'wetlands_house_dens_tbl', 'DATA', 'MEAN')
[f.name for f in arcpy.ListFields(temp_wet)]
arcpy.AddField_management(temp_wet, 'HousingStr', 'DOUBLE')
arcpy.CalculateField_management(temp_wet, 'HousingStr', '!MEAN!', 'PYTHON3')
arcpy.JoinField_management(wet, 'OBJECTID', temp_wet, 'OBJECTID_1', ['HousingStr'])
[f.name for f in arcpy.ListFields(wet)]

# PHreatophytes
temp_phr = arcpy.sa.ZonalStatisticsAsTable(phr2, 'OBJECTID', hd, 'phr_house_dens_tbl', 'DATA', 'MEAN')
[f.name for f in arcpy.ListFields(temp_phr)]
arcpy.AddField_management(temp_phr, 'HousingStr', 'DOUBLE')
arcpy.CalculateField_management(temp_phr, 'HousingStr', '!MEAN!', 'PYTHON3')
arcpy.JoinField_management(phr2, 'OBJECTID', temp_phr, 'OBJECTID_1', ['HousingStr'])
[f.name for f in arcpy.ListFields(phr2)]

# Lakes/playas
temp_lp = arcpy.sa.ZonalStatisticsAsTable(lp, 'OBJECTID', hd, 'lp_house_dens_tbl', 'DATA', 'MEAN')
[f.name for f in arcpy.ListFields(temp_lp)]
arcpy.AddField_management(temp_lp, 'HousingStr', 'DOUBLE')
arcpy.CalculateField_management(temp_lp, 'HousingStr', '!MEAN!', 'PYTHON3')
arcpy.JoinField_management(lp, 'OBJECTID', temp_lp, 'OBJECTID_1', ['HousingStr'])
[f.name for f in arcpy.ListFields(lp)]

# Rivers/streams
temp_rs = arcpy.sa.ZonalStatisticsAsTable(rs, 'OBJECTID', hd, 'rs_house_dens_tbl', 'DATA', 'MEAN')
[f.name for f in arcpy.ListFields(temp_rs)]
arcpy.AddField_management(temp_rs, 'HousingStr', 'DOUBLE')
arcpy.CalculateField_management(temp_rs, 'HousingStr', '!MEAN!', 'PYTHON3')
arcpy.JoinField_management(rs, 'OBJECTID', temp_rs, 'OBJECTID_1', ['HousingStr'])
[f.name for f in arcpy.ListFields(rs)]


##### Not sure if using housing AND roads density would be redundant ####
#
## Roads just because
#roads = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\road_density_norm.tif'
#
#### CHECK TO MAKE SURE NULL VALUES AREN'T MESSING WITH THE SUMMARY STATS ###
#
## Springs
#temp_spr = arcpy.sa.ExtractValuesToPoints(springs, roads, 'springs_rds_dens')
#[f.name for f in arcpy.ListFields(temp_spr)]
#arcpy.AddField_management(temp_spr, 'Roads_Str', 'DOUBLE')
#arcpy.CalculateField_management(temp_spr, 'Roads_Str', '!RASTERVALU!', 'PYTHON3')
#arcpy.JoinField_management(springs, 'OBJECTID', temp_spr, 'OBJECTID', ['Roads_Str'])
#[f.name for f in arcpy.ListFields(springs)]
#
## Wetlands
#temp_wet = arcpy.sa.ZonalStatisticsAsTable(wet, 'OBJECTID', roads, 'wetlands_rds_dens_tbl', 'DATA', 'MEAN')
#[f.name for f in arcpy.ListFields(temp_wet)]
#arcpy.AddField_management(temp_wet, 'Roads_Str', 'DOUBLE')
#arcpy.CalculateField_management(temp_wet, 'Roads_Str', '!MEAN!', 'PYTHON3')
#arcpy.JoinField_management(wet, 'OBJECTID', temp_wet, 'OBJECTID_1', ['Roads_Str'])
#[f.name for f in arcpy.ListFields(wet)]
#
## PHreatophytes
#temp_phr = arcpy.sa.ZonalStatisticsAsTable(phr2, 'OBJECTID', roads, 'phr_rds_dens_tbl', 'DATA', 'MEAN')
#[f.name for f in arcpy.ListFields(temp_phr)]
#arcpy.AddField_management(temp_phr, 'Roads_Str', 'DOUBLE')
#arcpy.CalculateField_management(temp_phr, 'Roads_Str', '!MEAN!', 'PYTHON3')
#arcpy.JoinField_management(phr2, 'OBJECTID', temp_phr, 'OBJECTID_1', ['Roads_Str'])
#[f.name for f in arcpy.ListFields(phr2)]
#
## Lakes/playas
#temp_lp = arcpy.sa.ZonalStatisticsAsTable(lp, 'OBJECTID', roads, 'lp_rds_dens_tbl', 'DATA', 'MEAN')
#[f.name for f in arcpy.ListFields(temp_lp)]
#arcpy.AddField_management(temp_lp, 'Roads_Str', 'DOUBLE')
#arcpy.CalculateField_management(temp_lp, 'Roads_Str', '!MEAN!', 'PYTHON3')
#arcpy.JoinField_management(lp, 'OBJECTID', temp_lp, 'OBJECTID_1', ['Roads_Str'])
#[f.name for f in arcpy.ListFields(lp)]
#
## Rivers/streams
#temp_rs = arcpy.sa.ZonalStatisticsAsTable(rs, 'OBJECTID', roads, 'rs_rds_dens_tbl', 'DATA', 'MEAN')
#[f.name for f in arcpy.ListFields(temp_rs)]
#arcpy.AddField_management(temp_rs, 'Roads_Str', 'DOUBLE')
#arcpy.CalculateField_management(temp_rs, 'Roads_Str', '!MEAN!', 'PYTHON3')
#arcpy.JoinField_management(rs, 'OBJECTID', temp_rs, 'OBJECTID_1', ['Roads_Str'])
#[f.name for f in arcpy.ListFields(rs)]
#
#
## Done
#
