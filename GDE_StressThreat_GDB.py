#-------------------------------------------------------------------------------
# Name:        Create new GDB for iGDEs with stressor and threat data
#              Build in metadata and field names/aliases
#              Also copy to shapefiles for individual layer distribution
#
# Author:      sarah.byer
#
# Created:     March 2022

#-------------------------------------------------------------------------------

# Import ArcGIS modules and check out spatial analyst extension
import arcpy, os
import pandas as pd
import numpy as np
from arcpy import env
from arcpy.sa import *
arcpy.CheckOutExtension("spatial")

# Path to temporary geodatabase
path = r"K:\GIS3\Projects\GDE\Maps\GDE_Threats\GDE_Threats.gdb"

# Environment settings
env.workspace = path
env.overwriteOutput = True
env.outputCoordinateSystem = arcpy.SpatialReference(26911) # Spatial reference NAD 1983 UTM Zone 11N. The code is '26911'

# Create copy of the OG iGDE database
# want to retain (as much as possible) og metadata and feature ids
#og_gde = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_050919.gdb'
og_gde = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_022120.gdb'
#gdes = arcpy.Copy_management(og_gde, r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322.gdb')
gdes = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322.gdb'
arcpy.ListFeatureClasses(gdes)

# database with stress/threat codes
stats = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_assess2.gdb'

##############################################################
 # For field defs and names, refer to E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\gde_field_defs_gdb.xlsx
 # Order, New_Field and New_Alias columns
##############################################################

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Springs
spr_og = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322.gdb\springs'
print(arcpy.GetCount_management(spr_og))
print([f.name for f in arcpy.ListFields(spr_og)])
spr_st = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_assess2.gdb\springs'
print(arcpy.GetCount_management(spr_st))
print([f.name for f in arcpy.ListFields(spr_st)])

# Good!

# Copy stats to table
# Join doesn't like drawing right from the fc
tbl = arcpy.TableToTable_conversion(spr_st, env.workspace, 'spring_stats_tbl')

# Before join - change order
#fields = ['Falling_GW', 'PumpStr', 'WW_Norm_Str', 'WW_Wgt_Str', 'Shallow_GW', 'CommitThr', 'WW_Norm_Thr', 'WW_Wgt_Thr',
#          'Recharge', 'Cl_Wgt_Str', 'RenormThr', 'Cl_Wgt_Thr',
#          'UnguStr', 'Ung_Wgt_Str', 'Ung_Wgt_Thr', 
#          'InvSpp_Str', 'NN_Wgt_Str', 'InvSpp_Thr', 'NN_Wgt_Thr',
#          'HousingStr', 'SurfaceDiv', 'HA_Norm_Str', 'HA_Wgt_Str', 'HousingThr', 'HA_Wgt_Thr',
#          'ScoreStr', 'ScoreThr']
fields_nowgt = ['Falling_GW', 'PumpStr', 'WW_Norm_Str', 'Shallow_GW', 'CommitThr', 'WW_Norm_Thr', 
          'Recharge', 'RenormThr', 
          'Ung_Wgt_Str', 'Ung_Wgt_Thr', 
          'NN_Wgt_Str', 'NN_Wgt_Thr',
          'HousingStr', 'SurfaceDiv', 'HA_Wgt_Str', 'HA_Wgt_Thr',
          'ScoreStr', 'ScoreThr']


arcpy.JoinField_management(spr_og, 'SPRING_ID', tbl, 'SPRING_ID', fields_nowgt)
[f.name for f in arcpy.ListFields(spr_og)]

# Last-minute, add intermediate grazing fields
ung_fields = ['GRAZE', 'HORSE_BURRO', 'ELK']
arcpy.JoinField_management(spr_og, 'SPRING_ID', tbl, 'SPRING_ID', ung_fields)
[f.name for f in arcpy.ListFields(spr_og)]

 # Alter fields - update field alias an maybe field name
 # Can use for other feature classes
fc = spr_og
arcpy.AlterField_management(fc, 'Falling_GW', 'Falling_GW', 'Falling Groundwater Level Score (groundwater withdrawal stress)')
arcpy.AlterField_management(fc, 'PumpStr', 'Pump_GW', 'Pumping Status Score (groundwater withdrawal stress)')
arcpy.AlterField_management(fc, 'WW_Norm_Str', 'GW_Str', 'Groundwater Withdrawal Stress Score')
arcpy.AlterField_management(fc, 'Shallow_GW', 'Shallow_GW', 'Shallow Groundwater Score (groundwater withdrawal threat)')
arcpy.AlterField_management(fc, 'CommitThr', 'Commit_GW', 'Appropriation Status Score (groundwater withdrawal threat)')
arcpy.AlterField_management(fc, 'WW_Norm_Thr', 'GW_Thr', 'Groundwater Withdrawal Threat Score')

arcpy.AlterField_management(fc, 'Recharge', 'Clim_Str', 'Climate Stress Score')
arcpy.AlterField_management(fc, 'RenormThr', 'Clim_Thr', 'Climate Threat Score')

arcpy.AlterField_management(fc, 'GRAZE', 'Ung_Graze', 'Cattle and Sheep Potential')
arcpy.AlterField_management(fc, 'HORSE_BURRO', 'Ung_WHB', 'Wild Horse and Burro Potential')
arcpy.AlterField_management(fc, 'ELK', 'Ung_Elk', 'Elk Potential')
arcpy.AlterField_management(fc, 'Ung_Wgt_Str', 'Ung_Str', 'Ungulate Stress Score')
arcpy.AlterField_management(fc, 'Ung_Wgt_Thr', 'Ung_Thr', 'Ungulate Threat Score')

arcpy.AlterField_management(fc, 'NN_Wgt_Str', 'NN_Str', 'Non-native Species Stress Score')
arcpy.AlterField_management(fc, 'NN_Wgt_Thr', 'NN_Thr', 'Non-native Species Threat Score')

arcpy.AlterField_management(fc, 'HousingStr', 'Hous_OHD', 'Housing Density Stress Score (other human development stress)')
arcpy.AlterField_management(fc, 'SurfaceDiv', 'SurfDiv_OHD', 'Surface Water Diversions (other human development stress)')
arcpy.AlterField_management(fc, 'HA_Wgt_Str', 'OHD_Str', 'Other Human Development Stress Score')
arcpy.AlterField_management(fc, 'HA_Wgt_Thr', 'OHD_Thr', 'Other Human Development Threat Score')

arcpy.AlterField_management(fc, 'ScoreStr', 'Score_Str', 'Overall Stress Score')
arcpy.AlterField_management(fc, 'ScoreThr', 'Score_Thr', 'Overall Threat Score')

[f.name for f in arcpy.ListFields(fc)]
# Need to manually move fields



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Wetlands
wet_og = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322.gdb\wetlands'
print(arcpy.GetCount_management(wet_og))
print([f.name for f in arcpy.ListFields(wet_og)])
wet_st = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_assess2.gdb\wetlands'
print(arcpy.GetCount_management(wet_st))
print([f.name for f in arcpy.ListFields(wet_st)])

# good!

# Copy stats to table
# Join doesn't like drawing right from the fc
tbl = arcpy.TableToTable_conversion(wet_st, env.workspace, 'wetland_stats_tbl')

fields_nowgt = ['Falling_GW', 'PumpStr', 'WW_Norm_Str', 'Shallow_GW', 'CommitThr', 'WW_Norm_Thr', 
          'Recharge', 'RenormThr', 
          'Ung_Wgt_Str', 'Ung_Wgt_Thr', 
          'NN_Wgt_Str', 'NN_Wgt_Thr',
          'HousingStr', 'SurfaceDiv', 'HA_Wgt_Str', 'HA_Wgt_Thr',
          'ScoreStr', 'ScoreThr']


arcpy.JoinField_management(wet_og, 'OBJECTID', tbl, 'OBJECTID', fields_nowgt)
[f.name for f in arcpy.ListFields(wet_og)]

# Last-minute, add intermediate grazing fields
tbl = 'wetland_stats_tbl'
ung_fields = ['GRAZE', 'HORSE_BURRO', 'ELK']
arcpy.JoinField_management(wet_og, 'OBJECTID', tbl, 'OBJECTID', ung_fields)
[f.name for f in arcpy.ListFields(wet_og)]

 # Alter fields - update field alias an maybe field name
 # Can use for other feature classes
fc = wet_og
arcpy.AlterField_management(fc, 'Falling_GW', 'Falling_GW', 'Falling Groundwater Level Score (groundwater withdrawal stress)')
arcpy.AlterField_management(fc, 'PumpStr', 'Pump_GW', 'Pumping Status Score (groundwater withdrawal stress)')
arcpy.AlterField_management(fc, 'WW_Norm_Str', 'GW_Str', 'Groundwater Withdrawal Stress Score')

arcpy.AlterField_management(fc, 'Shallow_GW', 'Shallow_GW', 'Shallow Groundwater Score (groundwater withdrawal threat)')
arcpy.AlterField_management(fc, 'CommitThr', 'Commit_GW', 'Appropriation Status Score (groundwater withdrawal threat)')
arcpy.AlterField_management(fc, 'WW_Norm_Thr', 'GW_Thr', 'Groundwater Withdrawal Threat Score')

arcpy.AlterField_management(fc, 'Recharge', 'Clim_Str', 'Climate Stress Score')
arcpy.AlterField_management(fc, 'RenormThr', 'Clim_Thr', 'Climate Threat Score')

arcpy.AlterField_management(fc, 'GRAZE', 'Ung_Graze', 'Cattle and Sheep Potential')
arcpy.AlterField_management(fc, 'HORSE_BURRO', 'Ung_WHB', 'Wild Horse and Burro Potential')
arcpy.AlterField_management(fc, 'ELK', 'Ung_Elk', 'Elk Potential')
arcpy.AlterField_management(fc, 'Ung_Wgt_Str', 'Ung_Str', 'Ungulate Stress Score')
arcpy.AlterField_management(fc, 'Ung_Wgt_Thr', 'Ung_Thr', 'Ungulate Threat Score')

arcpy.AlterField_management(fc, 'NN_Wgt_Str', 'NN_Str', 'Non-native Species Stress Score')
arcpy.AlterField_management(fc, 'NN_Wgt_Thr', 'NN_Thr', 'Non-native Species Threat Score')

arcpy.AlterField_management(fc, 'HousingStr', 'Hous_OHD', 'Housing Density Stress Score (other human development stress)')
arcpy.AlterField_management(fc, 'SurfaceDiv', 'SurfDiv_OHD', 'Surface Water Diversions (other human development stress)')
arcpy.AlterField_management(fc, 'HA_Wgt_Str', 'OHD_Str', 'Other Human Development Stress Score')
arcpy.AlterField_management(fc, 'HA_Wgt_Thr', 'OHD_Thr', 'Other Human Development Threat Score')

arcpy.AlterField_management(fc, 'ScoreStr', 'Score_Str', 'Overall Stress Score')
arcpy.AlterField_management(fc, 'ScoreThr', 'Score_Thr', 'Overall Threat Score')

[f.name for f in arcpy.ListFields(fc)]
# Need to manually move fields


# Wetlands
wet_og = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322.gdb\wetlands'
print(arcpy.GetCount_management(wet_og))
print([f.name for f in arcpy.ListFields(wet_og)])
wet_st = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_assess2.gdb\wetlands'
print(arcpy.GetCount_management(wet_st))
print([f.name for f in arcpy.ListFields(wet_st)])

# good!


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Rivers and Streams
rs_og = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322.gdb\Rivers_Streams'
print(arcpy.GetCount_management(rs_og))
print([f.name for f in arcpy.ListFields(rs_og)])
rs_st = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_assess2.gdb\Rivers_Streams'
print(arcpy.GetCount_management(rs_st))
print([f.name for f in arcpy.ListFields(rs_st)])

# 1 more river/stream in old version...
st_ids = list()
with arcpy.da.SearchCursor(rs_st, 'OBJECTID') as cursor:
    for row in cursor:
        st_ids.append(row[0])
del cursor

new_ids = list()
with arcpy.da.SearchCursor(rs_og, 'OBJECTID') as cursor:
    for row in cursor:
        if(row[0] not in st_ids):
            new_ids.append(row[0])
del cursor

print(new_ids)

# Copy stats to table
# Join doesn't like drawing right from the fc
tbl = arcpy.TableToTable_conversion(rs_st, env.workspace, 'riverstreams_stats_tbl')

fields_nowgt = ['Falling_GW', 'PumpStr', 'WW_Norm_Str', 'Shallow_GW', 'CommitThr', 'WW_Norm_Thr', 
          'Recharge', 'RenormThr', 
          'Ung_Wgt_Str', 'Ung_Wgt_Thr', 
          'NN_Wgt_Str', 'NN_Wgt_Thr',
          'HousingStr', 'SurfaceDiv', 'HA_Wgt_Str', 'HA_Wgt_Thr',
          'ScoreStr', 'ScoreThr']


arcpy.JoinField_management(rs_og, 'OBJECTID', tbl, 'OBJECTID', fields_nowgt)
[f.name for f in arcpy.ListFields(rs_og)]

# Last-minute, add intermediate grazing fields
tbl = 'riverstreams_stats_tbl'
ung_fields = ['GRAZE', 'HORSE_BURRO', 'ELK']
arcpy.JoinField_management(rs_og, 'OBJECTID', tbl, 'OBJECTID', ung_fields)
[f.name for f in arcpy.ListFields(rs_og)]

# Alter fields - update field alias an maybe field name
# Can use for other feature classes
fc = rs_og
arcpy.AlterField_management(fc, 'Falling_GW', 'Falling_GW', 'Falling Groundwater Level Score (groundwater withdrawal stress)')
arcpy.AlterField_management(fc, 'PumpStr', 'Pump_GW', 'Pumping Status Score (groundwater withdrawal stress)')
arcpy.AlterField_management(fc, 'WW_Norm_Str', 'GW_Str', 'Groundwater Withdrawal Stress Score')

arcpy.AlterField_management(fc, 'Shallow_GW', 'Shallow_GW', 'Shallow Groundwater Score (groundwater withdrawal threat)')
arcpy.AlterField_management(fc, 'CommitThr', 'Commit_GW', 'Appropriation Status Score (groundwater withdrawal threat)')
arcpy.AlterField_management(fc, 'WW_Norm_Thr', 'GW_Thr', 'Groundwater Withdrawal Threat Score')

arcpy.AlterField_management(fc, 'Recharge', 'Clim_Str', 'Climate Stress Score')
arcpy.AlterField_management(fc, 'RenormThr', 'Clim_Thr', 'Climate Threat Score')

arcpy.AlterField_management(fc, 'GRAZE', 'Ung_Graze', 'Cattle and Sheep Potential')
arcpy.AlterField_management(fc, 'HORSE_BURRO', 'Ung_WHB', 'Wild Horse and Burro Potential')
arcpy.AlterField_management(fc, 'ELK', 'Ung_Elk', 'Elk Potential')
arcpy.AlterField_management(fc, 'Ung_Wgt_Str', 'Ung_Str', 'Ungulate Stress Score')
arcpy.AlterField_management(fc, 'Ung_Wgt_Thr', 'Ung_Thr', 'Ungulate Threat Score')

arcpy.AlterField_management(fc, 'NN_Wgt_Str', 'NN_Str', 'Non-native Species Stress Score')
arcpy.AlterField_management(fc, 'NN_Wgt_Thr', 'NN_Thr', 'Non-native Species Threat Score')

arcpy.AlterField_management(fc, 'HousingStr', 'Hous_OHD', 'Housing Density Stress Score (other human development stress)')
arcpy.AlterField_management(fc, 'SurfaceDiv', 'SurfDiv_OHD', 'Surface Water Diversions (other human development stress)')
arcpy.AlterField_management(fc, 'HA_Wgt_Str', 'OHD_Str', 'Other Human Development Stress Score')
arcpy.AlterField_management(fc, 'HA_Wgt_Thr', 'OHD_Thr', 'Other Human Development Threat Score')

arcpy.AlterField_management(fc, 'ScoreStr', 'Score_Str', 'Overall Stress Score')
arcpy.AlterField_management(fc, 'ScoreThr', 'Score_Thr', 'Overall Threat Score')

[f.name for f in arcpy.ListFields(fc)]
# Need to manually move fields


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Lakes and Playas
lp_og = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322.gdb\Lakes_Playas'
print(arcpy.GetCount_management(lp_og))
print([f.name for f in arcpy.ListFields(lp_og)])
lp_st = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_assess2.gdb\Lakes_Playas'
print(arcpy.GetCount_management(lp_st))
print([f.name for f in arcpy.ListFields(lp_st)])

# Good!

# Copy stats to table
# Join doesn't like drawing right from the fc
tbl = arcpy.TableToTable_conversion(lp_st, env.workspace, 'lakesplayas_stats_tbl')

fields_nowgt = ['Falling_GW', 'PumpStr', 'WW_Norm_Str', 'Shallow_GW', 'CommitThr', 'WW_Norm_Thr', 
          'Recharge', 'RenormThr', 
          'Ung_Wgt_Str', 'Ung_Wgt_Thr', 
          'NN_Wgt_Str', 'NN_Wgt_Thr',
          'HousingStr', 'SurfaceDiv', 'HA_Wgt_Str', 'HA_Wgt_Thr',
          'ScoreStr', 'ScoreThr']


arcpy.JoinField_management(lp_og, 'OBJECTID', tbl, 'OBJECTID', fields_nowgt)
[f.name for f in arcpy.ListFields(lp_og)]

# Last-minute, add intermediate grazing fields
tbl = 'lakesplayas_stats_tbl'
# Need to create ungulate fields - they're all zero
arcpy.AddField_management(tbl, 'GRAZE', 'LONG')
arcpy.AddField_management(tbl, 'HORSE_BURRO', 'LONG')
arcpy.AddField_management(tbl, 'ELK', 'LONG')
ung_fields = ['GRAZE', 'HORSE_BURRO', 'ELK']
arcpy.CalculateField_management(tbl, ung_fields[0], 0, 'PYTHON3')
arcpy.CalculateField_management(tbl, ung_fields[1], 0, 'PYTHON3')
arcpy.CalculateField_management(tbl, ung_fields[2], 0, 'PYTHON3')

ung_fields = ['GRAZE', 'HORSE_BURRO', 'ELK']
arcpy.JoinField_management(lp_og, 'OBJECTID', tbl, 'OBJECTID', ung_fields)
[f.name for f in arcpy.ListFields(lp_og)]

# Alter fields - update field alias an maybe field name
# Can use for other feature classes
fc = lp_og
arcpy.AlterField_management(fc, 'Falling_GW', 'Falling_GW', 'Falling Groundwater Level Score (groundwater withdrawal stress)')
arcpy.AlterField_management(fc, 'PumpStr', 'Pump_GW', 'Pumping Status Score (groundwater withdrawal stress)')
arcpy.AlterField_management(fc, 'WW_Norm_Str', 'GW_Str', 'Groundwater Withdrawal Stress Score')

arcpy.AlterField_management(fc, 'Shallow_GW', 'Shallow_GW', 'Shallow Groundwater Score (groundwater withdrawal threat)')
arcpy.AlterField_management(fc, 'CommitThr', 'Commit_GW', 'Appropriation Status Score (groundwater withdrawal threat)')
arcpy.AlterField_management(fc, 'WW_Norm_Thr', 'GW_Thr', 'Groundwater Withdrawal Threat Score')

arcpy.AlterField_management(fc, 'Recharge', 'Clim_Str', 'Climate Stress Score')
arcpy.AlterField_management(fc, 'RenormThr', 'Clim_Thr', 'Climate Threat Score')

arcpy.AlterField_management(fc, 'GRAZE', 'Ung_Graze', 'Cattle and Sheep Potential')
arcpy.AlterField_management(fc, 'HORSE_BURRO', 'Ung_WHB', 'Wild Horse and Burro Potential')
arcpy.AlterField_management(fc, 'ELK', 'Ung_Elk', 'Elk Potential')
arcpy.AlterField_management(fc, 'Ung_Wgt_Str', 'Ung_Str', 'Ungulate Stress Score')
arcpy.AlterField_management(fc, 'Ung_Wgt_Thr', 'Ung_Thr', 'Ungulate Threat Score')

arcpy.AlterField_management(fc, 'NN_Wgt_Str', 'NN_Str', 'Non-native Species Stress Score')
arcpy.AlterField_management(fc, 'NN_Wgt_Thr', 'NN_Thr', 'Non-native Species Threat Score')

arcpy.AlterField_management(fc, 'HousingStr', 'Hous_OHD', 'Housing Density Stress Score (other human development stress)')
arcpy.AlterField_management(fc, 'SurfaceDiv', 'SurfDiv_OHD', 'Surface Water Diversions (other human development stress)')
arcpy.AlterField_management(fc, 'HA_Wgt_Str', 'OHD_Str', 'Other Human Development Stress Score')
arcpy.AlterField_management(fc, 'HA_Wgt_Thr', 'OHD_Thr', 'Other Human Development Threat Score')

arcpy.AlterField_management(fc, 'ScoreStr', 'Score_Str', 'Overall Stress Score')
arcpy.AlterField_management(fc, 'ScoreThr', 'Score_Thr', 'Overall Threat Score')

[f.name for f in arcpy.ListFields(fc)]
# Need to manually move fields


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Phreatophytes
# better to just replace in new GDB
# Exploded to calculate stress & threat scores
phr_og = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322.gdb\Phreatophytes'
print(arcpy.GetCount_management(phr_og))
print([f.name for f in arcpy.ListFields(phr_og)])
phr_st = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_assess2.gdb\Phreatophytes_Explode'
print(arcpy.GetCount_management(phr_st))
print([f.name for f in arcpy.ListFields(phr_st)])

# Delete Phreatophytes feature in new GDB
arcpy.Delete_management(phr_og)

# Copy exploded phreatophytes to new GDB
# Only copy over necessary fields
phr_fields = ['SOURCE_CODE', 'PHR_TYPE', 'PHR_GROUP', 'COMMENTS',
            'Falling_GW', 'PumpStr', 'WW_Norm_Str', 'Shallow_GW', 'CommitThr', 'WW_Norm_Thr', 
            'Recharge', 'RenormThr', 
            'Ung_Wgt_Str', 'Ung_Wgt_Thr', 
            'NN_Wgt_Str', 'NN_Wgt_Thr',
            'HousingStr', 'SurfaceDiv', 'HA_Wgt_Str', 'HA_Wgt_Thr',
            'ScoreStr', 'ScoreThr']

# Field mapping object
fmap = arcpy.FieldMappings()
for f in phr_fields:
    print(f)
    fm = arcpy.FieldMap()
    fm.addInputField(phr_st, f)
    fmap.addFieldMap(fm)
    
# Copy to GDB with only named fields - will take a minute
arcpy.FeatureClassToFeatureClass_conversion(phr_st, r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322.gdb', 'Phreatophytes',
                                            '', fmap)
[f.name for f in arcpy.ListFields(phr_st)]

# Last-minute, add intermediate grazing fields
phr_st2 = arcpy.CopyFeatures_management(r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_assess2.gdb\Phreatophytes_Explode', 
                                'K:\GIS3\Projects\GDE\Maps\GDE_Threats\GDE_Threats.gdb\temp_phr_explode')
[f.name for f in arcpy.ListFields(phr_st2)]
ung_fields = ['GRAZE', 'HORSE_BURR', 'ELK'] # Oops, typo on HORSE_BURRO
arcpy.JoinField_management(phr_og, 'OBJECTID', phr_st2, 'OBJECTID', ung_fields)
[f.name for f in arcpy.ListFields(phr_og)]

# Alter fields - update field alias an maybe field name
# Can use for other feature classes
phr_og = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322.gdb\Phreatophytes'
fc = phr_og
arcpy.AlterField_management(fc, 'Falling_GW', 'Falling_GW', 'Falling Groundwater Level Score (groundwater withdrawal stress)')
arcpy.AlterField_management(fc, 'PumpStr', 'Pump_GW', 'Pumping Status Score (groundwater withdrawal stress)')
arcpy.AlterField_management(fc, 'WW_Norm_Str', 'GW_Str', 'Groundwater Withdrawal Stress Score')

arcpy.AlterField_management(fc, 'Shallow_GW', 'Shallow_GW', 'Shallow Groundwater Score (groundwater withdrawal threat)')
arcpy.AlterField_management(fc, 'CommitThr', 'Commit_GW', 'Appropriation Status Score (groundwater withdrawal threat)')
arcpy.AlterField_management(fc, 'WW_Norm_Thr', 'GW_Thr', 'Groundwater Withdrawal Threat Score')

arcpy.AlterField_management(fc, 'Recharge', 'Clim_Str', 'Climate Stress Score')
arcpy.AlterField_management(fc, 'RenormThr', 'Clim_Thr', 'Climate Threat Score')

arcpy.AlterField_management(fc, 'GRAZE', 'Ung_Graze', 'Cattle and Sheep Potential')
arcpy.AlterField_management(fc, 'HORSE_BURR', 'Ung_WHB', 'Wild Horse and Burro Potential')
arcpy.AlterField_management(fc, 'ELK', 'Ung_Elk', 'Elk Potential')
arcpy.AlterField_management(fc, 'Ung_Wgt_Str', 'Ung_Str', 'Ungulate Stress Score')
arcpy.AlterField_management(fc, 'Ung_Wgt_Thr', 'Ung_Thr', 'Ungulate Threat Score')

arcpy.AlterField_management(fc, 'NN_Wgt_Str', 'NN_Str', 'Non-native Species Stress Score')
arcpy.AlterField_management(fc, 'NN_Wgt_Thr', 'NN_Thr', 'Non-native Species Threat Score')

arcpy.AlterField_management(fc, 'HousingStr', 'Hous_OHD', 'Housing Density Stress Score (other human development stress)')
arcpy.AlterField_management(fc, 'SurfaceDiv', 'SurfDiv_OHD', 'Surface Water Diversions (other human development stress)')
arcpy.AlterField_management(fc, 'HA_Wgt_Str', 'OHD_Str', 'Other Human Development Stress Score')
arcpy.AlterField_management(fc, 'HA_Wgt_Thr', 'OHD_Thr', 'Other Human Development Threat Score')

arcpy.AlterField_management(fc, 'ScoreStr', 'Score_Str', 'Overall Stress Score')
arcpy.AlterField_management(fc, 'ScoreThr', 'Score_Thr', 'Overall Threat Score')

[f.name for f in arcpy.ListFields(fc)]
# Need to manually move fields


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Copy to shapefiles

# Folder to write shapefiles to
shps = r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322_shps'

arcpy.ListFeatures(gdes)
arcpy.CopyFeatures_management(spr_og, r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322_shps\NV_iGDE_Springs.shp')
arcpy.CopyFeatures_management(wet_og, r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322_shps\NV_iGDE_Wetlands.shp')
arcpy.CopyFeatures_management(rs_og, r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322_shps\NV_iGDE_Rivers_Streams.shp')
arcpy.CopyFeatures_management(lp_og, r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322_shps\NV_iGDE_Lakes_Playas.shp')
arcpy.CopyFeatures_management(phr_og, r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322_shps\NV_iGDE_Phreatohpytes.shp')

arcpy.CopyFeatures_management('K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322.gdb\Species', 
                              r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322_shps\NV_iGDE_Species.shp')
arcpy.TableToTable_conversion('K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322.gdb\Species_tbl', 
                              r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322_shps', 'NV_iGDE_Species_tbl.dbf')
arcpy.TableToTable_conversion('K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322.gdb\Source_tbl', 
                              r'K:\GIS3\Projects\GDE\Geospatial\NV_iGDE_032322_shps', 'NV_iGDE_Source_tbl.dbf')







