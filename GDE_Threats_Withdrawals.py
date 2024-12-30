#-------------------------------------------------------------------------------
# Name:        Consolidate well level trends and water rights data around GDEs
# Purpose:     Water level data from NWIS (USGS) and WellNet (NDWR), water rights data from NDWR hydro abstracts
#
# Author:      sarah.byer
#
# Created:     October 2021

#-------------------------------------------------------------------------------

# Import ArcGIS modules and check out spatial analyst extension
import arcpy, os
from arcpy import env
from arcpy.sa import *
arcpy.CheckOutExtension("spatial")

# Path to temporary geodatabase
path =  r"E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Maps\GDE_Threats.gdb"

# Environment settings
env.workspace = path
env.overwriteOutput = True
env.outputCoordinateSystem = arcpy.SpatialReference(26911) # Spatial reference NAD 1983 UTM Zone 11N. The code is '26911'


# Groundwater levels
gw = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\gwlevels_stats_110321.shp'

# Groundwater levels - 20 years of pre-irrigation data
gw = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\gwlevels_newrules_051021.shp'

# Pumping - water rights data
#wr = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\ndwr_abstracts_keepapps.shp'

# appropriation/pumping status
# Wilson 2020: https://tnc.box.com/s/ns19ih9cka3g4at8qrrcs4j314mmp74c
basins = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NDWR\hydrographic_basin_boundaries.shp'
[f.name for f in arcpy.ListFields(basins)]

# Copy GDE database layers to a new GDB - attribute these ones
gde2 = arcpy.Copy_management(r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_022120.gdb', 
                             r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb')
arcpy.ListFeatureClasses(r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb')

# GDE layers
springs = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Springs'
wet = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Wetlands'
phr = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Phreatophytes' #**need to explode**#
lp = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Lakes_Playas'
rs = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Rivers_Streams'

# Explode Phretophytes so each one can be attributed separately
phr2 = arcpy.MultipartToSinglepart_management(phr, r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Phreatophytes_Explode')
arcpy.GetCount_management(phr2)

phr2 = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Phreatophytes_Explode'
[f.name for f in arcpy.ListFields(phr2)]

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Coding:

#------------------------------------------------------------------------------- 
# Threats to GDEs = if GDE falls within half-mile of 50ft below-surface (shallow) groundwater
# Lopes et al 2006 shallow gw buffered by a half-mile
shallow = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Maps\GDE_Threats.gdb\shallow_gw_dissolve'

# For each GDE type - attribute whether it falls within half-mile of shallow gw
# Springs
arcpy.AddField_management(springs, 'Shallow_GW', 'LONG')
[f.name for f in arcpy.ListFields(springs)]
spr_select = arcpy.SelectLayerByLocation_management(springs, 'INTERSECT', shallow)
arcpy.CalculateField_management(spr_select, 'Shallow_GW', 1)

# Wetlands
arcpy.AddField_management(wet, 'Shallow_GW', 'LONG')
[f.name for f in arcpy.ListFields(wet)]
wet_select = arcpy.SelectLayerByLocation_management(wet, 'INTERSECT', shallow)
arcpy.CalculateField_management(wet_select, 'Shallow_GW', 1)

# Phreatophytes
arcpy.AddField_management(phr2, 'Shallow_GW', 'LONG')
[f.name for f in arcpy.ListFields(phr2)]
phr_select = arcpy.SelectLayerByLocation_management(phr2, 'INTERSECT', shallow)
arcpy.CalculateField_management(phr_select, 'Shallow_GW', 1)

# Lakes/playas
arcpy.AddField_management(lp, 'Shallow_GW', 'LONG')
[f.name for f in arcpy.ListFields(lp)]
lp_select = arcpy.SelectLayerByLocation_management(lp, 'INTERSECT', shallow)
arcpy.CalculateField_management(lp_select, 'Shallow_GW', 1)

# Rivers/streams
arcpy.AddField_management(rs, 'Shallow_GW', 'LONG')
[f.name for f in arcpy.ListFields(rs)]
rs_select = arcpy.SelectLayerByLocation_management(rs, 'INTERSECT', shallow)
arcpy.CalculateField_management(rs_select, 'Shallow_GW', 1)


#-------------------------------------------------------------------------------
# Stressor for GDE = if significant falling trend observed "near" GDE
# Does GDE have to land on the well? What if well lands on No GDEs?
# "Near" == within half mile? Like the threat?
# "Significant" falling trend == New (corrected) p-value <= 0.05 and negative Sens_Slope value

# Need to fix sens-slope values at very negative locations at Ruby Lake
# Change in vertical datum drastically reduced water level below surface
[f.name for f in arcpy.ListFields(gw)]

# Half mile buffer around falling gw levels
gw_fall = arcpy.SelectLayerByAttribute_management(gw, 'NEW_SELECTION', "SigFall = 1")
arcpy.GetCount_management(gw_fall)
gw_buff = arcpy.Buffer_analysis(gw_fall, 'gw_falling_halfmile', '0.5 Mile', '', '', 'ALL')

gw_buff = r'gw_falling_halfmile'
arcpy.GetCount_management(gw_buff)

# Select and attribute GDEs, as above
# Springs
arcpy.AddField_management(springs, 'Falling_GW', 'LONG')
[f.name for f in arcpy.ListFields(springs)]
spr_select = arcpy.SelectLayerByLocation_management(springs, 'INTERSECT', gw_buff)
arcpy.CalculateField_management(spr_select, 'Falling_GW', 1)

# Wetlands
arcpy.AddField_management(wet, 'Falling_GW', 'LONG')
[f.name for f in arcpy.ListFields(wet)]
wet_select = arcpy.SelectLayerByLocation_management(wet, 'INTERSECT', gw_buff)
arcpy.CalculateField_management(wet_select, 'Falling_GW', 1)

# Phreatophytes
arcpy.AddField_management(phr2, 'Falling_GW', 'LONG')
[f.name for f in arcpy.ListFields(phr2)]
phr_select = arcpy.SelectLayerByLocation_management(phr2, 'INTERSECT', gw_buff)
arcpy.CalculateField_management(phr_select, 'Falling_GW', 1)

# Lakes/playas
arcpy.AddField_management(lp, 'Falling_GW', 'LONG')
[f.name for f in arcpy.ListFields(lp)]
lp_select = arcpy.SelectLayerByLocation_management(lp, 'INTERSECT', gw_buff)
arcpy.CalculateField_management(lp_select, 'Falling_GW', 1)

# Rivers/streams
arcpy.AddField_management(rs, 'Falling_GW', 'LONG')
[f.name for f in arcpy.ListFields(rs)]
rs_select = arcpy.SelectLayerByLocation_management(rs, 'INTERSECT', gw_buff)
arcpy.CalculateField_management(rs_select, 'Falling_GW', 1)

# Not attributing with specific stressor/threat data
# I.e. Can tell if a GDE is stressed, but don't know the exact trend of the well it is getting stressed by
# Need to go to source data (NDWR/NWIS trend sites) for that information

#-------------------------------------------------------------------------------
# Threat or stress to GDEs from overcommitment or pumping
# Stressor to GDEs from over pumping

# Select and attribute GDEs, as above
# Springs
spr_commit = arcpy.SpatialJoin_analysis(springs, basins, 'temp_springs_join')
[f.name for f in arcpy.ListFields(spr_commit)]
arcpy.GetCount_management(spr_commit)
arcpy.GetCount_management(springs)

arcpy.AddField_management(spr_commit, 'CommitThr', 'DOUBLE')
with arcpy.da.UpdateCursor(spr_commit, ['OverCommit', 'CommitThr']) as cursor:
    for row in cursor:
        if row[0] == 100:
            row[1] = 0.5
            cursor.updateRow(row)
        elif row[0] == 200:
            row[1] = 1
            cursor.updateRow(row)
        elif row[0] == 0:
            row[1] = 0
            cursor.updateRow(row)
del cursor

arcpy.AddField_management(spr_commit, 'PumpStr', 'DOUBLE')
arcpy.CalculateField_management(spr_commit, 'PumpStr', '!OverPump!', 'PYTHON3')

arcpy.JoinField_management(springs, 'OBJECTID', spr_commit, 'OBJECTID', ['CommitThr', 'PumpStr']) # Also joining overpump (stress)
[f.name for f in arcpy.ListFields(springs)]
# Manually deleted/edited some springs that fall just on the border or outside the basin boundaries


# Phreatophytes
phr_commit = arcpy.SpatialJoin_analysis(phr2, basins, 'temp_phr_join')
[f.name for f in arcpy.ListFields(phr_commit)]
arcpy.GetCount_management(phr_commit)
arcpy.GetCount_management(phr2)

arcpy.AddField_management(phr_commit, 'CommitThr', 'DOUBLE')
with arcpy.da.UpdateCursor(phr_commit, ['OverCommit', 'CommitThr']) as cursor:
    for row in cursor:
        if row[0] == 100:
            row[1] = 0.5
            cursor.updateRow(row)
        elif row[0] == 200:
            row[1] = 1
            cursor.updateRow(row)
        elif row[0] == 0:
            row[1] = 0
            cursor.updateRow(row)
del cursor

arcpy.AddField_management(phr_commit, 'PumpStr', 'DOUBLE')
arcpy.CalculateField_management(phr_commit, 'PumpStr', '!OverPump!', 'PYTHON3')

arcpy.JoinField_management(phr2, 'OBJECTID', phr_commit, 'OBJECTID', ['CommitThr', 'PumpStr']) # Also joining overpump (stress)
[f.name for f in arcpy.ListFields(phr2)]


# Wetlands
wet_commit = arcpy.SpatialJoin_analysis(wet, basins, 'temp_wet_join')
[f.name for f in arcpy.ListFields(wet_commit)]
arcpy.GetCount_management(wet_commit)
arcpy.GetCount_management(wet)

arcpy.AddField_management(wet_commit, 'CommitThr', 'DOUBLE')
with arcpy.da.UpdateCursor(wet_commit, ['OverCommit', 'CommitThr']) as cursor:
    for row in cursor:
        if row[0] == 100:
            row[1] = 0.5
            cursor.updateRow(row)
        elif row[0] == 200:
            row[1] = 1
            cursor.updateRow(row)
        elif row[0] == 0:
            row[1] = 0
            cursor.updateRow(row)
del cursor

arcpy.AddField_management(wet_commit, 'PumpStr', 'DOUBLE')
arcpy.CalculateField_management(wet_commit, 'PumpStr', '!OverPump!', 'PYTHON3')

arcpy.JoinField_management(wet, 'OBJECTID', wet_commit, 'OBJECTID', ['CommitThr', 'PumpStr']) # Also joining overpump (stress)
[f.name for f in arcpy.ListFields(wet)]


# Rivers/streams
rs_commit = arcpy.SpatialJoin_analysis(rs, basins, 'temp_rs_join')
[f.name for f in arcpy.ListFields(rs_commit)]
arcpy.GetCount_management(rs_commit)
arcpy.GetCount_management(rs)

arcpy.AddField_management(rs_commit, 'CommitThr', 'DOUBLE')
with arcpy.da.UpdateCursor(rs_commit, ['OverCommit', 'CommitThr']) as cursor:
    for row in cursor:
        if row[0] == 100:
            row[1] = 0.5
            cursor.updateRow(row)
        elif row[0] == 200:
            row[1] = 1
            cursor.updateRow(row)
        elif row[0] == 0:
            row[1] = 0
            cursor.updateRow(row)
del cursor

arcpy.AddField_management(rs_commit, 'PumpStr', 'DOUBLE')
arcpy.CalculateField_management(rs_commit, 'PumpStr', '!OverPump!', 'PYTHON3')

arcpy.JoinField_management(rs, 'OBJECTID', rs_commit, 'OBJECTID', ['CommitThr', 'PumpStr']) # Also joining overpump (stress)
[f.name for f in arcpy.ListFields(rs)]


# Lakes/playas
lp_commit = arcpy.SpatialJoin_analysis(lp, basins, 'temp_lp_join')
[f.name for f in arcpy.ListFields(lp_commit)]
arcpy.GetCount_management(lp_commit)
arcpy.GetCount_management(lp)

arcpy.AddField_management(lp_commit, 'CommitThr', 'DOUBLE')
with arcpy.da.UpdateCursor(lp_commit, ['OverCommit', 'CommitThr']) as cursor:
    for row in cursor:
        if row[0] == 100:
            row[1] = 0.5
            cursor.updateRow(row)
        elif row[0] == 200:
            row[1] = 1
            cursor.updateRow(row)
        elif row[0] == 0:
            row[1] = 0
            cursor.updateRow(row)
del cursor

arcpy.AddField_management(lp_commit, 'PumpStr', 'DOUBLE')
arcpy.CalculateField_management(lp_commit, 'PumpStr', '!OverPump!', 'PYTHON3')

arcpy.JoinField_management(lp, 'OBJECTID', lp_commit, 'OBJECTID', ['CommitThr', 'PumpStr']) # Also joining overpump (stress)
[f.name for f in arcpy.ListFields(lp)]
