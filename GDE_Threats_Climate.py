#-------------------------------------------------------------------------------
# Name:        Climate stressors and threats to GDEs
# Purpose:     How does climate present a threat to GDEs? May need to wait for formula from Kevin and Louis
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
# GDE-specific threats
# Threat == 1 if GDE in local recharge area
# Threat == 0 if GDE not in local recharge area

# Recharge areas digitized from Mifflin 1988, Plate 3
rch = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\USGS\GDE_Recharge06172020.shp'
arcpy.GetCount_management(rch)

# For each GDE type - attribute whether it intersects with a recharge area
# Springs
arcpy.AddField_management(springs, 'Recharge', 'LONG')
arcpy.CalculateField_management(springs, 'Recharge', 0)
[f.name for f in arcpy.ListFields(springs)]
spr_select = arcpy.SelectLayerByLocation_management(springs, 'INTERSECT', rch)
arcpy.CalculateField_management(spr_select, 'Recharge', 1)

# Wetlands
arcpy.AddField_management(wet, 'Recharge', 'LONG')
arcpy.CalculateField_management(wet, 'Recharge', 0)
[f.name for f in arcpy.ListFields(wet)]
wet_select = arcpy.SelectLayerByLocation_management(wet, 'INTERSECT', rch)
arcpy.CalculateField_management(wet_select, 'Recharge', 1)

# Phreatophytes
arcpy.AddField_management(phr2, 'Recharge', 'LONG')
arcpy.CalculateField_management(phr2, 'Recharge', 0)
[f.name for f in arcpy.ListFields(phr2)]
phr_select = arcpy.SelectLayerByLocation_management(phr2, 'INTERSECT', rch)
arcpy.CalculateField_management(phr_select, 'Recharge', 1)

# Lakes/playas
arcpy.AddField_management(lp, 'Recharge', 'LONG')
arcpy.CalculateField_management(lp, 'Recharge', 0)
[f.name for f in arcpy.ListFields(lp)]
lp_select = arcpy.SelectLayerByLocation_management(lp, 'INTERSECT', rch)
arcpy.CalculateField_management(lp_select, 'Recharge', 1)

# Rivers/streams
arcpy.AddField_management(rs, 'Recharge', 'LONG')
arcpy.CalculateField_management(rs, 'Recharge', 0)
[f.name for f in arcpy.ListFields(rs)]
rs_select = arcpy.SelectLayerByLocation_management(rs, 'INTERSECT', rch)
arcpy.CalculateField_management(rs_select, 'Recharge', 1)



#-------------------------------------------------------------------------------
# Basin-wide threats
# Attribute all GDEs in basin with that basin's relative climate threat value

clim = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\BasinClimateThr.shp'
[f.name for f in arcpy.ListFields(clim)]

# Springs
spr_clim = arcpy.SpatialJoin_analysis(springs, clim, 'temp_springs_join')
[f.name for f in arcpy.ListFields(spr_clim)]
arcpy.GetCount_management(spr_clim)
arcpy.GetCount_management(springs)

arcpy.JoinField_management(springs, 'OBJECTID', spr_clim, 'OBJECTID', ['ClmThrR'])
[f.name for f in arcpy.ListFields(springs)]

# Wetlands
wet_clim = arcpy.SpatialJoin_analysis(wet, clim, 'temp_wet_join')
[f.name for f in arcpy.ListFields(wet_clim)]
arcpy.GetCount_management(wet_clim)
arcpy.GetCount_management(wet)

arcpy.JoinField_management(wet, 'OBJECTID', wet_clim, 'OBJECTID', ['ClmThrR'])
[f.name for f in arcpy.ListFields(wet)]

# Phreatophytes
phr_clim = arcpy.SpatialJoin_analysis(phr2, clim, 'temp_phr_join')
[f.name for f in arcpy.ListFields(phr_clim)]
arcpy.GetCount_management(phr_clim)
arcpy.GetCount_management(phr2)

arcpy.JoinField_management(phr2, 'OBJECTID', phr_clim, 'OBJECTID', ['ClmThrR'])
[f.name for f in arcpy.ListFields(phr2)]

# rivers/streams
rs_clim = arcpy.SpatialJoin_analysis(rs, clim, 'temp_rs_join')
[f.name for f in arcpy.ListFields(rs_clim)]
arcpy.GetCount_management(rs_clim)
arcpy.GetCount_management(rs)

arcpy.JoinField_management(rs, 'OBJECTID', rs_clim, 'OBJECTID', ['ClmThrR'])
[f.name for f in arcpy.ListFields(rs)]

# Lakes/playas
lp_clim = arcpy.SpatialJoin_analysis(lp, clim, 'temp_lp_join')
[f.name for f in arcpy.ListFields(lp_clim)]
arcpy.GetCount_management(lp_clim)
arcpy.GetCount_management(lp)

arcpy.JoinField_management(lp, 'OBJECTID', lp_clim, 'OBJECTID', ['ClmThrR'])
[f.name for f in arcpy.ListFields(lp)]

# Done

#-------------------------------------------------------------------------------
# Cannot determine stressors since climate change is purely future/potential