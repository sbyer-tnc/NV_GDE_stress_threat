#-------------------------------------------------------------------------------
# Name:        Score stressors and threats to GDEs
# Purpose:     Sum up scores of stressors and threats to GDEs
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

# Function to fill rows with zeroes if values are null
def null2zero(feature, field):
    with arcpy.da.UpdateCursor(feature, field) as cursor:
        for row in cursor:
            if row[0] is None:
                row[0] = 0
                cursor.updateRow(row)
    del cursor
    return feature


# Function to get numeric maximum of any feature's fields
def getMax(feature, field):
    all_rows = [i[0] for i in arcpy.da.SearchCursor(feature, field)]
    max_val = max(all_rows)
    return max_val


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Springs
# Combine and normalize grouped scores
[f.name for f in arcpy.ListFields(springs)]
null2zero(springs, 'PumpStr')
null2zero(springs, 'CommitThr')
null2zero(springs, 'Shallow_GW')
null2zero(springs, 'Falling_GW')
null2zero(springs, 'Recharge')
null2zero(springs, 'ClmThrR')
null2zero(springs, 'InvSpp_Thr')
null2zero(springs, 'InvSpp_Str')
null2zero(springs, 'HousingThr')
null2zero(springs, 'HousingStr')

# Water withdrawal weighting - stress
max1 = getMax(springs, 'Falling_GW')
max2 = getMax(springs, 'PumpStr')
max_norm = max1 + max2
print(max_norm)
arcpy.AddField_management(springs, 'WW_Norm_Str', 'DOUBLE')
where = '(!Falling_GW! + !PumpStr!) / ' + "%s" %max_norm
arcpy.CalculateField_management(springs, 'WW_Norm_Str', where, 'PYTHON3')
arcpy.AddField_management(springs, 'WW_Wgt_Str', 'DOUBLE')
where = '(!WW_Norm_Str!) * 5'
arcpy.CalculateField_management(springs, 'WW_Wgt_Str', where, 'PYTHON3')

# Water withdrawal weighting - threat
max1 = getMax(springs, 'Shallow_GW')
max2 = getMax(springs, 'CommitThr')
max_norm = max1 + max2
print(max_norm)
arcpy.AddField_management(springs, 'WW_Norm_Thr', 'DOUBLE')
where = '(!Shallow_GW! + !CommitThr!) / ' + "%s" %max_norm
arcpy.CalculateField_management(springs, 'WW_Norm_Thr', where, 'PYTHON3') 
arcpy.AddField_management(springs, 'WW_Wgt_Thr', 'DOUBLE')
where = '(!WW_Norm_Thr!) * 5'
arcpy.CalculateField_management(springs, 'WW_Wgt_Thr', where, 'PYTHON3')

# Climate weighting - stress
arcpy.AddField_management(springs, 'Cl_Wgt_Str', 'DOUBLE')
where = '(!Recharge!) * 4'
arcpy.CalculateField_management(springs, 'Cl_Wgt_Str', where, 'PYTHON3')

# Climate weighting - threat
arcpy.AddField_management(springs, 'Cl_Wgt_Thr', 'DOUBLE')
where = '(!ClmThrR!) * 4'
arcpy.CalculateField_management(springs, 'Cl_Wgt_Thr', where, 'PYTHON3')

# Human activity weighting - stress
max1 = getMax(springs, 'HousingStr')
max2 = getMax(springs, 'SurfaceDiv')
max_norm = max1 + max2
print(max_norm)
arcpy.AddField_management(springs, 'HA_Norm_Str', 'DOUBLE')
where = '(!HousingStr! + !SurfaceDiv!) / ' + "%s" %max_norm
arcpy.CalculateField_management(springs, 'HA_Norm_Str', where, 'PYTHON3')
arcpy.AddField_management(springs, 'HA_Wgt_Str', 'DOUBLE')
where = '(!HA_Norm_Str!) * 3'
arcpy.CalculateField_management(springs, 'HA_Wgt_Str', where, 'PYTHON3')

# Human activity weighting - threat
arcpy.AddField_management(springs, 'HA_Wgt_Thr', 'DOUBLE')
where = '(!HousingThr!) * 3'
arcpy.CalculateField_management(springs, 'HA_Wgt_Thr', where, 'PYTHON3')

# Non-natice species weighting - stress
arcpy.AddField_management(springs, 'NN_Wgt_Str', 'DOUBLE')
where = '(!InvSpp_Str!) * 2'
arcpy.CalculateField_management(springs, 'NN_Wgt_Str', where, 'PYTHON3')

# Non-native species weighting - threat
arcpy.AddField_management(springs, 'NN_Wgt_Thr', 'DOUBLE')
where = '(!InvSpp_Thr!) * 2'
arcpy.CalculateField_management(springs, 'NN_Wgt_Thr', where, 'PYTHON3')

# Ungulate weighting - stress (same as threat)
arcpy.AddField_management(springs, 'Ung_Wgt_Str', 'DOUBLE')
where = '(!UnguStr!) * 1'
arcpy.CalculateField_management(springs, 'Ung_Wgt_Str', where, 'PYTHON3')

# Ungulate weighting - threat (same as stress)
arcpy.AddField_management(springs, 'Ung_Wgt_Thr', 'DOUBLE')
where = '(!UnguStr!) * 1'
arcpy.CalculateField_management(springs, 'Ung_Wgt_Thr', where, 'PYTHON3')


# Sum all stress and threat scores
[f.name for f in arcpy.ListFields(springs)]
arcpy.AddField_management(springs, 'ScoreStr', 'DOUBLE')
arcpy.CalculateField_management(springs, 'ScoreStr', '!WW_Wgt_Str! + !Cl_Wgt_Str! + !HA_Wgt_Str! + !NN_Wgt_Str! + !Ung_Wgt_Str!', 'PYTHON3')
arcpy.AddField_management(springs, 'ScoreThr', 'DOUBLE')
arcpy.CalculateField_management(springs, 'ScoreThr', '!WW_Wgt_Thr! + !Cl_Wgt_Thr! + !HA_Wgt_Thr! + !NN_Wgt_Thr! + !Ung_Wgt_Thr!', 'PYTHON3')


#arcpy.DeleteField_management(springs, ['ScoreStr', 'ScoreThr', 'WW_Wgt_Str', 'Cl_Wgt_Str', 'HA_Wgt_Str', 
#                                       'NN_Wgt_Str', 'Ung_Wgt_Str', 'WW_Wgt_Thr', 'Cl_Wgt_Thr', 'HA_Wgt_Thr', 
#                                       'NN_Wgt_Thr', 'Ung_Wgt_Thr', 'WW_Norm_Str', 'WW_Norm_Thr', 'HA_Norm_Str'])


#-------------------------------------------------------------------------------
# Wetlands
[f.name for f in arcpy.ListFields(wet)]
null2zero(wet, 'PumpStr')
null2zero(wet, 'CommitThr')
null2zero(wet, 'Shallow_GW')
null2zero(wet, 'Falling_GW')
null2zero(wet, 'Recharge')
null2zero(wet, 'ClmThrR')
null2zero(wet, 'InvSpp_Thr')
null2zero(wet, 'InvSpp_Str')
null2zero(wet, 'HousingThr')
null2zero(wet, 'HousingStr')

# Water withdrawal weighting - stress
max1 = getMax(wet, 'Falling_GW')
max2 = getMax(wet, 'PumpStr')
max_norm = max1 + max2
print(max_norm)
arcpy.AddField_management(wet, 'WW_Norm_Str', 'DOUBLE')
where = '(!Falling_GW! + !PumpStr!) / ' + "%s" %max_norm
arcpy.CalculateField_management(wet, 'WW_Norm_Str', where, 'PYTHON3')
arcpy.AddField_management(wet, 'WW_Wgt_Str', 'DOUBLE')
where = '(!WW_Norm_Str!) * 5'
arcpy.CalculateField_management(wet, 'WW_Wgt_Str', where, 'PYTHON3')

# Water withdrawal weighting - threat
max1 = getMax(wet, 'Shallow_GW')
max2 = getMax(wet, 'CommitThr')
max_norm = max1 + max2
print(max_norm)
arcpy.AddField_management(wet, 'WW_Norm_Thr', 'DOUBLE')
where = '(!Shallow_GW! + !CommitThr!) / ' + "%s" %max_norm
arcpy.CalculateField_management(wet, 'WW_Norm_Thr', where, 'PYTHON3') 
arcpy.AddField_management(wet, 'WW_Wgt_Thr', 'DOUBLE')
where = '(!WW_Norm_Thr!) * 5'
arcpy.CalculateField_management(wet, 'WW_Wgt_Thr', where, 'PYTHON3')

# Climate weighting - stress
arcpy.AddField_management(wet, 'Cl_Wgt_Str', 'DOUBLE')
where = '(!Recharge!) * 4'
arcpy.CalculateField_management(wet, 'Cl_Wgt_Str', where, 'PYTHON3')

# Climate weighting - threat
arcpy.AddField_management(wet, 'Cl_Wgt_Thr', 'DOUBLE')
where = '(!ClmThrR!) * 4'
arcpy.CalculateField_management(wet, 'Cl_Wgt_Thr', where, 'PYTHON3')

# Human activity weighting - stress
max1 = getMax(wet, 'HousingStr')
max2 = getMax(wet, 'SurfaceDiv')
max_norm = max1 + max2
print(max_norm)
arcpy.AddField_management(wet, 'HA_Norm_Str', 'DOUBLE')
where = '(!HousingStr! + !SurfaceDiv!) / ' + "%s" %max_norm
arcpy.CalculateField_management(wet, 'HA_Norm_Str', where, 'PYTHON3')
arcpy.AddField_management(wet, 'HA_Wgt_Str', 'DOUBLE')
where = '(!HA_Norm_Str!) * 3'
arcpy.CalculateField_management(wet, 'HA_Wgt_Str', where, 'PYTHON3')

# Human activity weighting - threat
arcpy.AddField_management(wet, 'HA_Wgt_Thr', 'DOUBLE')
where = '(!HousingThr!) * 3'
arcpy.CalculateField_management(wet, 'HA_Wgt_Thr', where, 'PYTHON3')

# Non-native species weighting - stress
arcpy.AddField_management(wet, 'NN_Wgt_Str', 'DOUBLE')
where = '(!InvSpp_Str!) * 2'
arcpy.CalculateField_management(wet, 'NN_Wgt_Str', where, 'PYTHON3')

# Non-native species weighting - threat
arcpy.AddField_management(wet, 'NN_Wgt_Thr', 'DOUBLE')
where = '(!InvSpp_Thr!) * 2'
arcpy.CalculateField_management(wet, 'NN_Wgt_Thr', where, 'PYTHON3')

# Ungulate weighting - stress (same as threat)
arcpy.AddField_management(wet, 'Ung_Wgt_Str', 'DOUBLE')
where = '(!UnguStr!) * 1'
arcpy.CalculateField_management(wet, 'Ung_Wgt_Str', where, 'PYTHON3')

# Ungulate weighting - threat (same as stress)
arcpy.AddField_management(wet, 'Ung_Wgt_Thr', 'DOUBLE')
where = '(!UnguStr!) * 1'
arcpy.CalculateField_management(wet, 'Ung_Wgt_Thr', where, 'PYTHON3')


# Sum all stress and threat scores
[f.name for f in arcpy.ListFields(wet)]
arcpy.AddField_management(wet, 'ScoreStr', 'DOUBLE')
arcpy.CalculateField_management(wet, 'ScoreStr', '!WW_Wgt_Str! + !Cl_Wgt_Str! + !HA_Wgt_Str! + !NN_Wgt_Str! + !Ung_Wgt_Str!', 'PYTHON3')
arcpy.AddField_management(wet, 'ScoreThr', 'DOUBLE')
arcpy.CalculateField_management(wet, 'ScoreThr', '!WW_Wgt_Thr! + !Cl_Wgt_Thr! + !HA_Wgt_Thr! + !NN_Wgt_Thr! + !Ung_Wgt_Thr!', 'PYTHON3')




#-------------------------------------------------------------------------------
# Phreatophytes
[f.name for f in arcpy.ListFields(phr2)]
null2zero(phr2, 'PumpStr')
null2zero(phr2, 'CommitThr')
null2zero(phr2, 'Shallow_GW')
null2zero(phr2, 'Falling_GW')
null2zero(phr2, 'Recharge')
null2zero(phr2, 'ClmThrR')
null2zero(phr2, 'InvSpp_Thr')
null2zero(phr2, 'InvSpp_Str')
null2zero(phr2, 'HousingThr')
null2zero(phr2, 'HousingStr')

# Water withdrawal weighting - stress
max1 = getMax(phr2, 'Falling_GW')
max2 = getMax(phr2, 'PumpStr')
max_norm = max1 + max2
print(max_norm)
arcpy.AddField_management(phr2, 'WW_Norm_Str', 'DOUBLE')
where = '(!Falling_GW! + !PumpStr!) / ' + "%s" %max_norm
arcpy.CalculateField_management(phr2, 'WW_Norm_Str', where, 'PYTHON3')
arcpy.AddField_management(phr2, 'WW_Wgt_Str', 'DOUBLE')
where = '(!WW_Norm_Str!) * 5'
arcpy.CalculateField_management(phr2, 'WW_Wgt_Str', where, 'PYTHON3')

# Water withdrawal weighting - threat
max1 = getMax(phr2, 'Shallow_GW')
max2 = getMax(phr2, 'CommitThr')
max_norm = max1 + max2
print(max_norm)
arcpy.AddField_management(phr2, 'WW_Norm_Thr', 'DOUBLE')
where = '(!Shallow_GW! + !CommitThr!) / ' + "%s" %max_norm
arcpy.CalculateField_management(phr2, 'WW_Norm_Thr', where, 'PYTHON3') 
arcpy.AddField_management(phr2, 'WW_Wgt_Thr', 'DOUBLE')
where = '(!WW_Norm_Thr!) * 5'
arcpy.CalculateField_management(phr2, 'WW_Wgt_Thr', where, 'PYTHON3')

# Climate weighting - stress
arcpy.AddField_management(phr2, 'Cl_Wgt_Str', 'DOUBLE')
where = '(!Recharge!) * 4'
arcpy.CalculateField_management(phr2, 'Cl_Wgt_Str', where, 'PYTHON3')

# Climate weighting - threat
arcpy.AddField_management(phr2, 'Cl_Wgt_Thr', 'DOUBLE')
where = '(!ClmThrR!) * 4'
arcpy.CalculateField_management(phr2, 'Cl_Wgt_Thr', where, 'PYTHON3')

# Human activity weighting - stress
max1 = getMax(phr2, 'HousingStr')
max2 = getMax(phr2, 'SurfaceDiv')
max_norm = max1 + max2
print(max_norm)
arcpy.AddField_management(phr2, 'HA_Norm_Str', 'DOUBLE')
where = '(!HousingStr! + !SurfaceDiv!) / ' + "%s" %max_norm
arcpy.CalculateField_management(phr2, 'HA_Norm_Str', where, 'PYTHON3')
arcpy.AddField_management(phr2, 'HA_Wgt_Str', 'DOUBLE')
where = '(!HA_Norm_Str!) * 3'
arcpy.CalculateField_management(phr2, 'HA_Wgt_Str', where, 'PYTHON3')

# Human activity weighting - threat
arcpy.AddField_management(phr2, 'HA_Wgt_Thr', 'DOUBLE')
where = '(!HousingThr!) * 3'
arcpy.CalculateField_management(phr2, 'HA_Wgt_Thr', where, 'PYTHON3')

# Non-native species weighting - stress
arcpy.AddField_management(phr2, 'NN_Wgt_Str', 'DOUBLE')
where = '(!InvSpp_Str!) * 2'
arcpy.CalculateField_management(phr2, 'NN_Wgt_Str', where, 'PYTHON3')

# Non-native species weighting - threat
arcpy.AddField_management(phr2, 'NN_Wgt_Thr', 'DOUBLE')
where = '(!InvSpp_Thr!) * 2'
arcpy.CalculateField_management(phr2, 'NN_Wgt_Thr', where, 'PYTHON3')

# Ungulate weighting - stress (same as threat)
arcpy.AddField_management(phr2, 'Ung_Wgt_Str', 'DOUBLE')
where = '(!UnguStr!) * 1'
arcpy.CalculateField_management(phr2, 'Ung_Wgt_Str', where, 'PYTHON3')

# Ungulate weighting - threat (same as stress)
arcpy.AddField_management(phr2, 'Ung_Wgt_Thr', 'DOUBLE')
where = '(!UnguStr!) * 1'
arcpy.CalculateField_management(phr2, 'Ung_Wgt_Thr', where, 'PYTHON3')


# Sum all stress and threat scores
[f.name for f in arcpy.ListFields(phr2)]
arcpy.AddField_management(phr2, 'ScoreStr', 'DOUBLE')
arcpy.CalculateField_management(phr2, 'ScoreStr', '!WW_Wgt_Str! + !Cl_Wgt_Str! + !HA_Wgt_Str! + !NN_Wgt_Str! + !Ung_Wgt_Str!', 'PYTHON3')
arcpy.AddField_management(phr2, 'ScoreThr', 'DOUBLE')
arcpy.CalculateField_management(phr2, 'ScoreThr', '!WW_Wgt_Thr! + !Cl_Wgt_Thr! + !HA_Wgt_Thr! + !NN_Wgt_Thr! + !Ung_Wgt_Thr!', 'PYTHON3')




#-------------------------------------------------------------------------------
# Rivers/streams
[f.name for f in arcpy.ListFields(rs)]
null2zero(rs, 'PumpStr')
null2zero(rs, 'CommitThr')
null2zero(rs, 'Shallow_GW')
null2zero(rs, 'Falling_GW')
null2zero(rs, 'Recharge')
null2zero(rs, 'ClmThrR')
null2zero(rs, 'InvSpp_Thr')
null2zero(rs, 'InvSpp_Str')
null2zero(rs, 'HousingThr')
null2zero(rs, 'HousingStr')

# Water withdrawal weighting - stress
max1 = getMax(rs, 'Falling_GW')
max2 = getMax(rs, 'PumpStr')
max_norm = max1 + max2
print(max_norm)
arcpy.AddField_management(rs, 'WW_Norm_Str', 'DOUBLE')
where = '(!Falling_GW! + !PumpStr!) / ' + "%s" %max_norm
arcpy.CalculateField_management(rs, 'WW_Norm_Str', where, 'PYTHON3')
arcpy.AddField_management(rs, 'WW_Wgt_Str', 'DOUBLE')
where = '(!WW_Norm_Str!) * 5'
arcpy.CalculateField_management(rs, 'WW_Wgt_Str', where, 'PYTHON3')

# Water withdrawal weighting - threat
max1 = getMax(rs, 'Shallow_GW')
max2 = getMax(rs, 'CommitThr')
max_norm = max1 + max2
print(max_norm)
arcpy.AddField_management(rs, 'WW_Norm_Thr', 'DOUBLE')
where = '(!Shallow_GW! + !CommitThr!) / ' + "%s" %max_norm
arcpy.CalculateField_management(rs, 'WW_Norm_Thr', where, 'PYTHON3') 
arcpy.AddField_management(rs, 'WW_Wgt_Thr', 'DOUBLE')
where = '(!WW_Norm_Thr!) * 5'
arcpy.CalculateField_management(rs, 'WW_Wgt_Thr', where, 'PYTHON3')

# Climate weighting - stress
arcpy.AddField_management(rs, 'Cl_Wgt_Str', 'DOUBLE')
where = '(!Recharge!) * 4'
arcpy.CalculateField_management(rs, 'Cl_Wgt_Str', where, 'PYTHON3')

# Climate weighting - threat
arcpy.AddField_management(rs, 'Cl_Wgt_Thr', 'DOUBLE')
where = '(!ClmThrR!) * 4'
arcpy.CalculateField_management(rs, 'Cl_Wgt_Thr', where, 'PYTHON3')

# Human activity weighting - stress
max1 = getMax(rs, 'HousingStr')
max2 = getMax(rs, 'SurfaceDiv')
max_norm = max1 + max2
print(max_norm)
arcpy.AddField_management(rs, 'HA_Norm_Str', 'DOUBLE')
where = '(!HousingStr! + !SurfaceDiv!) / ' + "%s" %max_norm
arcpy.CalculateField_management(rs, 'HA_Norm_Str', where, 'PYTHON3')
arcpy.AddField_management(rs, 'HA_Wgt_Str', 'DOUBLE')
where = '(!HA_Norm_Str!) * 3'
arcpy.CalculateField_management(rs, 'HA_Wgt_Str', where, 'PYTHON3')

# Human activity weighting - threat
arcpy.AddField_management(rs, 'HA_Wgt_Thr', 'DOUBLE')
where = '(!HousingThr!) * 3'
arcpy.CalculateField_management(rs, 'HA_Wgt_Thr', where, 'PYTHON3')

# Non-native species weighting - stress
arcpy.AddField_management(rs, 'NN_Wgt_Str', 'DOUBLE')
where = '(!InvSpp_Str!) * 2'
arcpy.CalculateField_management(rs, 'NN_Wgt_Str', where, 'PYTHON3')

# Non-native species weighting - threat
arcpy.AddField_management(rs, 'NN_Wgt_Thr', 'DOUBLE')
where = '(!InvSpp_Thr!) * 2'
arcpy.CalculateField_management(rs, 'NN_Wgt_Thr', where, 'PYTHON3')

# Ungulate weighting - stress (same as threat)
arcpy.AddField_management(rs, 'Ung_Wgt_Str', 'DOUBLE')
where = '(!UnguStr!) * 1'
arcpy.CalculateField_management(rs, 'Ung_Wgt_Str', where, 'PYTHON3')

# Ungulate weighting - threat (same as stress)
arcpy.AddField_management(rs, 'Ung_Wgt_Thr', 'DOUBLE')
where = '(!UnguStr!) * 1'
arcpy.CalculateField_management(rs, 'Ung_Wgt_Thr', where, 'PYTHON3')


# Sum all stress and threat scores
[f.name for f in arcpy.ListFields(rs)]
arcpy.AddField_management(rs, 'ScoreStr', 'DOUBLE')
arcpy.CalculateField_management(rs, 'ScoreStr', '!WW_Wgt_Str! + !Cl_Wgt_Str! + !HA_Wgt_Str! + !NN_Wgt_Str! + !Ung_Wgt_Str!', 'PYTHON3')
arcpy.AddField_management(rs, 'ScoreThr', 'DOUBLE')
arcpy.CalculateField_management(rs, 'ScoreThr', '!WW_Wgt_Thr! + !Cl_Wgt_Thr! + !HA_Wgt_Thr! + !NN_Wgt_Thr! + !Ung_Wgt_Thr!', 'PYTHON3')





#-------------------------------------------------------------------------------
# Lakes/playas
[f.name for f in arcpy.ListFields(lp)]
null2zero(lp, 'PumpStr')
null2zero(lp, 'CommitThr')
null2zero(lp, 'Shallow_GW')
null2zero(lp, 'Falling_GW')
null2zero(lp, 'Recharge')
null2zero(lp, 'ClmThrR')
null2zero(lp, 'InvSpp_Thr')
null2zero(lp, 'InvSpp_Str')
null2zero(lp, 'HousingThr')
null2zero(lp, 'HousingStr')

# Water withdrawal weighting - stress
max1 = getMax(lp, 'Falling_GW')
max2 = getMax(lp, 'PumpStr')
max_norm = max1 + max2
print(max_norm)
arcpy.AddField_management(lp, 'WW_Norm_Str', 'DOUBLE')
where = '(!Falling_GW! + !PumpStr!) / ' + "%s" %max_norm
arcpy.CalculateField_management(lp, 'WW_Norm_Str', where, 'PYTHON3')
arcpy.AddField_management(lp, 'WW_Wgt_Str', 'DOUBLE')
where = '(!WW_Norm_Str!) * 5'
arcpy.CalculateField_management(lp, 'WW_Wgt_Str', where, 'PYTHON3')

# Water withdrawal weighting - threat
max1 = getMax(lp, 'Shallow_GW')
max2 = getMax(lp, 'CommitThr')
max_norm = max1 + max2
print(max_norm)
arcpy.AddField_management(lp, 'WW_Norm_Thr', 'DOUBLE')
where = '(!Shallow_GW! + !CommitThr!) / ' + "%s" %max_norm
arcpy.CalculateField_management(lp, 'WW_Norm_Thr', where, 'PYTHON3') 
arcpy.AddField_management(lp, 'WW_Wgt_Thr', 'DOUBLE')
where = '(!WW_Norm_Thr!) * 5'
arcpy.CalculateField_management(lp, 'WW_Wgt_Thr', where, 'PYTHON3')


# Climate weighting - stress
arcpy.AddField_management(lp, 'Cl_Wgt_Str', 'DOUBLE')
where = '(!Recharge!) * 4'
arcpy.CalculateField_management(lp, 'Cl_Wgt_Str', where, 'PYTHON3')

# Climate weighting - threat
arcpy.AddField_management(lp, 'Cl_Wgt_Thr', 'DOUBLE')
where = '(!ClmThrR!) * 4'
arcpy.CalculateField_management(lp, 'Cl_Wgt_Thr', where, 'PYTHON3')

# Human activity weighting - stress
max1 = getMax(lp, 'HousingStr')
max2 = getMax(lp, 'SurfaceDiv')
max_norm = max1 + max2
print(max_norm)
arcpy.AddField_management(lp, 'HA_Norm_Str', 'DOUBLE')
where = '(!HousingStr! + !SurfaceDiv!) / ' + "%s" %max_norm
arcpy.CalculateField_management(lp, 'HA_Norm_Str', where, 'PYTHON3')
arcpy.AddField_management(lp, 'HA_Wgt_Str', 'DOUBLE')
where = '(!HA_Norm_Str!) * 3'
arcpy.CalculateField_management(lp, 'HA_Wgt_Str', where, 'PYTHON3')

# Human activity weighting - threat
arcpy.AddField_management(lp, 'HA_Wgt_Thr', 'DOUBLE')
where = '(!HousingThr!) * 3'
arcpy.CalculateField_management(lp, 'HA_Wgt_Thr', where, 'PYTHON3')

# Non-native species weighting - stress
arcpy.AddField_management(lp, 'NN_Wgt_Str', 'DOUBLE')
where = '(!InvSpp_Str!) * 2'
arcpy.CalculateField_management(lp, 'NN_Wgt_Str', where, 'PYTHON3')

# Non-native species weighting - threat
arcpy.AddField_management(lp, 'NN_Wgt_Thr', 'DOUBLE')
where = '(!InvSpp_Thr!) * 2'
arcpy.CalculateField_management(lp, 'NN_Wgt_Thr', where, 'PYTHON3')

# Ungulate weighting - stress (same as threat)
arcpy.AddField_management(lp, 'Ung_Wgt_Str', 'DOUBLE')
where = '(!UnguStr!) * 1'
arcpy.CalculateField_management(lp, 'Ung_Wgt_Str', where, 'PYTHON3')

# Ungulate weighting - threat (same as stress)
arcpy.AddField_management(lp, 'Ung_Wgt_Thr', 'DOUBLE')
where = '(!UnguStr!) * 1'
arcpy.CalculateField_management(lp, 'Ung_Wgt_Thr', where, 'PYTHON3')


# Sum all stress and threat scores
[f.name for f in arcpy.ListFields(lp)]
arcpy.AddField_management(lp, 'ScoreStr', 'DOUBLE')
arcpy.CalculateField_management(lp, 'ScoreStr', '!WW_Wgt_Str! + !Cl_Wgt_Str! + !HA_Wgt_Str! + !NN_Wgt_Str! + !Ung_Wgt_Str!', 'PYTHON3')
arcpy.AddField_management(lp, 'ScoreThr', 'DOUBLE')
arcpy.CalculateField_management(lp, 'ScoreThr', '!WW_Wgt_Thr! + !Cl_Wgt_Thr! + !HA_Wgt_Thr! + !NN_Wgt_Thr! + !Ung_Wgt_Thr!', 'PYTHON3')



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Bivariate symbology for some reason not possible for phreatophytes or streams/rivers
# Can 

# Phreatophytes
[f.name for f in arcpy.ListFields(phr2)]
arcpy.AddField_management(phr2, 'Bivariate', 'TEXT')
with arcpy.da.UpdateCursor(phr2, ['ScoreStr', 'ScoreThr', 'Bivariate']) as cursor:
    for row in cursor:
        if row[0] < 1 and row[1] < 1:
            row[2] = 'Both low'
            cursor.updateRow(row)
        elif (row[0] >= 1 and row[0] < 5) and (row[1] >= 1 and row[1] < 5):
            row[2] = 'Both med'
            cursor.updateRow(row)
        elif (row[0] >= 5) and (row[1] >= 5):
            row[2] = 'Both high'
            cursor.updateRow(row)
        elif (row[0] < 1) and (row[1] >= 1 and row[1] < 5):
            row[2] = 'Low stress, med threat'
            cursor.updateRow(row)
        elif (row[0] >= 1 and row[0] < 5) and (row[1] < 1):
            row[2] = 'Med stress, low threat'
            cursor.updateRow(row)            
        elif (row[0] < 1) and (row[1] >= 5):
            row[2] = 'Low stress, high threat'
            cursor.updateRow(row)   
        elif (row[0] >= 5) and (row[1] < 1):
            row[2] = 'High stress, low threat'
            cursor.updateRow(row)
        elif (row[0] >= 1 and row[0] < 5) and (row[1] >= 5):
            row[2] = 'Med stress, high threat'
            cursor.updateRow(row) 
        elif (row[0] >= 5) and (row[1] >= 1 and row[1] < 5):
            row[2] = 'High stress, med threat'
            cursor.updateRow(row)
        print(row[2])    
del cursor            
            
            
# Rivers/streams
[f.name for f in arcpy.ListFields(rs)]
arcpy.AddField_management(rs, 'Bivariate', 'TEXT')
with arcpy.da.UpdateCursor(rs, ['ScoreStr', 'ScoreThr', 'Bivariate']) as cursor:
    for row in cursor:
        if row[0] < 1 and row[1] < 1:
            row[2] = 'Both low'
            cursor.updateRow(row)
        elif (row[0] >= 1 and row[0] < 5) and (row[1] >= 1 and row[1] < 5):
            row[2] = 'Both med'
            cursor.updateRow(row)
        elif (row[0] >= 5) and (row[1] >= 5):
            row[2] = 'Both high'
            cursor.updateRow(row)
        elif (row[0] < 1) and (row[1] >= 1 and row[1] < 5):
            row[2] = 'Low stress, med threat'
            cursor.updateRow(row)
        elif (row[0] >= 1 and row[0] < 5) and (row[1] < 1):
            row[2] = 'Med stress, low threat'
            cursor.updateRow(row)            
        elif (row[0] < 1) and (row[1] >= 5):
            row[2] = 'Low stress, high threat'
            cursor.updateRow(row)   
        elif (row[0] >= 5) and (row[1] < 1):
            row[2] = 'High stress, low threat'
            cursor.updateRow(row)
        elif (row[0] >= 1 and row[0] < 5) and (row[1] >= 5):
            row[2] = 'Med stress, high threat'
            cursor.updateRow(row) 
        elif (row[0] >= 5) and (row[1] >= 1 and row[1] < 5):
            row[2] = 'High stress, med threat'
            cursor.updateRow(row)
        print(row[2])    
del cursor            






#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Summarize by basin
# Only for priority basins - which we have all data inputs for

basins = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Maps\GDE_Threats.gdb\hydrographic_basins_priority




