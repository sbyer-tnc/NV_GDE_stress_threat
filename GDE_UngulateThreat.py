#-------------------------------------------------------------------------------
# Name:        Threats to GDEs due to ungulates
# Purpose:     Determine whether GDEs are at risk of damage from ungulates
#              Depends on GDE type and ungulate potential presence
#              
# Refer to spreadsheet for notes on affected systems:  D:\GDE_Threats\OtherData\Grazing\grazing_impacts_by_GDE.xlsx
# GDE methods report has descriptions of systems for each layer: https://www.conservationgateway.org/ConservationByGeography/NorthAmerica/UnitedStates/nevada/water/Documents/NV_iGDE_MethodsReport.pdf
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


def table_to_data_frame(in_table, input_fields=None, where_clause=None):
    """Function will convert an arcgis table into a pandas dataframe with an object ID index, and the selected
    input fields using an arcpy.da.SearchCursor."""
    OIDFieldName = arcpy.Describe(in_table).OIDFieldName
    if input_fields:
        final_fields = [OIDFieldName] + input_fields
    else:
        final_fields = [field.name for field in arcpy.ListFields(in_table)]
    data = [row for row in arcpy.da.SearchCursor(in_table, final_fields, where_clause=where_clause)]
    fc_dataframe = pd.DataFrame(data, columns=final_fields)
    fc_dataframe = fc_dataframe.set_index(OIDFieldName, drop=True)
    return fc_dataframe

#-------------------------------------------------------------------------------

# Ungulate potential presence data

graze = arcpy.CopyFeatures_management(r'D:\GDE_Threats\OtherData\Grazing\BLM_National_Grazing_Allotments\gra.gdb\gra_allot_poly', 'graze')
hma = arcpy.CopyFeatures_management(r'D:\GDE_Threats\OtherData\BLM_National_Wild_Horse_and_Burro\blm_natl_whb_geocortex.gdb\whb_hma_pop_poly', 'hma')
elk = arcpy.CopyFeatures_management(r'D:\GDE_Threats\OtherData\Occupied_Elk_Distribution.shp', 'elk')

#-------------------------------------------------------------------------------

# GDE data

phr = r'D:\GDE_Threats\Hydrology\NV_iGDE_022120.gdb\Phreatophytes'
wet = r'D:\GDE_Threats\Hydrology\NV_iGDE_022120.gdb\Wetlands'
springs = r'D:\GDE_Threats\Hydrology\NV_iGDE_022120.gdb\Springs'
# No impacts to lakes, playas

# Streams/rivers?????

#-------------------------------------------------------------------------------

# Hydro areas and hexagons for summarizing
# Only copy once, then read in

ha = arcpy.CopyFeatures_management(r'D:\GDE_Threats\Hydrology\NV_iGDE_Story_022120_shp\NV_HydrographicAreas.shp', 'ha_ungulates')
hexagons = arcpy.CopyFeatures_management(r'D:\GDE_Threats\Hydrology\NV_iGDE_Story_022120_shp\NV_Hexagons.shp', 'hex_ungulates')

ha = 'ha_ungulates'
hexagons = 'hex_ungulates'

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Rules for groups of/individual GDE types

##########################################
# Phreatophytes
##########################################

# Aspen Woodland and Aspen Mixed-conifer
# Affected by:
# Cows = graze allot - cows
# Sheep = graze allot - sheep
# Elk = elk distribution

# Isolate aspen features
aspen = arcpy.Copy_management(phr, 'phr_aspen')
[f.name for f in arcpy.ListFields(aspen)]
with arcpy.da.UpdateCursor(aspen, ['PHR_TYPE']) as cursor:
    for row in cursor:
        if "Aspen" not in row[0]:
            print(row[0])
            cursor.deleteRow()
del cursor
arcpy.GetCount_management(aspen)

# Attribute aspen/forested features affected by ungulates - set values to 0
arcpy.AddField_management(aspen, 'GRAZE', 'LONG')
arcpy.AddField_management(aspen, 'HORSE_BURRO', 'LONG')
arcpy.AddField_management(aspen, 'ELK', 'LONG')
arcpy.CalculateField_management(aspen, 'GRAZE', 0, 'PYTHON3')
arcpy.CalculateField_management(aspen, 'HORSE_BURRO', 0, 'PYTHON3')
arcpy.CalculateField_management(aspen, 'ELK', 0, 'PYTHON3')

# Separate into individual features
aspen = arcpy.MultipartToSinglepart_management(aspen, 'aspen_explode')

# Impacts from grazing
aspen_select = arcpy.SelectLayerByLocation_management(aspen, 'INTERSECT', graze, '', 'NEW_SELECTION')
int(arcpy.GetCount_management(aspen_select)[0]) 
arcpy.CalculateField_management(aspen_select, 'GRAZE', 1, 'PYTHON3')

# Elk distribution
aspen_select = arcpy.SelectLayerByLocation_management(aspen, 'INTERSECT', elk, '', 'NEW_SELECTION')
int(arcpy.GetCount_management(aspen_select)[0]) 
arcpy.CalculateField_management(aspen_select, 'ELK', 1, 'PYTHON3')

# Sum up ungulate impacts
arcpy.AddField_management(aspen, 'UNGULATES', 'LONG')
arcpy.CalculateField_management(aspen, 'UNGULATES', '!GRAZE! + !HORSE_BURRO! + !ELK!', 'PYTHON3')

#------------------------------------------------------
# Mesquite
# Affected by:
# Horse/burros
# Cows - in the Mojave if grazing permits are not retired

# MORMON PEAK NV01044
# ROX-TULE NV01063
# GARDEN SPRING NV01065
# GOURD SPRING NV01071
# SNOW SPRINGS NV01074
# SUMMIT SPRING NV01077
# WHITE ROCK NV01078
# PAHRANAGAT WEST NV01081
# DELAMAR NV01083
# BREEDLOVE NV11010
# LOWER LAKE WEST NV11013
# GRAPEVINE NV11032
# HENRIE COMPLEX NV11034
# BALD MTN NV21003
# LIME MOUNTAIN NV21005
# LOWER LAKE EAST NV21022
# Hidden Valley NV15412
# CRESCENT (N-4) NV00128
# RAZORBACK NV00093
# MACGRUDER MTN NV00099

active = ['NV01044', 'NV01063', 'NV01065', 'NV01071', 'NV01074', 'NV01077',
'NV01078', 'NV01081', 'NV01083', 'NV11010', 'NV11013', 'NV11032',
'NV11034', 'NV21003', 'NV21005', 'NV21022', 'NV15412', 'NV00128',
'NV00093', 'NV00099']

# Isolate mesquite features
mesquite = arcpy.Copy_management(phr, 'phr_mesquite')
[f.name for f in arcpy.ListFields(mesquite)]
with arcpy.da.UpdateCursor(mesquite, ['PHR_TYPE']) as cursor:
    for row in cursor:
        if "Mesquite" not in row[0]:
            print(row[0])
            cursor.deleteRow()
del cursor
arcpy.GetCount_management(mesquite)

# Attribute mesquite features affected by ungulates - set values to 0
arcpy.AddField_management(mesquite, 'GRAZE', 'LONG')
arcpy.AddField_management(mesquite, 'HORSE_BURRO', 'LONG')
arcpy.AddField_management(mesquite, 'ELK', 'LONG')
arcpy.CalculateField_management(mesquite, 'GRAZE', 0, 'PYTHON3')
arcpy.CalculateField_management(mesquite, 'HORSE_BURRO', 0, 'PYTHON3')
arcpy.CalculateField_management(mesquite, 'ELK', 0, 'PYTHON3')

# Separate into individual features
mesquite = arcpy.MultipartToSinglepart_management(mesquite, 'mesquite_explode')

# Horse/burro distribution
mesquite_select = arcpy.SelectLayerByLocation_management(mesquite, 'INTERSECT', hma, '', 'NEW_SELECTION')
int(arcpy.GetCount_management(mesquite_select)[0]) 
arcpy.CalculateField_management(mesquite_select, 'HORSE_BURRO', 1, 'PYTHON3')


# Sum up ungulate impacts
arcpy.AddField_management(mesquite, 'UNGULATES', 'LONG')
arcpy.CalculateField_management(mesquite, 'UNGULATES', '!GRAZE! + !HORSE_BURRO! + !ELK!', 'PYTHON3')

#------------------------------------------------------
# Summarize By Basin

############################
# Aspen
aspen = 'aspen_explode'
aspen_ha = arcpy.Union_analysis([aspen, r'D:\GDE_Threats\Hydrology\NV_iGDE_Story_022120_shp\NV_HydrographicAreas.shp'], 'aspen_ha')
arcpy.GetCount_management(aspen_ha)

# Area of each meadow polygon
arcpy.AddField_management(aspen_ha, "ACRES", "DOUBLE")
arcpy.CalculateGeometryAttributes_management(aspen_ha, [["ACRES", "AREA"]], "", "ACRES", env.outputCoordinateSystem)
[f.name for f in arcpy.ListFields(aspen_ha)]

# Create binary attribute that will be used to ID polygons as GDE and whether they are ungulate impacted
arcpy.AddField_management(aspen_ha, 'GDE_UNG', 'LONG')
arcpy.CalculateField_management(aspen_ha, 'GDE_UNG', 0, 'PYTHON3') # Only non GDEs will have GDE_UNG = 0
with arcpy.da.UpdateCursor(aspen_ha, ['UNGULATES', 'PHR_TYPE', 'GDE_UNG']) as cursor:
    for row in cursor:
        if (row[1] is not None) and (row[0] > 0):
            row[2] = 10
        elif (row[1] != "") and (row[0] == 0):
            row[2] = 1
        cursor.updateRow(row)
del cursor
# Impacted GDEs (10 = impacted GDE)            
# Total GDEs (1 = just GDE, may or may not be impacted)
# Non-GDEs (0 = not GDE, obviously not impacted)

# Calculate proportion of GDEs in basin that are/are not impacted
# Dissolve by GDE_UNG AND HYD_AREA; sum acres that are/are not GDE
aspen_dis = arcpy.Dissolve_management(aspen_ha, 'aspen_ha_dis', ['HYD_AREA', 'GDE_UNG'], [['ACRES', 'SUM']])
aspen_dis = 'aspen_ha_dis'

# For every HYD_AREA, create attributes for total GDE area, impacted GDE area, and proportion GDEs that are impacted
# Create pandas dataframe from attribute table
aspen_pd = table_to_data_frame(in_table = aspen_dis)
aspen_pd.columns

# Long to wide format on pandas dataframe - each hyd area has one row with attributes for GDE-impacted, GDE-Total, Non-GDE
aspen_df = aspen_pd.pivot(index = 'HYD_AREA', columns = 'GDE_UNG', values = 'SUM_ACRES')
aspen_df.columns
adf = aspen_df.rename(columns = {0: "NoGDE", 1: "GDENoImp", 10: "GDEImpact"})
adf = adf.fillna(0)

# Calculate proportion of GDEs in hyd area that have impact potential
#adf['ImpAspen'] = adf.apply(lambda row: row.GDEImpact / (row.GDENoImp + row.GDEImpact), axis=1)
adf['VulAspen'] = adf.apply(lambda row: row.GDEImpact, axis=1) # Raw acres of vulnerable aspen
adf['AllAspen'] = adf.apply(lambda row: (row.GDENoImp + row.GDEImpact), axis=1)
adf.to_csv(r'D:\GDE_Threats\OtherData\Grazing\aspen_ha_impact_prop.csv')

# Join proportion values back to original hexes
adf = arcpy.TableToTable_conversion(r'D:\GDE_Threats\OtherData\Grazing\aspen_ha_impact_prop.csv', env.workspace, 'adf_temp')
arcpy.JoinField_management(ha, 'HYD_AREA', adf, 'HYD_AREA', ['AllAspen', 'VulAspen'])

############################
# Mesquite
mesquite = 'mesquite_explode'
mesquite_ha = arcpy.Union_analysis([mesquite, r'D:\GDE_Threats\Hydrology\NV_iGDE_Story_022120_shp\NV_HydrographicAreas.shp'], 'mesquite_ha')
arcpy.GetCount_management(mesquite_ha)

# Area of each mesquite polygon
arcpy.AddField_management(mesquite_ha, "ACRES", "DOUBLE")
arcpy.CalculateGeometryAttributes_management(mesquite_ha, [["ACRES", "AREA"]], "", "ACRES", env.outputCoordinateSystem)
[f.name for f in arcpy.ListFields(mesquite_ha)]

# Create binary attribute that will be used to ID polygons as GDE and whether they are ungulate impacted
arcpy.AddField_management(mesquite_ha, 'GDE_UNG', 'LONG')
arcpy.CalculateField_management(mesquite_ha, 'GDE_UNG', 0, 'PYTHON3') # Only non GDEs will have GDE_UNG = 0
with arcpy.da.UpdateCursor(mesquite_ha, ['UNGULATES', 'PHR_TYPE', 'GDE_UNG']) as cursor:
    for row in cursor:
        if (row[1] is not None) and (row[0] > 0):
            row[2] = 10
        elif (row[1] != "") and (row[0] == 0):
            row[2] = 1
        cursor.updateRow(row)
del cursor
# Impacted GDEs (10 = impacted GDE)            
# Total GDEs (1 = just GDE, may or may not be impacted)
# Non-GDEs (0 = not GDE, obviously not impacted)

# Calculate proportion of GDEs in basin that are/are not impacted
# Dissolve by GDE_UNG AND HYD_AREA; sum acres that are/are not GDE
mesquite_dis = arcpy.Dissolve_management(mesquite_ha, 'mesquite_ha_dis', ['HYD_AREA', 'GDE_UNG'], [['ACRES', 'SUM']])
mesquite_dis = 'mesquite_ha_dis'

# For every HYD_AREA, create attributes for total GDE area, impacted GDE area, and proportion GDEs that are impacted
# Create pandas dataframe from attribute table
mesquite_pd = table_to_data_frame(in_table = mesquite_dis)
mesquite_pd.columns

# Long to wide format on pandas dataframe - each hyd area has one row with attributes for GDE-impacted, GDE-Total, Non-GDE
mesquite_df = mesquite_pd.pivot(index = 'HYD_AREA', columns = 'GDE_UNG', values = 'SUM_ACRES')
mesquite_df.columns
mdf = mesquite_df.rename(columns = {0: "NoGDE", 1: "GDENoImp", 10: "GDEImpact"})
mdf = mdf.fillna(0)

# Calculate acres of GDEs in hyd area that have impact potential
#mdf['ImpMesquite'] = mdf.apply(lambda row: row.GDEImpact / (row.GDENoImp + row.GDEImpact), axis=1) # Proportion
mdf['VulMesquite'] = mdf.apply(lambda row: row.GDEImpact, axis=1) # Raw acres of vulnerable mesquite
mdf['AllMesquite'] = mdf.apply(lambda row: (row.GDENoImp + row.GDEImpact), axis=1)
mdf.to_csv(r'D:\GDE_Threats\OtherData\Grazing\mesquite_ha_impact_prop.csv')

# Join proportion values back to original HAs
mdf = arcpy.TableToTable_conversion(r'D:\GDE_Threats\OtherData\Grazing\mesquite_ha_impact_prop.csv', env.workspace, 'msdf_temp')
arcpy.JoinField_management(ha, 'HYD_AREA', mdf, 'HYD_AREA', ['VulMesquite', 'AllMesquite'])

############################
# Create field to hold sum of acres that are impactable
ha = 'ha_ungulates'
[f.name for f in arcpy.ListFields(ha)]

# Acres of vulnerable GDEs
arcpy.AddField_management(ha, 'Phr_Vul', 'DOUBLE')
arcpy.CalculateField_management(ha, 'Phr_Vul', '(!VulAspen! + !VulMesquite!)', 'PYTHON3')

# % of impactable phreatophytes that are vulnerable - divide by sum of aspen and mesquite acres! Not AREA_PHR
arcpy.AddField_management(ha, 'Phr_VulPer', 'DOUBLE')
arcpy.CalculateField_management(ha, 'Phr_VulPer', '(!VulAspen! + !VulMesquite!) / (!AllAspen! + !AllMesquite!)', 'PYTHON3')

#------------------------------------------------------
# Summarize By Hexagon

############################
# Aspen
aspen_hex = arcpy.Union_analysis([aspen, hexagons], 'aspen_hex')
arcpy.GetCount_management(aspen_hex)
[f.name for f in arcpy.ListFields(aspen_hex)]

# (Re)Calculate acres of chunked GDE polygons
arcpy.AddField_management(aspen_hex, "ACRES", "DOUBLE")
arcpy.CalculateGeometryAttributes_management(aspen_hex, [["ACRES", "AREA"]], "", "ACRES", env.outputCoordinateSystem)

# Create binary attribute that will be used to ID polygons as GDE and whether they are ungulate impacted
arcpy.AddField_management(aspen_hex, 'GDE_UNG', 'LONG')
arcpy.CalculateField_management(aspen_hex, 'GDE_UNG', 0, 'PYTHON3') # Only non GDEs will have GDE_UNG = 0
with arcpy.da.UpdateCursor(aspen_hex, ['UNGULATES', 'PHR_TYPE', 'GDE_UNG']) as cursor:
    for row in cursor:
        if (row[1] is not None) and (row[0] > 0):
            row[2] = 10
        elif (row[1] != "") and (row[0] == 0):
            row[2] = 1
        cursor.updateRow(row)
del cursor

# Calculate proportion of GDEs in hex that are/are not impacted
# Dissolve by GDE_UNG AND Hex_ID; sum acres that are/are not GDE
aspen_dis = arcpy.Dissolve_management(aspen_hex, 'aspen_hex_dis', ['Hex_ID', 'GDE_UNG'], [['ACRES', 'SUM']])
aspen_dis = 'aspen_hex_dis'

# For every Hex_ID, create attributes for total GDE area, impacted GDE area, and proportion GDEs that are impacted
# Create pandas dataframe from attribute table
aspen_pd = table_to_data_frame(in_table = aspen_dis)
aspen_pd.columns

# Long to wide format on pandas dataframe - each hyd area has one row with attributes for GDE-impacted, GDE-Total, Non-GDE
aspen_df = aspen_pd.pivot(index = 'Hex_ID', columns = 'GDE_UNG', values = 'SUM_ACRES')
aspen_df.columns
adf = aspen_df.rename(columns = {0: "NoGDE", 1: "GDENoImp", 10: "GDEImpact"})
adf = adf.fillna(0)

# Calculate proportion of GDEs in hyd area that have impact potential
adf['VulAspen'] = adf.apply(lambda row: row.GDEImpact, axis=1) # Raw acres of vulnerable aspen
adf['AllAspen'] = adf.apply(lambda row: (row.GDENoImp + row.GDEImpact), axis=1)
adf.to_csv(r'D:\GDE_Threats\OtherData\Grazing\aspen_hex_impact_prop.csv')

# Join proportion values back to original hexes
adf = arcpy.TableToTable_conversion(r'D:\GDE_Threats\OtherData\Grazing\aspen_hex_impact_prop.csv', env.workspace, 'adf_temp')
arcpy.JoinField_management(hexagons, 'Hex_ID', adf, 'Hex_ID', ['VulAspen', 'AllAspen'])

############################
# Mesquite
mesquite_hex = arcpy.Union_analysis([mesquite, hexagons], 'aspen_hex')
arcpy.GetCount_management(mesquite_hex)
[f.name for f in arcpy.ListFields(mesquite_hex)]

# (Re)Calculate acres of chunked GDE polygons
arcpy.AddField_management(mesquite_hex, "ACRES", "DOUBLE")
arcpy.CalculateGeometryAttributes_management(mesquite_hex, [["ACRES", "AREA"]], "", "ACRES", env.outputCoordinateSystem)

# Create binary attribute that will be used to ID polygons as GDE and whether they are ungulate impacted
arcpy.AddField_management(mesquite_hex, 'GDE_UNG', 'LONG')
arcpy.CalculateField_management(mesquite_hex, 'GDE_UNG', 0, 'PYTHON3') # Only non GDEs will have GDE_UNG = 0
with arcpy.da.UpdateCursor(mesquite_hex, ['UNGULATES', 'PHR_TYPE', 'GDE_UNG']) as cursor:
    for row in cursor:
        if (row[1] is not None) and (row[0] > 0):
            row[2] = 10
        elif (row[1] != "") and (row[0] == 0):
            row[2] = 1
        cursor.updateRow(row)
del cursor

# Calculate proportion of GDEs in hex that are/are not impacted
# Dissolve by GDE_UNG AND Hex_ID; sum acres that are/are not GDE
mesquite_dis = arcpy.Dissolve_management(mesquite_hex, 'mesquite_hex_dis', ['Hex_ID', 'GDE_UNG'], [['ACRES', 'SUM']])
mesquite_dis = 'mesquite_hex_dis'

# For every Hex_ID, create attributes for total GDE area, impacted GDE area, and proportion GDEs that are impacted
# Create pandas dataframe from attribute table
mesquite_pd = table_to_data_frame(in_table = mesquite_dis)
mesquite_pd.columns

# Long to wide format on pandas dataframe - each hyd area has one row with attributes for GDE-impacted, GDE-Total, Non-GDE
mesquite_df = mesquite_pd.pivot(index = 'Hex_ID', columns = 'GDE_UNG', values = 'SUM_ACRES')
mesquite_df.columns
mdf = mesquite_df.rename(columns = {0: "NoGDE", 1: "GDENoImp", 10: "GDEImpact"})
mdf = mdf.fillna(0)

# Calculate proportion of GDEs in hyd area that have impact potential
mdf['VulMesquite'] = mdf.apply(lambda row: row.GDEImpact, axis=1) # Raw acres of vulnerable mesquite
mdf['AllMesquite'] = mdf.apply(lambda row: (row.GDENoImp + row.GDEImpact), axis=1)
mdf.to_csv(r'D:\GDE_Threats\OtherData\Grazing\mesquite_hex_impact_prop.csv')

# Join proportion values back to original hexes
mdf = arcpy.TableToTable_conversion(r'D:\GDE_Threats\OtherData\Grazing\mesquite_hex_impact_prop.csv', env.workspace, 'msdf_temp')
arcpy.JoinField_management(hexagons, 'Hex_ID', mdf, 'Hex_ID', ['VulMesquite', 'AllMesquite'])

############################
# Create field to hold sum of acres that are impactable
hexagons = 'hex_ungulates'
[f.name for f in arcpy.ListFields(hexagons)]

# Acres of vulnerable GDEs
arcpy.AddField_management(hexagons, 'Phr_Vul', 'DOUBLE')
arcpy.CalculateField_management(hexagons, 'Phr_Vul', '(!VulAspen! + !VulMesquite!)', 'PYTHON3')

# % of impactable phreatophytes that are vulnerable - divide by sum of aspen and mesquite acres! Not AREA_PHR
arcpy.AddField_management(hexagons, 'Phr_VulPer', 'DOUBLE')
arcpy.CalculateField_management(hexagons, 'Phr_VulPer', '(!VulAspen! + !VulMesquite!) / (!AllAspen! + !AllMesquite!)', 'PYTHON3')



##########################################
# Wetlands
##########################################
# Wet meadows and montane riparian == meadow; palustrine forests == pf
# Affected by:
# Horse/burros
# Cows/sheep

# Currently we have ovevrlap between meadow and pf features
# Need to refine meadow selection - can't include 'palustrine forest'

# Refer to as variable 'meadow'

# Isolate meadow and montane riparian features
meadow = arcpy.Copy_management(wet, 'wet_meadow_riparian')
subtypes = ["meadow", "aquatic", "emergent", "montane"]
[f.name for f in arcpy.ListFields(meadow)]
with arcpy.da.UpdateCursor(meadow, ['WET_TYPE', 'WET_SUBTYPE']) as cursor:
    for row in cursor:
        if row[1] not in subtypes:
            print(row[1])
            cursor.deleteRow()
            if (row[0] != "Palustrine" or row[0] != "Riparian"):
                print(row[0])
                cursor.deleteRow()
            if (row[0] != "Palustrine" and row[0] != "Forest"):
                print(row[0])
                cursor.deleteRow()
del cursor
arcpy.GetCount_management(meadow)

# Attribute mesquite features affected by ungulates - set values to 0
arcpy.AddField_management(meadow, 'GRAZE', 'LONG')
arcpy.AddField_management(meadow, 'HORSE_BURRO', 'LONG')
arcpy.AddField_management(meadow, 'ELK', 'LONG')
arcpy.CalculateField_management(meadow, 'GRAZE', 0, 'PYTHON3')
arcpy.CalculateField_management(meadow, 'HORSE_BURRO', 0, 'PYTHON3')
arcpy.CalculateField_management(meadow, 'ELK', 0, 'PYTHON3')

# Separate into individual features
meadow = arcpy.MultipartToSinglepart_management(meadow, 'meadow_explode')

# Impacts from grazing
meadow_select = arcpy.SelectLayerByLocation_management(meadow, 'INTERSECT', graze, '', 'NEW_SELECTION')
int(arcpy.GetCount_management(meadow_select)[0]) 
arcpy.CalculateField_management(meadow_select, 'GRAZE', 1, 'PYTHON3')

# Horse/burro
meadow_select = arcpy.SelectLayerByLocation_management(meadow, 'INTERSECT', hma, '', 'NEW_SELECTION')
int(arcpy.GetCount_management(meadow_select)[0]) 
arcpy.CalculateField_management(meadow_select, 'HORSE_BURRO', 1, 'PYTHON3')

# Sum up ungulate impacts
arcpy.AddField_management(meadow, 'UNGULATES', 'LONG')
arcpy.CalculateField_management(meadow, 'UNGULATES', '!GRAZE! + !HORSE_BURRO! + !ELK!', 'PYTHON3')

#-----------------------------------------------------------------
# Palustrine/wet forest impacts
# Share same impact-by-ungulates as aspen systems
# No, move this to a wetland-impact section instead!!!!
# Add Palustrine-forest features from Wetlands layer as "aspen"
pf = arcpy.Copy_management(wet, 'wet_pf')
[f.name for f in arcpy.ListFields(pf)]
with arcpy.da.UpdateCursor(pf, ['WET_TYPE', 'WET_SUBTYPE']) as cursor:
    for row in cursor:
        if (row[0] != "Palustrine") and (row[1] != "Forest"):
            print(row[0])
            cursor.deleteRow()
del cursor
arcpy.GetCount_management(pf)

# Attribute aspen/forested features affected by ungulates - set values to 0
arcpy.AddField_management(pf, 'GRAZE', 'LONG')
arcpy.AddField_management(pf, 'HORSE_BURRO', 'LONG')
arcpy.AddField_management(pf, 'ELK', 'LONG')
arcpy.CalculateField_management(pf, 'GRAZE', 0, 'PYTHON3')
arcpy.CalculateField_management(pf, 'HORSE_BURRO', 0, 'PYTHON3')
arcpy.CalculateField_management(pf, 'ELK', 0, 'PYTHON3')

# Separate into individual features
pf = arcpy.MultipartToSinglepart_management(pf, 'pf_explode')

# Impacts from grazing
pf_select = arcpy.SelectLayerByLocation_management(pf, 'INTERSECT', graze, '', 'NEW_SELECTION')
int(arcpy.GetCount_management(pf_select)[0]) 
arcpy.CalculateField_management(pf_select, 'GRAZE', 1, 'PYTHON3')

# Elk distribution
pf_select = arcpy.SelectLayerByLocation_management(pf, 'INTERSECT', elk, '', 'NEW_SELECTION')
int(arcpy.GetCount_management(pf_select)[0]) 
arcpy.CalculateField_management(pf_select, 'ELK', 1, 'PYTHON3')

# Sum up ungulate impacts
arcpy.AddField_management(pf, 'UNGULATES', 'LONG')
arcpy.CalculateField_management(pf, 'UNGULATES', '!GRAZE! + !HORSE_BURRO! + !ELK!', 'PYTHON3')

#-----------------------------------------------------------------
# Summarize by basin

# Meadow
meadow = 'meadow_explode'
meadow_ha = arcpy.Union_analysis([meadow, r'D:\GDE_Threats\Hydrology\NV_iGDE_Story_022120_shp\NV_HydrographicAreas.shp'], 'meadow_ha')
arcpy.GetCount_management(meadow_ha)

# Area of each meadow polygon
arcpy.AddField_management(meadow_ha, "ACRES", "DOUBLE")
arcpy.CalculateGeometryAttributes_management(meadow_ha, [["ACRES", "AREA"]], "", "ACRES", env.outputCoordinateSystem)
[f.name for f in arcpy.ListFields(meadow_ha)]

# Create binary attribute that will be used to ID polygons as GDE and whether they are ungulate impacted
arcpy.AddField_management(meadow_ha, 'GDE_UNG', 'LONG')
arcpy.CalculateField_management(meadow_ha, 'GDE_UNG', 0, 'PYTHON3') # Only non GDEs will have GDE_UNG = 0
with arcpy.da.UpdateCursor(meadow_ha, ['UNGULATES', 'WET_TYPE', 'GDE_UNG']) as cursor:
    for row in cursor:
        if (row[1] is not None) and (row[0] > 0):
            row[2] = 10
        elif (row[1] != "") and (row[0] == 0):
            row[2] = 1
        cursor.updateRow(row)
del cursor
# Impacted GDEs (10 = impacted GDE)            
# Total GDEs (1 = just GDE, may or may not be impacted)
# Non-GDEs (0 = not GDE, obviously not impacted)

# Calculate proportion of GDEs in basin that are/are not impacted
# Dissolve by GDE_UNG AND HYD_AREA; sum acres that are/are not GDE
meadow_dis = arcpy.Dissolve_management(meadow_ha, 'meadow_ha_dis', ['HYD_AREA', 'GDE_UNG'], [['ACRES', 'SUM']])
meadow_dis = 'meadow_ha_dis'

# For every HYD_AREA, create attributes for total GDE area, impacted GDE area, and proportion GDEs that are impacted
# Create pandas dataframe from attribute table
meadow_pd = table_to_data_frame(in_table = meadow_dis)
meadow_pd.columns

# Long to wide format on pandas dataframe - each hyd area has one row with attributes for GDE-impacted, GDE-Total, Non-GDE
meadow_df = meadow_pd.pivot(index = 'HYD_AREA', columns = 'GDE_UNG', values = 'SUM_ACRES')
meadow_df.columns
mdf = meadow_df.rename(columns = {0: "NoGDE", 1: "GDENoImp", 10: "GDEImpact"})
mdf = mdf.fillna(0)

# Calculate proportion of GDEs in hyd area that have impact potential
mdf['VulMeadow'] = mdf.apply(lambda row: row.GDEImpact, axis=1) # Raw acres of vulnerable meadow
mdf['AllMeadow'] = mdf.apply(lambda row: (row.GDENoImp + row.GDEImpact), axis=1)
mdf.head(10)
mdf.to_csv(r'D:\GDE_Threats\OtherData\Grazing\meadow_ha_impact_prop.csv')

# Join proportion values and total meadow acres back to original HAs
mdf = arcpy.TableToTable_conversion(r'D:\GDE_Threats\OtherData\Grazing\meadow_ha_impact_prop.csv', env.workspace, 'mdf_temp')
arcpy.JoinField_management(ha, 'HYD_AREA', mdf, 'HYD_AREA', ['VulMeadow', 'AllMeadow'])



# Palustrine forest
pf = 'pf_explode'
pf_ha = arcpy.Union_analysis([pf, r'D:\GDE_Threats\Hydrology\NV_iGDE_Story_022120_shp\NV_HydrographicAreas.shp'], 'pf_ha')
arcpy.GetCount_management(pf_ha)

# Area of each meadow polygon
arcpy.AddField_management(pf_ha, "ACRES", "DOUBLE")
arcpy.CalculateGeometryAttributes_management(pf_ha, [["ACRES", "AREA"]], "", "ACRES", env.outputCoordinateSystem)
[f.name for f in arcpy.ListFields(pf_ha)]

# Create binary attribute that will be used to ID polygons as GDE and whether they are ungulate impacted
arcpy.AddField_management(pf_ha, 'GDE_UNG', 'LONG')
arcpy.CalculateField_management(pf_ha, 'GDE_UNG', 0, 'PYTHON3') # Only non GDEs will have GDE_UNG = 0
with arcpy.da.UpdateCursor(pf_ha, ['UNGULATES', 'WET_TYPE', 'GDE_UNG']) as cursor:
    for row in cursor:
        if (row[1] is not None) and (row[0] > 0):
            row[2] = 10
        elif (row[1] != "") and (row[0] == 0):
            row[2] = 1
        cursor.updateRow(row)
del cursor
# Impacted GDEs (10 = impacted GDE)            
# Total GDEs (1 = just GDE, may or may not be impacted)
# Non-GDEs (0 = not GDE, obviously not impacted)

# Calculate proportion of GDEs in basin that are/are not impacted
# Dissolve by GDE_UNG AND HYD_AREA; sum acres that are/are not GDE
pf_dis = arcpy.Dissolve_management(pf_ha, 'pf_ha_dis', ['HYD_AREA', 'GDE_UNG'], [['ACRES', 'SUM']])
pf_dis = 'pf_ha_dis'

# For every HYD_AREA, create attributes for total GDE area, impacted GDE area, and proportion GDEs that are impacted
# Create pandas dataframe from attribute table
pf_pd = table_to_data_frame(in_table = pf_dis)
pf_pd.columns

# Long to wide format on pandas dataframe - each hyd area has one row with attributes for GDE-impacted, GDE-Total, Non-GDE
pf_df = pf_pd.pivot(index = 'HYD_AREA', columns = 'GDE_UNG', values = 'SUM_ACRES')
pf_df.columns
pdf = pf_df.rename(columns = {0: "NoGDE", 1: "GDENoImp", 10: "GDEImpact"})
pdf = pdf.fillna(0)

# Calculate proportion of GDEs in hyd area that have impact potential
pdf['VulPalFrst'] = pdf.apply(lambda row: row.GDEImpact, axis=1) # Raw acres of vulnerable meadow
pdf['AllPalFrst'] = pdf.apply(lambda row: (row.GDENoImp + row.GDEImpact), axis=1)
pdf.head(10)
pdf.to_csv(r'D:\GDE_Threats\OtherData\Grazing\palforest_ha_impact_prop.csv')

# Join proportion values and total meadow acres back to original HAs
pdf = arcpy.TableToTable_conversion(r'D:\GDE_Threats\OtherData\Grazing\palforest_ha_impact_prop.csv', env.workspace, 'pdf_temp')
[f.name for f in arcpy.ListFields(pdf)]
arcpy.JoinField_management(ha, 'HYD_AREA', pdf, 'HYD_AREA', ['VulPalFrst', 'AllPalFrst'])



############################
# Create field to hold sum of acres that are impactable
ha = 'ha_ungulates'
[f.name for f in arcpy.ListFields(ha)]

# Acres of vulnerable GDEs
arcpy.AddField_management(ha, 'Wet_Vul', 'DOUBLE')
arcpy.CalculateField_management(ha, 'Wet_Vul', '(!VulMeadow! + !VulPalFrst!)', 'PYTHON3')

# % of impactable phreatophytes that are vulnerable - divide by sum of aspen and mesquite acres! Not AREA_PHR
arcpy.AddField_management(ha, 'Wet_VulPer', 'DOUBLE')
arcpy.CalculateField_management(ha, 'Wet_VulPer', '(!VulMeadow! + !VulPalFrst!) / (!AllMeadow! + !AllPalFrst!)', 'PYTHON3')


#-----------------------------------------------------------------
# Summarize by hexagon

# Meadows
meadow = 'meadow_explode'
meadow_hex = arcpy.Union_analysis([meadow, r'D:\GDE_Threats\Hydrology\NV_iGDE_Story_022120_shp\NV_HydrographicAreas.shp'], 'meadow_hex')
arcpy.GetCount_management(meadow_hex)

# Area of each meadow polygon
arcpy.AddField_management(meadow_hex, "ACRES", "DOUBLE")
arcpy.CalculateGeometryAttributes_management(meadow_hex, [["ACRES", "AREA"]], "", "ACRES", env.outputCoordinateSystem)
[f.name for f in arcpy.ListFields(meadow_hex)]

# Create binary attribute that will be used to ID polygons as GDE and whether they are ungulate impacted
arcpy.AddField_management(meadow_hex, 'GDE_UNG', 'LONG')
arcpy.CalculateField_management(meadow_hex, 'GDE_UNG', 0, 'PYTHON3') # Only non GDEs will have GDE_UNG = 0
with arcpy.da.UpdateCursor(meadow_hex, ['UNGULATES', 'WET_TYPE', 'GDE_UNG']) as cursor:
    for row in cursor:
        if (row[1] is not None) and (row[0] > 0):
            row[2] = 10
        elif (row[1] != "") and (row[0] == 0):
            row[2] = 1
        cursor.updateRow(row)
del cursor
# Impacted GDEs (10 = impacted GDE)            
# Total GDEs (1 = just GDE, may or may not be impacted)
# Non-GDEs (0 = not GDE, obviously not impacted)

# Calculate proportion of GDEs in basin that are/are not impacted
# Dissolve by GDE_UNG AND HYD_AREA; sum acres that are/are not GDE
meadow_dis = arcpy.Dissolve_management(meadow_hex, 'meadow_hex_dis', ['HYD_AREA', 'GDE_UNG'], [['ACRES', 'SUM']])
meadow_dis = 'meadow_hex_dis'

# For every HYD_AREA, create attributes for total GDE area, impacted GDE area, and proportion GDEs that are impacted
# Create pandas dataframe from attribute table
meadow_pd = table_to_data_frame(in_table = meadow_dis)
meadow_pd.columns

# Long to wide format on pandas dataframe - each hyd area has one row with attributes for GDE-impacted, GDE-Total, Non-GDE
meadow_df = meadow_pd.pivot(index = 'Hex_ID', columns = 'GDE_UNG', values = 'SUM_ACRES')
meadow_df.columns
mdf = meadow_df.rename(columns = {0: "NoGDE", 1: "GDENoImp", 10: "GDEImpact"})
mdf = mdf.fillna(0)

# Calculate proportion of GDEs in hyd area that have impact potential
mdf['VulMeadow'] = mdf.apply(lambda row: row.GDEImpact, axis=1) # Raw acres of vulnerable meadow
mdf['AllMeadow'] = mdf.apply(lambda row: (row.GDENoImp + row.GDEImpact), axis=1)
mdf.head(10)
mdf.to_csv(r'D:\GDE_Threats\OtherData\Grazing\meadow_hex_impact_prop.csv')

# Join proportion values and total meadow acres back to original HAs
mdf = arcpy.TableToTable_conversion(r'D:\GDE_Threats\OtherData\Grazing\meadow_hex_impact_prop.csv', env.workspace, 'mdf_temp')
arcpy.JoinField_management(hexagons, 'Hex_ID', mdf, 'Hex_ID', ['VulMeadow', 'AllMeadow'])



# Pal forests
pf_hex = arcpy.Union_analysis([pf, hexagons], 'pf_hex')
arcpy.GetCount_management(pf_hex)
[f.name for f in arcpy.ListFields(pf_hex)]

# (Re)Calculate acres of chunked GDE polygons
arcpy.AddField_management(pf_hex, "ACRES", "DOUBLE")
arcpy.CalculateGeometryAttributes_management(pf_hex, [["ACRES", "AREA"]], "", "ACRES", env.outputCoordinateSystem)

# Create binary attribute that will be used to ID polygons as GDE and whether they are ungulate impacted
arcpy.AddField_management(pf_hex, 'GDE_UNG', 'LONG')
arcpy.CalculateField_management(pf_hex, 'GDE_UNG', 0, 'PYTHON3') # Only non GDEs will have GDE_UNG = 0
with arcpy.da.UpdateCursor(pf_hex, ['UNGULATES', 'WET_TYPE', 'GDE_UNG']) as cursor:
    for row in cursor:
        if (row[1] is not None) and (row[0] > 0):
            row[2] = 10
        elif (row[1] != "") and (row[0] == 0):
            row[2] = 1
        cursor.updateRow(row)
del cursor

# Calculate proportion of GDEs in hex that are/are not impacted
# Dissolve by GDE_UNG AND Hex_ID; sum acres that are/are not GDE
pf_dis = arcpy.Dissolve_management(pf_hex, 'pf_hex_dis', ['Hex_ID', 'GDE_UNG'], [['ACRES', 'SUM']])
pf_dis = 'pf_hex_dis'

# For every Hex_ID, create attributes for total GDE area, impacted GDE area, and proportion GDEs that are impacted
# Create pandas dataframe from attribute table
pf_pd = table_to_data_frame(in_table = pf_dis)
pf_pd.columns

# Long to wide format on pandas dataframe - each hyd area has one row with attributes for GDE-impacted, GDE-Total, Non-GDE
pf_df = pf_pd.pivot(index = 'Hex_ID', columns = 'GDE_UNG', values = 'SUM_ACRES')
pf_df.columns
pdf = pf_df.rename(columns = {0: "NoGDE", 1: "GDENoImp", 10: "GDEImpact"})
pdf = pdf.fillna(0)

# Calculate proportion of GDEs in hyd area that have impact potential
pdf['VulPalFrst'] = pdf.apply(lambda row: row.GDEImpact, axis=1) # Raw acres of vulnerable meadow
pdf['AllPalFrst'] = pdf.apply(lambda row: (row.GDENoImp + row.GDEImpact), axis=1)
pdf.head(10)
pdf.to_csv(r'D:\GDE_Threats\OtherData\Grazing\palforest_hex_impact_prop.csv')

# Join proportion values back to original hexes
pdf = arcpy.TableToTable_conversion(r'D:\GDE_Threats\OtherData\Grazing\palforest_hex_impact_prop.csv', env.workspace, 'pdf_temp')
[f.name for f in arcpy.ListFields(pdf)]
arcpy.JoinField_management(hexagons, 'Hex_ID', pdf, 'Hex_ID', ['VulPalFrst', 'AllPalFrst'])

############################
# Create field to hold sum of acres that are impactable
hexagons = 'hex_ungulates'
[f.name for f in arcpy.ListFields(hexagons)]

# Acres of vulnerable GDEs
arcpy.AddField_management(hexagons, 'Wet_Vul', 'DOUBLE')
arcpy.CalculateField_management(hexagons, 'Wet_Vul', '(!VulMeadow! + !VulPalFrst!)', 'PYTHON3')

# % of impactable phreatophytes that are vulnerable - divide by sum of aspen and mesquite acres! Not AREA_PHR
arcpy.AddField_management(hexagons, 'Wet_VulPer', 'DOUBLE')
arcpy.CalculateField_management(hexagons, 'Wet_VulPer', '(!VulMeadow! + !VulPalFrst!) / (!AllMeadow! + !AllPalFrst!)', 'PYTHON3')



##########################################
# Springs
# Affected by:
# Horse/burros
# Cows/sheep

# Not considering whether springs are fenced - data not available on statewide scales

spr = arcpy.Copy_management(springs, 'spr_ungulates')
[f.name for f in arcpy.ListFields(spr)]

# Attribute mesquite features affected by ungulates - set values to 0
arcpy.AddField_management(spr, 'GRAZE', 'LONG')
arcpy.AddField_management(spr, 'HORSE_BURRO', 'LONG')
arcpy.AddField_management(spr, 'ELK', 'LONG')
arcpy.CalculateField_management(spr, 'GRAZE', 0, 'PYTHON3')
arcpy.CalculateField_management(spr, 'HORSE_BURRO', 0, 'PYTHON3')
arcpy.CalculateField_management(spr, 'ELK', 0, 'PYTHON3')

# Impacts from grazing
spr_select = arcpy.SelectLayerByLocation_management(spr, 'INTERSECT', graze, '', 'NEW_SELECTION')
int(arcpy.GetCount_management(spr_select)[0]) 
arcpy.CalculateField_management(spr_select, 'GRAZE', 1, 'PYTHON3')

# Horse/burro
spr_select = arcpy.SelectLayerByLocation_management(spr, 'INTERSECT', hma, '', 'NEW_SELECTION')
int(arcpy.GetCount_management(spr_select)[0]) 
arcpy.CalculateField_management(spr_select, 'HORSE_BURRO', 1, 'PYTHON3')

# Sum up ungulate impacts
arcpy.AddField_management(spr, 'UNGULATES', 'LONG')
arcpy.CalculateField_management(spr, 'UNGULATES', '!GRAZE! + !HORSE_BURRO! + !ELK!', 'PYTHON3')


#################
# Basin and hexagon-wide stats - % of at-risk-from-ungulates GDEs

# By Basin
springs = 'spr_ungulates'
springs_ha = arcpy.SpatialJoin_analysis(ha, springs, "springs_ha", "JOIN_ONE_TO_MANY", "KEEP_ALL", "", "INTERSECT")
arcpy.GetCount_management(springs_ha)
[f.name for f in arcpy.ListFields(springs)]

# Count of springs in HA
spring_ha_count = arcpy.Statistics_analysis(springs_ha, "springs_ha_count", [["UNGULATES", "COUNT"]], ["HYD_AREA", "UNGULATES"])
[f.name for f in arcpy.ListFields(spring_ha_count)]

# For every HYD_AREA, create attributes for total GDE area, impacted GDE area, and proportion GDEs that are impacted
# Create pandas dataframe from attribute table
spring_pd = table_to_data_frame(in_table = spring_ha_count)
spring_pd.columns

# Long to wide format on pandas dataframe - each hyd area has one row with attributes for GDE-impacted, GDE-Total, Non-GDE
spring_df = spring_pd.pivot_table(index = ['HYD_AREA'], columns = 'UNGULATES', values = 'COUNT_UNGULATES')
spring_df.columns
mdf = spring_df.rename(columns = {0: "NoImpact", 1: "Impact1", 2: "Impact2"})
mdf = mdf.fillna(0)

# Calculate proportion of GDEs in hyd area that have impact potential
mdf['ImpSpring'] = mdf.apply(lambda row: (row.Impact1 + row.Impact2) / (row.NoImpact + row.Impact1 + row.Impact2), axis=1)
mdf.to_csv(r'D:\GDE_Threats\OtherData\Grazing\springs_ha_impact_prop.csv')

# Join proportion values back to original hexes
spdf = arcpy.TableToTable_conversion(r'D:\GDE_Threats\OtherData\Grazing\springs_ha_impact_prop.csv', env.workspace, 'spdf_temp')
arcpy.JoinField_management(ha, 'HYD_AREA', spdf, 'HYD_AREA', ['ImpSpring'])



# By Hexagon
springs = 'spr_ungulates'
springs_hex = arcpy.SpatialJoin_analysis(hexagons, springs, "springs_hex", "JOIN_ONE_TO_MANY", "KEEP_ALL", "", "INTERSECT")
arcpy.GetCount_management(springs_hex)
[f.name for f in arcpy.ListFields(springs_hex)]

# Count of springs in Hexagons
spring_hex_count = arcpy.Statistics_analysis(springs_hex, "springs_hex_count", [["UNGULATES", "COUNT"]], ["Hex_ID", "UNGULATES"])
[f.name for f in arcpy.ListFields(spring_hex_count)]

# For every Hex_ID, create attributes for total GDE area, impacted GDE area, and proportion GDEs that are impacted
# Create pandas dataframe from attribute table
spring_pd = table_to_data_frame(in_table = spring_hex_count)
spring_pd.columns

# Long to wide format on pandas dataframe - each hexagon has one row with attributes for GDE-impacted, GDE-Total, Non-GDE
spring_df = spring_pd.pivot_table(index = 'Hex_ID', columns = 'UNGULATES', values = 'COUNT_UNGULATES')
spring_df.columns
mdf = spring_df.rename(columns = {0: "NoImpact", 1: "Impact1", 2: "Impact2"})
mdf = mdf.fillna(0)

# Calculate proportion of GDEs in hyd area that have impact potential
mdf['ImpSpring'] = mdf.apply(lambda row: (row.Impact1 + row.Impact2) / (row.NoImpact + row.Impact1 + row.Impact2), axis=1)
mdf.to_csv(r'D:\GDE_Threats\OtherData\Grazing\springs_hex_impact_prop.csv')

# Join proportion values back to original hexes
spdf = arcpy.TableToTable_conversion(r'D:\GDE_Threats\OtherData\Grazing\springs_hex_impact_prop.csv', env.workspace, 'spdf_temp')
arcpy.JoinField_management(hexagons, 'Hex_ID', spdf, 'Hex_ID', ['ImpSpring'])


[f.name for f in arcpy.ListFields(ha)]
[f.name for f in arcpy.ListFields(hexagons)]





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Start here
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# GDE layers
springs = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Springs'
wet = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Wetlands'
lp = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Lakes_Playas'
rs = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Rivers_Streams'
phr2 = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Phreatophytes_Explode'

# Ungulate layers
graze = arcpy.CopyFeatures_management(r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\Grazing\BLM_National_Grazing_Allotments\gra.gdb\gra_allot_poly', 'graze')
hma = arcpy.CopyFeatures_management(r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\BLM_National_Wild_Horse_and_Burro\blm_natl_whb_geocortex.gdb\whb_hma_pop_poly', 'hma')
elk = arcpy.CopyFeatures_management(r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\Occupied_Elk_Distribution.shp', 'elk')

# Ungulate impacts by GDE type:
# E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\Grazing\grazing_impacts_by_GDE.xlsx

#-------------------------------------------------------------------------------
# Attribute phreatophytes - separate by type

# Aspen
# Affected by:
# Cows = graze allot - cows
# Sheep = graze allot - sheep
# Elk = elk distribution

arcpy.AddField_management(phr2, 'GRAZE', 'LONG')
arcpy.AddField_management(phr2, 'HORSE_BURRO', 'LONG')
arcpy.AddField_management(phr2, 'ELK', 'LONG')
arcpy.CalculateField_management(phr2, 'GRAZE', 0, 'PYTHON3')
arcpy.CalculateField_management(phr2, 'HORSE_BURRO', 0, 'PYTHON3')
arcpy.CalculateField_management(phr2, 'ELK', 0, 'PYTHON3')


# Impacts from sheep and cow grazing
where = "PHR_TYPE LIKE '%Aspen%'"
aspen = arcpy.SelectLayerByAttribute_management(phr2, 'NEW_SELECTION', where) # Select aspen
int(arcpy.GetCount_management(aspen)[0]) 
aspen_select = arcpy.SelectLayerByLocation_management(aspen, 'INTERSECT', graze, '', 'SUBSET_SELECTION') # sub-select grazing
int(arcpy.GetCount_management(aspen_select)[0]) 
arcpy.CalculateField_management(aspen_select, 'GRAZE', 1, 'PYTHON3')

# Elk distribution
aspen = arcpy.SelectLayerByAttribute_management(phr2, 'NEW_SELECTION', where) # Select aspen
aspen_select = arcpy.SelectLayerByLocation_management(aspen, 'INTERSECT', elk, '', 'NEW_SELECTION') # sub-select elk
int(arcpy.GetCount_management(aspen_select)[0]) 
arcpy.CalculateField_management(aspen_select, 'ELK', 1, 'PYTHON3')

# Sum up ungulate impacts
arcpy.AddField_management(phr2, 'UNGULATES', 'LONG')
arcpy.CalculateField_management(phr2, 'UNGULATES', '!GRAZE! + !HORSE_BURRO! + !ELK!', 'PYTHON3')

# Mesquite
# Affected by:
# Horse/burros
# Cows - in the Mojave if grazing permits are not retired and during summer...?

# Horse/burro distribution
where = "PHR_TYPE LIKE '%Mesquite%'"
mesquite = arcpy.SelectLayerByAttribute_management(phr2, 'NEW_SELECTION', where) # Select mesquite
int(arcpy.GetCount_management(mesquite)[0]) 
mesquite_select = arcpy.SelectLayerByLocation_management(mesquite, 'INTERSECT', hma, '', 'SUBSET_SELECTION')
int(arcpy.GetCount_management(mesquite_select)[0]) 
arcpy.CalculateField_management(mesquite_select, 'HORSE_BURRO', 1, 'PYTHON3')

# Cow grazing in the Mojave - need to subset grazing allotment data
# Subset non-retired grazing permits (mesquite only in the mojave anyway)
# Check authorization use by allotment forms for any active allotment that overlaps with mapped mesquite
mesquite = arcpy.SelectLayerByAttribute_management(phr2, 'NEW_SELECTION', where) # Select mesquite
graze_select = arcpy.SelectLayerByLocation_management(graze, 'INTERSECT', mesquite, '', 'NEW_SELECTION')
int(arcpy.GetCount_management(graze_select)[0]) 
allots = arcpy.CopyFeatures_management(graze_select, 'mesquite_graze_allots')
allots = 'mesquite_graze_allots'
[f.name for f in arcpy.ListFields(allots)]

active = ['NV01044', 'NV01063', 'NV01065', 'NV01071', 'NV01074', 'NV01077',
'NV01078', 'NV01081', 'NV01083', 'NV11010', 'NV11013', 'NV11032',
'NV11034', 'NV21003', 'NV21005', 'NV21022', 'NV15412', 'NV00128',
'NV00093', 'NV00099']
with arcpy.da.UpdateCursor(allots, ['ST_ALLOT']) as cursor:
    for row in cursor:
        if row[0] in active:
            print(row[0])
        else:
            cursor.deleteRow()
del cursor

where = "PHR_TYPE LIKE '%Mesquite%'"
mesquite = arcpy.SelectLayerByAttribute_management(phr2, 'NEW_SELECTION', where) # Select mesquite
int(arcpy.GetCount_management(mesquite)[0]) 
mesquite_select = arcpy.SelectLayerByLocation_management(mesquite, 'INTERSECT', allots, '', 'SUBSET_SELECTION')
int(arcpy.GetCount_management(mesquite_select)[0]) 
arcpy.CalculateField_management(mesquite_select, 'GRAZE', 1, 'PYTHON3')


# Sum up ungulate impacts in PHREATOPHYTES
#arcpy.AddField_management(phr2, 'UNGULATES', 'LONG')
arcpy.CalculateField_management(phr2, 'UNGULATES', '!GRAZE! + !HORSE_BURRO! + !ELK!', 'PYTHON3')



#-------------------------------------------------------------------------------
# Attribute wetlands - separate by type

arcpy.AddField_management(wet, 'GRAZE', 'LONG')
arcpy.AddField_management(wet, 'HORSE_BURRO', 'LONG')
arcpy.AddField_management(wet, 'ELK', 'LONG')
arcpy.CalculateField_management(wet, 'GRAZE', 0, 'PYTHON3')
arcpy.CalculateField_management(wet, 'HORSE_BURRO', 0, 'PYTHON3')
arcpy.CalculateField_management(wet, 'ELK', 0, 'PYTHON3')
[f.name for f in arcpy.ListFields(wet)]

# Palustrine - Aquatic, Emergent, Meadow
# Affected by cow, sheep, and horse/burrow
# Impacts from sheep and cow grazing
where = "WET_TYPE = 'Palustrine' AND WET_SUBTYPE IN ('aquatic', 'emergent', 'meadow')"
aqua = arcpy.SelectLayerByAttribute_management(wet, 'NEW_SELECTION', where) # Select aquatic, emergent, meadow
int(arcpy.GetCount_management(aqua)[0]) 
aqua_select = arcpy.SelectLayerByLocation_management(aqua, 'INTERSECT', graze, '', 'SUBSET_SELECTION') # sub-select grazing
int(arcpy.GetCount_management(aqua_select)[0]) 
arcpy.CalculateField_management(aqua_select, 'GRAZE', 1, 'PYTHON3')

# Impacts from horse/burrow
aqua = arcpy.SelectLayerByAttribute_management(wet, 'NEW_SELECTION', where) # Select aquatic, emergent, meadow
int(arcpy.GetCount_management(aqua)[0]) 
aqua_select = arcpy.SelectLayerByLocation_management(aqua, 'INTERSECT', hma, '', 'SUBSET_SELECTION') # sub-select HMAs
int(arcpy.GetCount_management(aqua_select)[0]) 
arcpy.CalculateField_management(aqua_select, 'HORSE_BURRO', 1, 'PYTHON3')


# Riparian - montane
# Cows and sheep
where = "WET_TYPE = 'Riparian' AND WET_SUBTYPE = 'montane'"
aqua = arcpy.SelectLayerByAttribute_management(wet, 'NEW_SELECTION', where) # Select montane riparian
int(arcpy.GetCount_management(aqua)[0]) 
aqua_select = arcpy.SelectLayerByLocation_management(aqua, 'INTERSECT', graze, '', 'SUBSET_SELECTION') # sub-select grazing
int(arcpy.GetCount_management(aqua_select)[0]) 
arcpy.CalculateField_management(aqua_select, 'GRAZE', 1, 'PYTHON3')

# Horse and burro
aqua = arcpy.SelectLayerByAttribute_management(wet, 'NEW_SELECTION', where) # Select montane riparian
aqua_select = arcpy.SelectLayerByLocation_management(aqua, 'INTERSECT', hma, '', 'SUBSET_SELECTION') # sub-select grazing
int(arcpy.GetCount_management(aqua_select)[0]) 
arcpy.CalculateField_management(aqua_select, 'HORSE_BURRO', 1, 'PYTHON3')


# Sum up ungulate impacts in WETLANDS
arcpy.AddField_management(wet, 'UNGULATES', 'LONG')
arcpy.CalculateField_management(wet, 'UNGULATES', '!GRAZE! + !HORSE_BURRO! + !ELK!', 'PYTHON3')


#-------------------------------------------------------------------------------
# Attribute springs

arcpy.AddField_management(springs, 'GRAZE', 'LONG')
arcpy.AddField_management(springs, 'HORSE_BURRO', 'LONG')
arcpy.AddField_management(springs, 'ELK', 'LONG')
arcpy.CalculateField_management(springs, 'GRAZE', 0, 'PYTHON3')
arcpy.CalculateField_management(springs, 'HORSE_BURRO', 0, 'PYTHON3')
arcpy.CalculateField_management(springs, 'ELK', 0, 'PYTHON3')
[f.name for f in arcpy.ListFields(springs)]

# Attribute springs affected by horses/burrows and/or sheep/cows
# Sheep/cows
spr_select = arcpy.SelectLayerByLocation_management(springs, 'INTERSECT', graze, '', 'NEW_SELECTION') # sub-select grazing
int(arcpy.GetCount_management(spr_select)[0]) 
arcpy.CalculateField_management(spr_select, 'GRAZE', 1, 'PYTHON3')

# Horses/burros
spr_select = arcpy.SelectLayerByLocation_management(springs, 'INTERSECT', hma, '', 'NEW_SELECTION') # sub-select HMAs
int(arcpy.GetCount_management(spr_select)[0]) 
arcpy.CalculateField_management(spr_select, 'HORSE_BURRO', 1, 'PYTHON3')

# Sum up ungulate impacts in SPRINGS
arcpy.AddField_management(springs, 'UNGULATES', 'LONG')
arcpy.CalculateField_management(springs, 'UNGULATES', '!GRAZE! + !HORSE_BURRO! + !ELK!', 'PYTHON3')


#-------------------------------------------------------------------------------
# Attribute rivers/streams

arcpy.AddField_management(rs, 'GRAZE', 'LONG')
arcpy.AddField_management(rs, 'HORSE_BURRO', 'LONG')
arcpy.AddField_management(rs, 'ELK', 'LONG')
arcpy.CalculateField_management(rs, 'GRAZE', 0, 'PYTHON3')
arcpy.CalculateField_management(rs, 'HORSE_BURRO', 0, 'PYTHON3')
arcpy.CalculateField_management(rs, 'ELK', 0, 'PYTHON3')
[f.name for f in arcpy.ListFields(rs)]

# Attribute rivers/streams affected by horses/burrows and/or sheep/cows
# Sheep/cows
rs_select = arcpy.SelectLayerByLocation_management(rs, 'INTERSECT', graze, '', 'NEW_SELECTION') # sub-select grazing
int(arcpy.GetCount_management(rs_select)[0]) 
arcpy.CalculateField_management(rs_select, 'GRAZE', 1, 'PYTHON3')

# Horses/burros
rs_select = arcpy.SelectLayerByLocation_management(rs, 'INTERSECT', hma, '', 'NEW_SELECTION') # sub-select HMAs
int(arcpy.GetCount_management(rs_select)[0]) 
arcpy.CalculateField_management(rs_select, 'HORSE_BURRO', 1, 'PYTHON3')

# Sum up ungulate impacts in SPRINGS
arcpy.AddField_management(rs, 'UNGULATES', 'LONG')
arcpy.CalculateField_management(rs, 'UNGULATES', '!GRAZE! + !HORSE_BURRO! + !ELK!', 'PYTHON3')

#-------------------------------------------------------------------------------
# Attribute lakes/playas
# No impacts

arcpy.AddField_management(lp, 'UNGULATES', 'LONG')
arcpy.CalculateField_management(lp, 'UNGULATES', 0, 'PYTHON3')


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Normalize Ungulate score (0 to 1)
# Ungulate score is a STRESSOR - it's current
# Divide UNGULATES value by 2 - the maximum score of ungulate impacts for any feature

# GDE layers
springs = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Springs'
wet = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Wetlands'
lp = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Lakes_Playas'
rs = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Rivers_Streams'
phr2 = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Phreatophytes_Explode'

arcpy.AddField_management(springs, 'UnguStr', 'DOUBLE')
arcpy.CalculateField_management(springs, 'UnguStr', '!UNGULATES!/2', 'PYTHON3')

arcpy.AddField_management(wet, 'UnguStr', 'DOUBLE')
arcpy.CalculateField_management(wet, 'UnguStr', '!UNGULATES!/2', 'PYTHON3')

arcpy.AddField_management(phr2, 'UnguStr', 'DOUBLE')
arcpy.CalculateField_management(phr2, 'UnguStr', '!UNGULATES!/2', 'PYTHON3')

arcpy.AddField_management(rs, 'UnguStr', 'DOUBLE')
arcpy.CalculateField_management(rs, 'UnguStr', '!UNGULATES!/2', 'PYTHON3')

arcpy.AddField_management(lp, 'UnguStr', 'DOUBLE')
arcpy.CalculateField_management(lp, 'UnguStr', '!UNGULATES!/2', 'PYTHON3')


##########################################
# Summarize ungulate proportional impact by: hydrographic area, hexagon

# not sure how to consolidate the springs with the acres of impacted GDEs...
# - 40 acres of 100 total GDE acres in HA are impacted = 40%
# - 5 springs of 8 springs in HA are impacted = 63%

# Alternative idea - one GDE-Ungulate score for each feature:
# All GDE types start with value == 1 (max = 5: phr, spr, wet, lakeplaya, stream)
# Lake/playa and rivers/streams are not impacted, will always == 1 - or could remove these altogether...
# Percent of impact to each type will be subtracted from 1
# Output value = ungulate-impact score

# For each feature we have: 
# ImpAspen  == 1
# ImpMesquite == 1
# ImpMeadow (wetland) == 1
# ImpSpring == 1

# Sinze we only care about impactable GDEs - 
# May also be useful to know propertion imactable/unimpactable GDEs per feature


##########################################
## How to turn into basin- or hexagon-wide stats?
## % of at-risk-from-ungulates GDEs
#
## By Basin
## Ex. Upper Reese River Valley meadows
#meadow056 = arcpy.Clip_analysis(meadow, r'D:\GDE_Threats\Maps\GDE_Threats.gdb\UpperReeseRiver056', 'meadow056')
#arcpy.GetCount_management(meadow056)
#
## Area of each meadow polygon
#arcpy.AddField_management(meadow056, "ACRES", "DOUBLE")
#arcpy.CalculateGeometryAttributes_management(meadow056, [["ACRES", "AREA"]], "", "ACRES", env.outputCoordinateSystem)
#
## Percent meadows in  with UNGULATES > 0
## Summary statistics with SUM(ACRES) by UNGULATES as case field
#meadow056_stats = arcpy.Statistics_analysis(meadow056, 'meadow056_statistics', [['ACRES', 'SUM']], 'UNGULATES')
#meadow056_stats = arcpy.TableToTable_conversion(meadow056_stats, r'D:\\GDE_Threats\\OtherData\\Grazing', 'meadow056_stats.csv')
#
## Create pandas dataframe
#tbl = pd.read_csv(r'D:\\GDE_Threats\\OtherData\\Grazing\\meadow056_stats.csv', header=0, index_col=False)
#
## Get sum of all acres
## These depend on dataframe indices - may not work coorectly all the time...
#acres = tbl.sum(axis = 0)[3]
#
## Get acres where ungulates = 0
#acres0 = tbl.iloc[0][3]
#
## Get percent of where ungulates != 0
#acres_impact = acres - acres0
#percent_impact = acres_impact/acres
#print(percent_impact)
#
## 95% of Upper Reese River Valley's meadowed areas are affected by ungulates in some form
#
## By hexagon
## Ex. Upper Reese River Valley meadows by hexgaons
#meadow056 = arcpy.Union_analysis([meadow, r'D:\GDE_Threats\Maps\GDE_Threats.gdb\UpperReeseRiver056_hex'], 'meadow056_hex')
#arcpy.GetCount_management(meadow056)
#
## (Re)Calculate acres of chunked GDE polygons
##arcpy.AddField_management(meadow056, "ACRES", "DOUBLE")
#arcpy.CalculateGeometryAttributes_management(meadow056, [["ACRES", "AREA"]], "", "ACRES", env.outputCoordinateSystem)
#
## Create binary attribute that will be used to ID polygons as GDE and whether they are ungulate impacted
#arcpy.AddField_management(meadow056, 'GDE_UNG', 'LONG')
#arcpy.CalculateField_management(meadow056, 'GDE_UNG', 0, 'PYTHON3') # Only non GDEs will have GDE_UNG = 0
#with arcpy.da.UpdateCursor(meadow056, ['UNGULATES', 'WET_TYPE', 'GDE_UNG']) as cursor:
#    for row in cursor:
#        if (row[1] is not None) and (row[0] > 0):
#            row[2] = 10
#        elif (row[1] != "") and (row[0] == 0):
#            row[2] = 1
#        cursor.updateRow(row)
#del cursor
## Impacted GDEs (10 = impacted GDE)            
## Total GDEs (1 = just GDE, may or may not be impacted)
## Non-GDEs (0 = not GDE, obviously not impacted)
#
#
## Calculate proportion of GDEs in hexagon that are/are not impacted
## Dissolve by GDE_UNG AND Hexagon ID; sum acres that are/are not GDE
#meadow_dis = arcpy.Dissolve_management(meadow056, 'meadow056_hex_dis', ['Hex_ID', 'GDE_UNG'], [['ACRES', 'SUM']])
#
## For every Hex_ID, create attributes for total GDE area, impacted GDE area, and proportion GDEs that are impacted
## Create pandas dataframe from attribute table
#meadow_pd = table_to_data_frame(in_table = meadow_dis)
#
## Long to wide format on pandas dataframe - each hex has one row with attributes for GDE-impacted, GDE-Total, Non-GDE
#meadow_df = meadow_pd.pivot(index = 'Hex_ID', columns = 'GDE_UNG', values = 'SUM_ACRES')
#meadow_df.columns
#mdf = meadow_df.rename(columns = {0: "NoGDE", 1: "GDENoImp", 10: "GDEImpact"})
#mdf = mdf.fillna(0)
#
## Calculate proportion of GDEs in hexagon that have impact potential
#mdf['PropImpact'] = mdf.apply(lambda row: row.GDEImpact / (row.GDENoImp + row.GDEImpact), axis=1)
#mdf = mdf.fillna(0)
#mdf.to_csv(r'D:\GDE_Threats\OtherData\Grazing\meadow056_gde_impact_prop.csv')
#
## Join proportion values back to original hexes
#mdf = arcpy.TableToTable_conversion(r'D:\GDE_Threats\OtherData\Grazing\meadow056_gde_impact_prop.csv', env.workspace, 'mdf_temp')
#arcpy.JoinField_management(r'D:\GDE_Threats\Maps\GDE_Threats.gdb\UpperReeseRiver056_hex', 'Hex_ID', mdf, 'Hex_ID', ['PropImpact'])

# Works!


##########################################
# Previous code, outdated:

## Percent meadows in  with UNGULATES > 0
## Summary statistics with SUM(ACRES) by UNGULATES and HYD_AREA as case field
#aspen_stats = arcpy.Statistics_analysis(aspen_ha, 'aspen_ha_stats', [['ACRES', 'SUM']], ['UNGULATES', 'HYD_AREA'])
#aspen_stats = arcpy.TableToTable_conversion(aspen_stats, r'D:\\GDE_Threats\\OtherData\\Grazing', 'aspen_ha_stats.csv')
#
## Create pandas dataframe
#tbl = pd.read_csv(r'D:\\GDE_Threats\\OtherData\\Grazing\\aspen_ha_stats.csv', header=0, index_col=False)
#
## Long to wide format on pandas dataframe - each hex has one row with attributes for GDE-impacted, GDE-Total, Non-GDE
#aspen_df = tbl.pivot(index = 'HYD_AREA', columns = 'UNGULATES', values = 'SUM_ACRES')
#aspen_df.columns
#adf = aspen_df.rename(columns = {0: "NoGDE", 1: "Impact1", 2: "Impact2"})
#adf = adf.fillna(0)
#
#adf['ImpactAcres'] = adf.apply(lambda row: row.Impact1 + row.Impact2, axis=1)
#
#tbl = tbl.groupby("UNGULATES").sum()
#
## Get sum of all acres
## These depend on dataframe indices - may not work coorectly all the time...
#acres_all = tbl.sum(axis = 0)[2]
#
## Get acres where ungulates = 0
#acres0 = tbl.iloc[0][2]
#
## Get percent of where ungulates != 0
#acres_impact = acres_all - acres0
#percent_impact = acres_impact/acres_all
#print(percent_impact)









