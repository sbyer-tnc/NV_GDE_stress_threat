#-------------------------------------------------------------------------------
# Name:        Compile invasive species data
# Purpose:     From full NAS database download, filter out the species we want to focus on as threats to NV GDEs
#              Clean and import occurrence data from EDDMaps.org
#              Find invasive species location data in the SSI database
#              Identify invasive species presence info from BLM AIM data
#
# Author:      sarah.byer
#
# Created:     March 2020

#-------------------------------------------------------------------------------

# Import ArcGIS modules and check out spatial analyst extension
import arcpy, os
from arcpy import env
from arcpy.sa import *
arcpy.CheckOutExtension("spatial")

# import other necessary modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import os
from matplotlib.dates import DateFormatter

# Path to temporary geodatabase
path =  r"E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Maps\GDE_Threats.gdb"

# Environment settings
env.workspace = path
env.overwriteOutput = True
env.outputCoordinateSystem = arcpy.SpatialReference(26911) # Spatial reference NAD 1983 UTM Zone 11N. The code is '26911'

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Format NAS data

# Read in file of ALL NAS occurrences
file = r"D:\GDE_Threats\Invasives\NAS\occurrence.csv"
nas = arcpy.TableToTable_conversion(file, path, 'nas_occurrence')
[f.name for f in arcpy.ListFields(nas)]

# List of species we want to focus on:
spp = ['Lithobates catesbeianus',
           'Procambarus clarkii',
           'Oreochromis aureus',
           'Pacifastacus leniusculus',
           'Faxonius virilis',
           'Esox lucius',
           'Micropterus dolomieu',
           'Dreissena rostriformis bugensis',
           'Dreissena polymorpha',
           'Potamopyrgus antipodarum',
#           'Tamarisk spp.',
#           'Lepidium draba',
#           'Centaurea spp.',
           'Myriophyllum spicatum']
len(spp)

# Filter for species in the list from the full NAS occurrence dataset
nas_spp = list()
for species in spp:
    with arcpy.da.SearchCursor(nas, 'scientificName') as cursor:
        for row in cursor:
            nas_spp.append(row[0])
#            if species in row[0]:
#                nas_spp.append(row[0])
    del cursor
nas_spp_u = set(nas_spp) # Unique list of all species in NAS dataset
for species in spp:
    if species in nas_spp:
        print(species)

#-------------------------------------------------------------------------------
# Filter out those species from NAS list

nas_filter = arcpy.Copy_management(nas, 'nas_occurrence_filter')
with arcpy.da.UpdateCursor(nas_filter, 'scientificName') as cursor:
    for row in cursor:
        if not(row[0] in spp):
            print("Deleting {}".format(row[0]))
            cursor.deleteRow()
del cursor
arcpy.GetCount_management(nas_filter)    

# Convert to point fc
nas_pt = arcpy.XYTableToPoint_management(nas_filter, "nas_filter", x_field = "decimalLongitude", y_field = "decimalLatitude")
arcpy.GetCount_management(nas_pt)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Read in datasets from EDDmaps.org

# Species that were not found in the NAS database
# Tamarisk
tam = arcpy.CopyFeatures_management(r"D:\GDE_Threats\Invasives\EDD\Tamarisk\points.shp", "edd_tamarisk")

# Lepidium (hoary cress... or white top???)
lep = arcpy.CopyFeatures_management(r"D:\GDE_Threats\Invasives\EDD\Lepidium\points.shp", "edd_lepidium")

# Cardaria (whitetop... or hoary cress???)
car = arcpy.CopyFeatures_management(r"D:\GDE_Threats\Invasives\EDD\Cardaria\points.shp", "edd_cardaria")

# Eurasian water mil-foil
eum = arcpy.CopyFeatures_management(r"D:\GDE_Threats\Invasives\EDD\EUMilfoil\points.shp", "edd_eumilfoil")

#-------------------------------------------------------------------------------
# Combine EDD datasets

# Need a dataset where all fields are long enough to store all records
fields = arcpy.ListFields(tam)
for field in fields:
    print("Field {0} is Type {1} with Length = {2}".format(field.name, field.type, field.length))

fms = arcpy.FieldMappings()
for field in fields:
    if (field.length < 50 and field.type == "String"):
        print("FIX MEEE!")
        fm = arcpy.FieldMap()
        fm.addInputField(tam, field.name)
        newfield = fm.outputField
        newfield.length = 255
        fm.outputField = newfield
        fms.addFieldMap(fm)
        
# Copy feature class with new field mappings
edd = arcpy.FeatureClassToFeatureClass_conversion(tam, path, 'edd_combo', field_mapping = fms)
        
# Append all datasets to create one
arcpy.Append_management([lep, car, eum], edd, "NO_TEST")
arcpy.GetCount_management(edd)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Combine EDD and NAS data

edd = 'edd_combo'
nas = 'nas_filter'

# Create text-version of the event date field from NAS data
arcpy.AddField_management(nas, "eventDateTxt", "TEXT")
arcpy.CalculateField_management(nas, "eventDateTxt", "str(!eventDate!)", "PYTHON3")

# Custom function for field mapping
def mapFields(inlayer, infield, mapfield_name, mapfield_type): # mapFields function
    fldMap = arcpy.FieldMap()
    #fldMap.addInputField(path + "\\" + str(inlayer), infield)
    fldMap.addInputField(inlayer, infield)
    mapOut = fldMap.outputField
    mapOut.name, mapOut.type = mapfield_name, mapfield_type
    fldMap.outputField = mapOut
    return fldMap

# Field mapping for EDD to NAS data
[f.name for f in arcpy.ListFields(edd)]
[f.name for f in arcpy.ListFields(nas)]
sciname_map = mapFields(edd, "SCINAME", "scientificName", "TEXT")
comname_map = mapFields(edd, "COMNAME", "vernacularName", "TEXT")
latmap = mapFields(edd, "LATITUDE", "decimalLatitude", "DOUBLE")
lonmap = mapFields(edd, "LONGITUDE", "decimalLongitude", "DOUBLE")
datemap = mapFields(edd, "OBSDATE", "eventDateTxt", "TEXT")
locmap = mapFields(edd, "LOCATION", "locality", "TEXT")
geomap = mapFields(edd, "RECSOURCE", "assocaitedReferences", "TEXT")
geomap2 = mapFields(edd, "RECSRCTYP", "georeferenceProtocol", "TEXT")
# MORE if needed....
# NAS data have good attributes

fldMap_list = [sciname_map, comname_map, latmap, lonmap, datemap, locmap, geomap, geomap2]
invFldMappings = arcpy.FieldMappings()
for fm in fldMap_list:
    invFldMappings.addFieldMap(fm)

# Make copy of nas layer and append EDD to it
combo = arcpy.Copy_management(nas, 'nas_edd_combo')
arcpy.Append_management(edd, combo, "NO_TEST", invFldMappings)
arcpy.GetCount_management(combo)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Get invasive species data from SSI database
# SCript originally form GDE database species data formatting...

# Format species data from SSI

springs = r"D:\GDE_Threats\Invasives\SSI_Apr_21_2019\Nevada_Springs_Apr_21_2019.gdb\Nevada_Springs_Apr_21_2019"

# Read in spring species tables to one table
spring_vert = arcpy.TableToTable_conversion(r"D:\GDE_Threats\Invasives\SSI_Apr_21_2019\Nevada_Springs_Apr_21_2019.gdb\\Nevada_Springs_Apr_21_2019_Summarized_TaxaVert_by_Site", env.workspace, "ssi_taxa_vert")


spring_invert = arcpy.TableToTable_conversion(r"D:\GDE_Threats\Invasives\SSI_Apr_21_2019\Nevada_Springs_Apr_21_2019.gdb\\Nevada_Springs_Apr_21_2019_Summarized_TaxaInvert_by_Site", env.workspace, "ssi_taxa_invert")


spring_flora = arcpy.TableToTable_conversion(r"D:\GDE_Threats\Invasives\SSI_Apr_21_2019\Nevada_Springs_Apr_21_2019.gdb\\Nevada_Springs_Apr_21_2019_Summarized_TaxaFlora_by_Site", env.workspace, "ssi_taxa_flora")


# Copy Vert table to which the Invert and Plant tables will be appended
spring_vert_copy = arcpy.Copy_management(spring_vert, "spring_species_merge")

# Map fullest species name field to the copied Vert species table and combine
def mapFields(inlayer, infield, mapfield_name, mapfield_alias, mapfield_type): # mapFields function
    fldMap = arcpy.FieldMap()
    fldMap.addInputField(inlayer, infield)
    mapOut = fldMap.outputField
    mapOut.name, mapOut.alias, mapOut.type = mapfield_name, mapfield_alias, mapfield_type
    fldMap.outputField = mapOut
    return fldMap
invert_name_map = mapFields(spring_invert, "FullName", "FaunaFullName", "FaunaFullName", "TEXT")
invert_comname_map = mapFields(spring_invert, "CommonName", "FaunaCommonName", "FaunaCommonName", "TEXT")
invert_site_map = mapFields(spring_invert, "SiteID", "SiteID", "SiteID", "TEXT")
invert_endemism_map = mapFields(spring_invert, "EndemismLevel", "EndemismLevel", "EndemismLevel", "LONG")
invert_order_map = mapFields(spring_invert, "Order_", "FaunaOrder", "FaunaOrder", "TEXT")
invert_family_map = mapFields(spring_invert, "Family", "FaunaFamily", "FaunaFamily", "TEXT")
invert_genus_map = mapFields(spring_invert, "Genus", "FaunaGenus", "FaunaGenus", "TEXT")
invert_species_map = mapFields(spring_invert, "Species", "FaunaSpecies", "FaunaSpecies", "TEXT")
plant_name_map = mapFields(spring_flora, "FloraSpecies", "FaunaFullName", "FaunaFullName", "TEXT")
plant_comname_map = mapFields(spring_flora, "FloraCommonName", "FaunaCommonName", "FaunaCommonName", "TEXT")
plant_site_map = mapFields(spring_flora, "SiteID", "SiteID", "SiteID", "TEXT")
plant_endemism_map = mapFields(spring_flora, "EndemismLevel", "EndemismLevel", "EndemismLevel", "LONG")
plant_family_map = mapFields(spring_flora, "Family", "FaunaFamily", "FaunaFamily", "TEXT")
plant_genus_map = mapFields(spring_flora, "Genus", "FaunaGenus", "FaunaGenus", "TEXT")
plant_species_map = mapFields(spring_flora, "Species", "FaunaSpecies", "FaunaSpecies", "TEXT")
fldMap_list = [invert_name_map, invert_site_map, invert_endemism_map, invert_order_map, invert_family_map, invert_genus_map, invert_species_map, plant_name_map, plant_site_map, plant_endemism_map, plant_family_map, plant_genus_map, plant_species_map]
allFldMappings = arcpy.FieldMappings()
for fm in fldMap_list:
    allFldMappings.addFieldMap(fm)
spring_species = arcpy.Append_management([spring_invert, spring_flora], spring_vert_copy, "NO_TEST", allFldMappings)
[f.name for f in arcpy.ListFields(spring_species)]
# Not all species have full names (ex. Fish, Magpie, Rattlesnake). These will not be included in final dataset


#-------------------------------------------------------------------------------
# Filter the species list for species of concern (above)
# Will likely need to use both scientific and common names, and search for parts of the names

species = "spring_species_merge"
[f.name for f in arcpy.ListFields(species)]

# Common names
com = ["bullfrog",
       "Red swamp crayfish",
       "Blue Tilapia",
       "Signal crayfish",
       "Northern crayfish",
       "Northern pike",
       "Small-mouth bass",
       "Quagga mussel",
       "Zebra mussel",
       "New Zealand mudsnail",
       "Tamarisk",
       "Hoary cress",
       "Tall whitetop",
       "Eurasian watermilfoil",
       "knapweed"]


# Filter common names from species list
species_common = arcpy.Copy_management(species, "species_commonnames")
with arcpy.da.UpdateCursor(species_common, ["FaunaCommonName"]) as cursor:
    for row in cursor:
        if row[0] in com:
            print(row[0])
#        else:
#            cursor.deleteRow()
del cursor
arcpy.GetCount_management(species_common)

# See if list contains any general species names from common name list
general = ["bullfrog", "crayfish", "mussel", "saltcedar", "salt cedar", "pike", "bass", "tilapia", "whitetop", "tamarix", "tamarisk"]
arcpy.AddField_management(species_common, "Invasive", "LONG")
for name in general:
    with arcpy.da.UpdateCursor(species_common, ["FaunaCommonName", "Invasive"]) as cursor:
        for row in cursor:
            if name.upper() in row[0].upper():
                print(row[0].upper()) 
                row[1] = 1
                cursor.updateRow(row)
#            else:
#                row[1] = 0
#                cursor.updateRow(row)
del cursor
# Only matches bullfrog, bass, and largemouth bass
# For now, consider all of these invasive
    


# Filter scientific names from species list

# Scientific names
spp = ['Lithobates catesbeianus', 'Rana catesbeiana',
           'Procambarus clarkii',
           'Oreochromis aureus',
           'Pacifastacus leniusculus',
           'Faxonius virilis',
           'Esox lucius',
           'Micropterus dolomieu',
           'Dreissena rostriformis bugensis',
           'Dreissena polymorpha',
           'Potamopyrgus antipodarum',
           'Tamarisk spp.',
           'Lepidium draba',
           'Centaurea spp.',
           'Myriophyllum spicatum']

### RULE ENACTED: Only species record with Genus + Species allowed to stay
# Create a new attribute to hold a logical scientific name for the species. FaunaFullName may = order + family + genus + species
# Standardize empty genus + species attributes
species_sci = arcpy.Copy_management(species_common, "species_scinames")
with arcpy.da.UpdateCursor(species_sci, ['FaunaGenus']) as cursor:
    for row in cursor:
        if (len(str(row[0])) < 2) or (row[0] is None):
            row[0] = None
            cursor.updateRow(row)
        else:
            print("{}".format(row[0]))
del cursor
with arcpy.da.UpdateCursor(species_sci, ['FaunaSpecies']) as cursor:
    for row in cursor:
        if (len(str(row[0])) < 2) or (row[0] is None):
            row[0] = None
            cursor.updateRow(row)
        else:
            print("{}".format(row[0]))
del cursor

# All "empty" rows in the genus/species attributes now have "Null"
arcpy.AddField_management(species_sci, "SciName", "TEXT")
with arcpy.da.UpdateCursor(species_sci, ['FaunaGenus', 'FaunaSpecies', 'SciName']) as cursor:
    for row in cursor:
        if row[1] is not None:
        # If Species field isn't empty, populate the sciname field    
            row[2] = str(row[0]) + " " + str(row[1])
            cursor.updateRow(row)
        else:
            row[2] = "NA"
            cursor.updateRow(row)
del cursor

# See if it contains any species from our list
for name in spp:
    with arcpy.da.UpdateCursor(species_sci, ["SciName", "Invasive"]) as cursor:
        for row in cursor:
            if name.upper() in row[0].upper():
                print(row[0].upper()) 
                row[1] = 1
                cursor.updateRow(row)
del cursor


#-------------------------------------------------------------------------------
# Join invasive species records to spring locations

invasives = arcpy.Copy_management(species_sci, "species_invasive")
with arcpy.da.UpdateCursor(invasives, 'Invasive') as cursor:
    for row in cursor:
        if row[0] is None:
            cursor.deleteRow()
del cursor
arcpy.GetCount_management(invasives)
[f.name for f in arcpy.ListFields(invasives)]

springs_invasive = arcpy.CopyFeatures_management(r"D:\GDE_Threats\Invasives\SSI_Apr_21_2019\Nevada_Springs_Apr_21_2019.gdb\Nevada_Springs_Apr_21_2019", "springs_invasive")

arcpy.JoinField_management(springs_invasive, 'SiteID', invasives, 'SiteID', ['FaunaCommonName', 'SciName'])
with arcpy.da.UpdateCursor(springs_invasive, ['SciName']) as cursor:
    for row in cursor:
        if row[0] is None:
            cursor.deleteRow()
        else:
            print(row[0])
del cursor
arcpy.GetCount_management(springs_invasive)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Join SSI invasive species points to EDD and NAS points

combo = arcpy.Copy_management('nas_edd_combo', 'nas_edd_ssi_combo')
[f.name for f in arcpy.ListFields(combo)]
arcpy.GetCount_management(combo)

ssi = 'springs_invasive'
[f.name for f in arcpy.ListFields(ssi)]

# Map SSI fields to EDD/NAS fields
def mapFields(inlayer, infield, mapfield_name, mapfield_type): # mapFields function
    fldMap = arcpy.FieldMap()
    fldMap.addInputField(inlayer, infield)
    mapOut = fldMap.outputField
    mapOut.name, mapOut.type = mapfield_name, mapfield_type
    fldMap.outputField = mapOut
    return fldMap

name_map = mapFields(ssi, "SciName", "scientificName", "TEXT")
comname_map = mapFields(ssi, "FaunaCommonName", "vernacularName", "TEXT")
ref_map = mapFields(ssi, "SiteID", "references", "TEXT")
#plant_site_map = mapFields(ssi, "SiteID", "SiteID", "SiteID", "TEXT")
#plant_endemism_map = mapFields(spring_flora, "EndemismLevel", "EndemismLevel", "EndemismLevel", "LONG")
#plant_family_map = mapFields(spring_flora, "Family", "FaunaFamily", "FaunaFamily", "TEXT")
#plant_genus_map = mapFields(spring_flora, "Genus", "FaunaGenus", "FaunaGenus", "TEXT")
#plant_species_map = mapFields(spring_flora, "Species", "FaunaSpecies", "FaunaSpecies", "TEXT")
fldMap_list = [name_map, comname_map, ref_map]
allFldMappings = arcpy.FieldMappings()
for fm in fldMap_list:
    allFldMappings.addFieldMap(fm)
    
combo = arcpy.Append_management([ssi], combo, "NO_TEST", allFldMappings)
arcpy.GetCount_management(combo)

# END

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Code layers for invasive species STRESSOR
# If NAS/EDD/SSI lands on GDE, value == 1, otherwise value == 0.1

# Again, half-mile buffer around species points
spp = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Maps\GDE_Threats.gdb\nas_edd_ssi_combo'
arcpy.GetCount_management(spp)
spp_buff = arcpy.Buffer_analysis(spp, 'spp_halfmile', '0.5 Mile', '', '', 'ALL')
spp_buff = r'spp_halfmile'
[f.name for f in arcpy.ListFields(spp_buff)]

# GDE layers
springs = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Springs'
wet = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Wetlands'
#phr = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Phreatophytes' #**need to explode**#
lp = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Lakes_Playas'
rs = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Rivers_Streams'
phr2 = r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\Hydrology\NV_iGDE_assess.gdb\Phreatophytes_Explode'


# Springs
arcpy.AddField_management(springs, 'InvSpp_Str', 'LONG')
arcpy.CalculateField_management(springs, 'InvSpp_Str', 0.1)
[f.name for f in arcpy.ListFields(springs)]
spr_select = arcpy.SelectLayerByLocation_management(springs, 'INTERSECT', spp_buff)
arcpy.CalculateField_management(spr_select, 'InvSpp_Str', 1)

# Wetlands
arcpy.AddField_management(wet, 'InvSpp_Str', 'LONG')
arcpy.CalculateField_management(wet, 'InvSpp_Str', 0.1)
[f.name for f in arcpy.ListFields(wet)]
wet_select = arcpy.SelectLayerByLocation_management(wet, 'INTERSECT', spp_buff)
arcpy.CalculateField_management(wet_select, 'InvSpp_Str', 1)

# Phreatophytes
arcpy.AddField_management(phr2, 'InvSpp_Str', 'LONG')
arcpy.CalculateField_management(phr2, 'InvSpp_Str', 0.1)
[f.name for f in arcpy.ListFields(phr2)]
phr_select = arcpy.SelectLayerByLocation_management(phr2, 'INTERSECT', spp_buff)
arcpy.CalculateField_management(phr_select, 'InvSpp_Str', 1)

# Lakes/playas
arcpy.AddField_management(lp, 'InvSpp_Str', 'LONG')
arcpy.CalculateField_management(lp, 'InvSpp_Str', 0.1)
[f.name for f in arcpy.ListFields(lp)]
lp_select = arcpy.SelectLayerByLocation_management(lp, 'INTERSECT', spp_buff)
arcpy.CalculateField_management(lp_select, 'InvSpp_Str', 1)

# Rivers/streams
arcpy.AddField_management(rs, 'InvSpp_Str', 'LONG')
arcpy.CalculateField_management(rs, 'InvSpp_Str', 0.1)
[f.name for f in arcpy.ListFields(rs)]
rs_select = arcpy.SelectLayerByLocation_management(rs, 'INTERSECT', spp_buff)
arcpy.CalculateField_management(rs_select, 'InvSpp_Str', 1)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Code layers for invasive species THREAT
# Road density as invasive species conduit

# This one doesn't cover entire state
#d = arcpy.Raster(r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\road_density.tif')

# Derived from 2014 TIGER road census data
# Available thru Esri
d = arcpy.Raster(r'E:\RCF Data Recovery\Recovered Files\P00\GDE_Threats\OtherData\road_density_tiger2014.tif')

# Normalize between 0 - 1
dmax = d.maximum
d_norm = d/dmax # divide by max value to normalize
d_norm.save('road_density_norm')
d_norm = arcpy.Raster(r'road_density_norm')

# Summarize for each GDE feature
# Springs
temp_spr = arcpy.sa.ExtractValuesToPoints(springs, d_norm, 'springs_road_density')
[f.name for f in arcpy.ListFields(temp_spr)]
arcpy.AddField_management(temp_spr, 'InvSpp_Thr', 'DOUBLE')
arcpy.CalculateField_management(temp_spr, 'InvSpp_Thr', '!RASTERVALU!', 'PYTHON3')
arcpy.JoinField_management(springs, 'OBJECTID', temp_spr, 'OBJECTID', ['InvSpp_Thr'])
[f.name for f in arcpy.ListFields(springs)]

# Wetlands
temp_wet = arcpy.sa.ZonalStatisticsAsTable(wet, 'OBJECTID', d_norm, 'wetlands_roads_density_tbl', 'DATA', 'MEAN')
[f.name for f in arcpy.ListFields(temp_wet)]
arcpy.AddField_management(temp_wet, 'InvSpp_Thr', 'DOUBLE')
arcpy.CalculateField_management(temp_wet, 'InvSpp_Thr', '!MEAN!', 'PYTHON3')
arcpy.JoinField_management(wet, 'OBJECTID', temp_wet, 'OBJECTID_1', ['InvSpp_Thr'])
[f.name for f in arcpy.ListFields(wet)]

# PHreatophytes
temp_phr = arcpy.sa.ZonalStatisticsAsTable(phr2, 'OBJECTID', d_norm, 'phr_roads_density_tbl', 'DATA', 'MEAN')
[f.name for f in arcpy.ListFields(temp_phr)]
arcpy.AddField_management(temp_phr, 'InvSpp_Thr', 'DOUBLE')
arcpy.CalculateField_management(temp_phr, 'InvSpp_Thr', '!MEAN!', 'PYTHON3')
arcpy.JoinField_management(phr2, 'OBJECTID', temp_phr, 'OBJECTID_1', ['InvSpp_Thr'])
[f.name for f in arcpy.ListFields(phr2)]

# Lakes/playas
temp_lp = arcpy.sa.ZonalStatisticsAsTable(lp, 'OBJECTID', d_norm, 'lp_roads_density_tbl', 'DATA', 'MEAN')
[f.name for f in arcpy.ListFields(temp_lp)]
arcpy.AddField_management(temp_lp, 'InvSpp_Thr', 'DOUBLE')
arcpy.CalculateField_management(temp_lp, 'InvSpp_Thr', '!MEAN!', 'PYTHON3')
arcpy.JoinField_management(lp, 'OBJECTID', temp_lp, 'OBJECTID_1', ['InvSpp_Thr'])
[f.name for f in arcpy.ListFields(lp)]

# Rivers/streams
temp_rs = arcpy.sa.ZonalStatisticsAsTable(rs, 'OBJECTID', d_norm, 'rs_roads_density_tbl', 'DATA', 'MEAN')
[f.name for f in arcpy.ListFields(temp_rs)]
arcpy.AddField_management(temp_rs, 'InvSpp_Thr', 'DOUBLE')
arcpy.CalculateField_management(temp_rs, 'InvSpp_Thr', '!MEAN!', 'PYTHON3')
arcpy.JoinField_management(rs, 'OBJECTID', temp_rs, 'OBJECTID_1', ['InvSpp_Thr'])
[f.name for f in arcpy.ListFields(rs)]









