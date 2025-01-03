# NV_GDE_stress_threat
**Scripts used to estimate stressors and threats to groundwater dependent ecosystems in Nevada**

Pre-requisites for running all scripts:
- Nevada iGDE database (geodatabase or shapefiles)
- R/RStudio
- Python interpreter
- arcpy module
- other packages/modules loaded in scripts

Recommended order for running scripts to replicate entire analysis:
1. GDE_threats_level_trends.R (calculate trends in groundwater levels for groundwater stress from withdrawal/pumping)
2. GDE_StressThreat_Withdrawals.ipynb (consoliate groundwater trend data and basin commitments to calculate stress and threat to GDEs from pumping/appropriation)
3. GDE_Threats_Climate.ipynb (apply modeled climate threat to GDEs)
4. GDE_Stress_Ungulates.ipynb (apply ungulate impacts to GDEs as a stressor)
5. GDE_StressThreat_Species.ipynb (apply non-native species distributions and potential to GDEs)
6. GDE_StressThreat_Humans.ipynb (apply human activity like urbanization and diversions to GDEs as stressors and threats)
7. GDE_StressThreat_Code.ipynb (consoliate all stress and threat values across GDEs; create bivariate field to visualize both stress and threat levels)
8. GDE_StressThreat_GDB.ipynb (format field names and aliases for stressor and threat data in a new geodatabase)