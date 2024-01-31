
#%% imports
#from re import sub
import sys
import os
import datetime
import traceback
import arcpy
import platform
from arcpy.sa import *

# connect to netid
netid = os.getlogin()
if netid == 'bkgelder':
    sys.path.append(os.path.join('C:\\Users', netid, 'Box\\Data_Sharing\\Scripts\\basics'))
else:
    sys.path.append(os.path.join('C:\\Users', netid, 'Box\\Scripts\\basics'))

import dem_functions2 as df

#checkout extensions and set environment variables
arcpy.CheckOutExtension("spatial")
arcpy.env.overwriteOutput = True

start_time=datetime.datetime.strftime(datetime.datetime.now(), '%Y_%m_%d_%H_%M_%S')
print(f"{start_time}")

arcpy.env.overwriteOutput = True
arcpy.env.compression = "NONE"
arcpy.CheckOutExtension("spatial")
arcpy.CheckOutExtension("ImageAnalyst")

#check if the terminal is an interacive jupyter one and only save the first value of sys.argv if it is
run_info=sys.argv

for each in sys.argv:
    if "jupyter" in each:
        run_info=[sys.argv[0]]


#%%
# ***FUNCTIONS***

# Generic function to check the projection on a raster or featureclass against one provide, project into the desired one 
# if needed and then return the correct location no mater which case is true
 
def projection_fixer(input_raster_or_fc, input_WKID_string, processing_geodatabase):
    
    base_name=os.path.basename(arcpy.Describe(input_raster_or_fc).catalogPath)
    output_name=base_name+f"_{input_WKID_string}_reprojection"
    output_path=os.path.join(processing_geodatabase,output_name)
    projection_method=arcpy.Describe(input_raster_or_fc).spatialReference
    current_WKID_string=str(projection_method.factoryCode)
    
    #check if a reprojection of the fc/raster has already been created in the given geodatabase
    if not arcpy.Exists(output_path):

        #if a reprojection hasn't been done, check if the current projection is the one you want, otherwise reproject it
        if current_WKID_string == input_WKID_string:
            return_location=input_raster_or_fc
            print(f"{input_raster_or_fc} already in WKID: {input_WKID_string}")
                
        else:
            return_location=output_path
            if arcpy.Describe(input_raster_or_fc).dataType == "FeatureClass":
                arcpy.Project_management(input_raster_or_fc, output_path, input_WKID_string)  
            elif arcpy.Describe(input_raster_or_fc).dataType == "RasterDataset":
                arcpy.management.ProjectRaster(input_raster_or_fc,output_path,input_WKID_string)    
        
    else:
        print(f"{output_path} already exitsts")
        return_location=output_path

        

    return(return_location)


# Preprocessing functions are run on HUC12 scale data 

# Flow preprocessing uses unfilled, hole-punched DEMs
def preprocessing_flow_rasters(input_DEM,grid_order_threshold):
    with arcpy.EnvManager(cellSize=input_DEM, snapRaster=input_DEM):
        HUC_12_flow_direction = arcpy.sa.FlowDirection(in_surface_raster=input_DEM, force_flow="NORMAL",out_drop_raster='', flow_direction_type="D8")
        HUC_12_flow_accumulation=arcpy.sa.FlowAccumulation(in_flow_direction_raster=HUC_12_flow_direction, in_weight_raster="", data_type="FLOAT", flow_direction_type="D8")
        HUC_12_grid_order=arcpy.sa.StreamOrder(HUC_12_flow_accumulation,HUC_12_flow_direction)
    
        #Here I make the raster for areas of higher grid order than the provided threshold so later I can use an con statement and isnull to flip it.
        #I did it this way because flow initiation cells (0 flow accumulation) are considered null by the stream order tool and so later for flowpath 
        # limiting the initiation cells were getting dropped out
        HUC_12_above_grid_order_threshold=arcpy.sa.Con(in_conditional_raster=HUC_12_grid_order, in_true_raster_or_constant=1, in_false_raster_or_constant="", where_clause=f"VALUE > {grid_order_threshold}")
    
    return(HUC_12_flow_direction, HUC_12_flow_accumulation, HUC_12_grid_order, HUC_12_above_grid_order_threshold)


# used in forest mode to clip out the current HUC 12 from the whole landfire raster 
# and return a raster with a spatial reference matching the DEM provided

def preprocessing_cut_landfire_to_HUC(canopy_cover_raster, HUC_12_DEM, HUC_12_boundary, processing_geodatabase):
    #wkid:26925 = UTM zone 15  wkid:5070 = Albers

    HUC_12_string=os.path.basename(HUC_12_DEM).replace(".tif","").replace("ep3m","")
    canopy_sr=arcpy.Describe(canopy_cover_raster).spatialReference
    HUC_12_DEM_sr=arcpy.Describe(HUC_12_DEM).spatialReference
    canopy_sr_string=str(canopy_sr.factoryCode)
    HUC_12_DEM_sr_string=str(HUC_12_DEM_sr.factoryCode)
    converted_boundary_path=os.path.join(processing_geodatabase,os.path.basename(HUC_12_boundary))+"_"+canopy_sr_string
    converted_boundary=arcpy.management.Project(HUC_12_boundary,converted_boundary_path,canopy_sr)
    clipped_canopy_cover_raster_path=os.path.join(processing_geodatabase,"canopy_cover_"+HUC_12_string)+"_"+canopy_sr_string
    clipped_canopy_cover_raster=arcpy.management.Clip(in_raster=canopy_cover_raster, out_raster=clipped_canopy_cover_raster_path, in_template_dataset=converted_boundary,clipping_geometry="NONE")
    reprojected_clipped_canopy_cover_raster_path=os.path.join(processing_geodatabase,"canopy_cover_"+HUC_12_string)+"_"+HUC_12_DEM_sr_string
    with arcpy.EnvManager(cellSize=HUC_12_DEM, snapRaster=HUC_12_DEM): 
        reprojected_clipped_canopy_cover_raster=arcpy.management.ProjectRaster(clipped_canopy_cover_raster,reprojected_clipped_canopy_cover_raster_path,HUC_12_DEM_sr)


    return(reprojected_clipped_canopy_cover_raster)

def preprocessing_make_flow_endpoints(field_boundaries_feature, channels_feature, flow_accumulation_raster, flow_endpoints_raster_name, thresholded_grid_order, processing_geodatabase):
    
    #set up names, environment variables and rasters for possible flow end point inputs
    HUC_12_channel_raster="channels"
    HUC_12_field_boundaries_raster="fld_bndries"
    
    with arcpy.EnvManager(cellSize=flow_accumulation_raster, snapRaster=flow_accumulation_raster):
        projection_method=arcpy.Describe(flow_accumulation_raster).spatialReference
        
        arcpy.conversion.FeatureToRaster(channels_feature, 'LINKNO', HUC_12_channel_raster)
        arcpy.conversion.FeatureToRaster(field_boundaries_feature, 'isAG', HUC_12_field_boundaries_raster)
        
        #combine all the possible enpoints: Null flow accumulation (flow spots that have weird issues where they form a circle), 
        # null land use (land that isn't coded as either ag nor non-ag), 
        # grid order threshold (land that is above the input threshold for distiguishing between sheet and rill and gulley erosion),
        # channel cells,
        # and cells coded as non-ag land

        null_flow_accumulation=arcpy.sa.Con(IsNull(flow_accumulation_raster), 1)
        null_land_use=arcpy.sa.Con(IsNull(HUC_12_field_boundaries_raster), 1)
        HUC_12_ag_land_raster=arcpy.sa.Con(in_conditional_raster=HUC_12_field_boundaries_raster, in_true_raster_or_constant=1, in_false_raster_or_constant="", where_clause="VALUE >= 1")
        HUC_12_non_ag_land_raster=arcpy.sa.Con(in_conditional_raster=HUC_12_field_boundaries_raster, in_true_raster_or_constant=1, in_false_raster_or_constant="", where_clause="VALUE = 0")
        arcpy.management.CreateMosaicDataset(processing_geodatabase, flow_endpoints_raster_name, projection_method)
        arcpy.management.AddRastersToMosaicDataset(in_mosaic_dataset=flow_endpoints_raster_name, raster_type= "Raster Dataset", input_path=HUC_12_channel_raster)
        arcpy.management.AddRastersToMosaicDataset(in_mosaic_dataset=flow_endpoints_raster_name, raster_type= "Raster Dataset", input_path=HUC_12_non_ag_land_raster)
        arcpy.management.AddRastersToMosaicDataset(in_mosaic_dataset=flow_endpoints_raster_name, raster_type= "Raster Dataset", input_path=null_flow_accumulation)
        arcpy.management.AddRastersToMosaicDataset(in_mosaic_dataset=flow_endpoints_raster_name, raster_type= "Raster Dataset", input_path=null_land_use)
        arcpy.management.AddRastersToMosaicDataset(in_mosaic_dataset=flow_endpoints_raster_name, raster_type= "Raster Dataset", input_path=thresholded_grid_order)

    return(flow_endpoints_raster_name, HUC_12_ag_land_raster, HUC_12_non_ag_land_raster)


def preprocessing_make_forest_flow_endpoints(HUC_12_canopy_cover_raster, HUC_12_boundary, channels_feature, OSM_roads_feature,OSM_railroads_feature,
                                             OSM_waterways_feature, OSM_water_feature, USFS_roads_feature, USFS_trails_feature,
                                             state_forest_roads_feature, flow_accumulation_raster, forest_flow_endpoints_raster_name,
                                             thresholded_grid_order, processing_geodatabase):
    
    #set up names, environment variables and rasters for possible flow end point inputs
    HUC_12_channel_raster="channels"


    files_to_clip=[OSM_roads_feature,OSM_railroads_feature,OSM_waterways_feature, OSM_water_feature, USFS_roads_feature, USFS_trails_feature,state_forest_roads_feature]

    with arcpy.EnvManager(scratchWorkspace=processing_geodatabase, workspace=processing_geodatabase,cellSize=flow_accumulation_raster, snapRaster=flow_accumulation_raster):
        
        projection_method=arcpy.Describe(flow_accumulation_raster).spatialReference

        # make and combine all the possible enpoints: 
        # roads (OSM, Forest service, state forest)
        # railways (OSM)
        # waterways (OSM)
        # water (OSM)
        # 
        # channel cells,
        # non-forested land,
        # Null flow accumulation (flow spots that have weird issues where they form a circle),  
        # grid order threshold (land that is above the input threshold for distiguishing between sheet and rill and gulley erosion),


        arcpy.conversion.FeatureToRaster(channels_feature, 'LINKNO', HUC_12_channel_raster)
        
        null_flow_accumulation=arcpy.sa.Con(IsNull(flow_accumulation_raster), 1)
        HUC_12_non_forest_land_raster=SetNull(HUC_12_canopy_cover_raster,1,"CC_PERCENT <> 'Non-Forested'")
        HUC_12_forest_land_raster=SetNull(HUC_12_canopy_cover_raster,1,"CC_PERCENT = 'Non-Forested'")
        
        arcpy.management.CreateMosaicDataset(processing_geodatabase, forest_flow_endpoints_raster_name, projection_method)
        arcpy.management.AddRastersToMosaicDataset(in_mosaic_dataset=forest_flow_endpoints_raster_name, raster_type= "Raster Dataset", input_path=HUC_12_channel_raster)
        arcpy.management.AddRastersToMosaicDataset(in_mosaic_dataset=forest_flow_endpoints_raster_name, raster_type= "Raster Dataset", input_path=HUC_12_non_forest_land_raster)
        arcpy.management.AddRastersToMosaicDataset(in_mosaic_dataset=forest_flow_endpoints_raster_name, raster_type= "Raster Dataset", input_path=null_flow_accumulation)
        arcpy.management.AddRastersToMosaicDataset(in_mosaic_dataset=forest_flow_endpoints_raster_name, raster_type= "Raster Dataset", input_path=thresholded_grid_order)

        for each_file in files_to_clip:
            
            fc_save_name=os.path.join(processing_geodatabase,os.path.basename(each_file)+"_clipped")
            raster_save_name=os.path.join(processing_geodatabase,os.path.basename(each_file)+"_raster")
            arcpy.analysis.Clip(each_file,HUC_12_boundary,fc_save_name)
            if int(arcpy.management.GetCount(fc_save_name)[0]) > 0:
                arcpy.conversion.FeatureToRaster(fc_save_name, "OBJECTID" , raster_save_name)
                arcpy.management.AddRastersToMosaicDataset(in_mosaic_dataset=forest_flow_endpoints_raster_name, raster_type= "Raster Dataset", input_path=raster_save_name)
            
            else:
                print(f"{each_file} has no intersection with {HUC_12_boundary}")

    return(forest_flow_endpoints_raster_name,HUC_12_forest_land_raster,HUC_12_non_forest_land_raster)
    
    #arcpy.management.GetCount()


def create_subcatchment_selection_list(subcatchments_feature, ag_land_raster, channels_feature, processing_geodatabase,sorting_type="TOTAL", subcatchments_upper_limit=100, min_percent_ag_land_threshold=None):
    with arcpy.EnvManager(workspace=processing_geodatabase,scratchWorkspace=processing_geodatabase):
        try:
        
            # calculate ag area in each subcatchment by watershed number
            subcatchment_ag_area = TabulateArea(subcatchments_feature, "WSNO", ag_land_raster, "Value", "ag_land_area")

            # create and join features for processing. subcatchment feature is needed for calculating ag land percent and channel feature for slope
            subcatchment_layer = arcpy.MakeFeatureLayer_management(subcatchments_feature, "subcatchment_lyr")
            channel_layer = arcpy.MakeFeatureLayer_management(channels_feature, "channel_lyr")
            arcpy.AddJoin_management(subcatchment_layer, "WSNO", subcatchment_ag_area, "WSNO")
            arcpy.AddJoin_management(subcatchment_layer, "WSNO", channel_layer, "WSNO")

            
            # selector sting names
            subcatchment_area = os.path.basename(subcatchments_feature) + ".Shape_Area"
            subcatchment_slope = os.path.basename(channels_feature) + ".Slope"
            subcatchment_watershed_number = os.path.basename(subcatchments_feature) + ".WSNO"

            # create selctor string to set limits on subcatchments

            selector_string = f"{subcatchment_area} >= 60702.9 AND {subcatchment_slope} > 0"
        
            #create empty lists and populate it with subcatchments that fit the criteria in the selector string
            subcatchment_sorting = []
            selected_subcatchments = []

            with arcpy.da.SearchCursor(subcatchment_layer, ["ag_land_area.VALUE_1", subcatchment_watershed_number, subcatchment_area], selector_string) as scur:
                for srow in scur:
                    subcatchment_sorting.append(srow)

            #get rid of subcatchments with no value for ag_land_area.VALUE_1
            subcatchment_sorting=[subcatchment for subcatchment in subcatchment_sorting if subcatchment[0] is not None]

            
            if sorting_type.upper() == "TOTAL":
                print("Subcatchments sorted by total ag land before thresholding")
                subcatchment_sorting=sorted(subcatchment_sorting, reverse=True)
                print(subcatchment_sorting)

                if len(subcatchment_sorting) > int(subcatchments_upper_limit):
                    subcatchment_sorting=subcatchment_sorting[:subcatchments_upper_limit]


            elif sorting_type.upper() == "PERCENT":
                print("Subcatchments sorted by percent ag land before thresholding")
                subcatchment_sorting=[[subcatchment[0]/subcatchment[2], subcatchment[1]] for subcatchment in subcatchment_sorting]    
                subcatchment_sorting=sorted(subcatchment_sorting, reverse=True)
                print(subcatchment_sorting)
                
                if min_percent_ag_land_threshold is not None:
                    print("Subcatchments selected based on percent ag land value being over threshold")
                    subcatchment_sorting=[subcatchment for subcatchment in subcatchment_sorting if subcatchment[0]>=min_percent_ag_land_threshold]

                            
                elif len(subcatchment_sorting)> int(subcatchments_upper_limit):
                    print("Subcatchments selected based on total number")
                    subcatchment_sorting=subcatchment_sorting[:subcatchments_upper_limit]
                    

            else:
                print("threshold type not properly defined")

            for subcatchment in subcatchment_sorting:
                selected_subcatchments.append(subcatchment[1])
            
            
        except:
            print("something went wrong with subcatchment selection")

    return(selected_subcatchments)


def create_forest_subcatchment_selection_list(subcatchments_feature, forest_land_raster, channels_feature, processing_geodatabase, sorting_type="TOTAL", subcatchments_upper_limit=100, min_percent_forest_land_threshold=None):
    with arcpy.EnvManager(workspace=processing_geodatabase,scratchWorkspace=processing_geodatabase):
        try:
        
            # calculate forest area in each subcatchment by watershed number
            subcatchment_forest_area = TabulateArea(subcatchments_feature, "WSNO", forest_land_raster, "Value", "forest_land_area")
            
            
            # create and join features for processing. subcatchment feature is needed for calculating ag land percent and channel feature for slope
            subcatchment_layer = arcpy.MakeFeatureLayer_management(subcatchments_feature, "subcatchment_lyr")
            channel_layer = arcpy.MakeFeatureLayer_management(channels_feature, "channel_lyr")
            arcpy.AddJoin_management(subcatchment_layer, "WSNO", subcatchment_forest_area, "WSNO")
            arcpy.AddJoin_management(subcatchment_layer, "WSNO", channel_layer, "WSNO")

            
            # selector sting names
            subcatchment_area = os.path.basename(subcatchments_feature) + ".Shape_Area"
            subcatchment_slope = os.path.basename(channels_feature) + ".Slope"
            subcatchment_watershed_number = os.path.basename(subcatchments_feature) + ".WSNO"

            # create selctor string to set limits on subcatchments

            selector_string = f"{subcatchment_area} >= 60702.9 AND {subcatchment_slope} > 0"

            #create empty lists and populate it with subcatchments that fit the criteria in the selector string
            subcatchment_sorting = []
            selected_subcatchments = []

            with arcpy.da.SearchCursor(subcatchment_layer, ["forest_land_area.VALUE_1", subcatchment_watershed_number, subcatchment_area], selector_string) as scur:
                for srow in scur:
                    subcatchment_sorting.append(srow)

            
            subcatchment_sorting=[subcatchment for subcatchment in subcatchment_sorting if subcatchment[0] is not None]

            
            if sorting_type.upper() == "TOTAL":
                print("Subcatchments sorted by total forest land before thresholding")
                subcatchment_sorting=sorted(subcatchment_sorting, reverse=True)
                print(subcatchment_sorting)

                if len(subcatchment_sorting) > int(subcatchments_upper_limit):
                    subcatchment_sorting=subcatchment_sorting[:int(subcatchments_upper_limit)]


            elif sorting_type.upper() == "PERCENT":
                print("Subcatchments sorted by percent forest land before thresholding")
                subcatchment_sorting=[[subcatchment[0]/subcatchment[2], subcatchment[1]] for subcatchment in subcatchment_sorting]    
                subcatchment_sorting=sorted(subcatchment_sorting, reverse=True)
                print(subcatchment_sorting)
                
                if min_percent_forest_land_threshold is not None:
                    print("Subcatchments selected based on percent ag land value being over threshold")
                    subcatchment_sorting=[subcatchment for subcatchment in subcatchment_sorting if subcatchment[0]>=min_percent_forest_land_threshold]

                            
                elif len(subcatchment_sorting)> int(subcatchments_upper_limit):
                    print("Subcatchments selected based on total number")
                    subcatchment_sorting=subcatchment_sorting[:int(subcatchments_upper_limit)]
                    

            else:
                print("threshold type not properly defined")

            for subcatchment in subcatchment_sorting:
                selected_subcatchments.append(subcatchment[1])
            
            
        except:
            print("something went wrong with subcatchment selection")

    return(selected_subcatchments)




#%%

processing_projection_dict={"agriculture":{"26915":["071000081505"]}}#"forest":{"26915":["090300010705","090300010707","090300010708"]}}#,"26914":["090300010705","090300010707","090300010708"]}}

for terrain_type, WKID_HUC_dict in processing_projection_dict.items():
    for current_WKID, current_HUC_12s_list in WKID_HUC_dict.items():
        for current_HUC_12 in current_HUC_12s_list:
            
            parameters = ["C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/pythonw.exe",
            "O:/DEP/Scripts/basics/cmd_FlowPath_v9.py",
            platform.node(),
            "2022",
            current_HUC_12,
            current_WKID,
            "mean18",
            "3"]

            if len(sys.argv) == 1:
                cleanup = False
                for i in parameters[2:]:
                    sys.argv.append(i)
                
            else:
                cleanup = True
                for counter,current_parameter in enumerate(parameters[2:]):
                    sys.argv[counter+1]=current_parameter
                    

            ####            node = sys.argv[1]
            ACPFyear = sys.argv[2]
            huc12 = sys.argv[3]
            srOutCode = sys.argv[4]
            interpType = sys.argv[5]
            cellSize = int(sys.argv[6])
            # print(f"parameters are {parameters}")
            print(f"sys.argvs are {sys.argv}")
            # print(f"srOutCode is {srOutCode}")

            nowYmd = datetime.datetime.strftime(datetime.datetime.now(), '%Y_%m_%d_%H_%M_%S')

            location_dictionary = {}
            for location_dictionary_key,location_dictionary_path in df.loadVariablesDict(platform.node(), ACPFyear, huc12, srOutCode, interpType, cellSize, nowYmd).items():
                if '\\\\EL3354-02\\O$' in location_dictionary_path:
                    location_dictionary[location_dictionary_key]=location_dictionary_path.replace('\\\\EL3354-02\\O$', 'M:')
                elif "\\\\EL3354-02\\D$" in location_dictionary_path:
                    location_dictionary[location_dictionary_key]=location_dictionary_path.replace('\\\\EL3354-02\\D$', 'D:')
                else:
                    location_dictionary[location_dictionary_key]=location_dictionary_path

            location_dictionary_versioned = {}
            for location_dictionary_key,location_dictionary_path in df.loadVariablesDict(platform.node(), ACPFyear, huc12, srOutCode, interpType, cellSize, nowYmd, 'tim_random_flowpath_update_testing').items():
                if '\\\\EL3354-02\\O$' in location_dictionary_path:
                    location_dictionary_versioned[location_dictionary_key]=location_dictionary_path.replace('\\\\EL3354-02\\O$', 'M:')
                elif "\\\\EL3354-02\\D$" in location_dictionary_path:
                    location_dictionary_versioned[location_dictionary_key]=location_dictionary_path.replace('\\\\EL3354-02\\D$', 'D:')
                else:
                    location_dictionary_versioned[location_dictionary_key]=location_dictionary_path
            #%%

            optional_parameters= [4,0.25,1985,"TOTAL"]

            #%%


            #check if a processing directory and geodatabase exist and if they don't, make them
             
            processing_directory_path=os.path.join(location_dictionary_versioned['depBase'],os.path.join(terrain_type.lower(),os.path.join(current_HUC_12[:8],current_HUC_12)))
            print(f"Processing in {processing_directory_path}")

            processing_gdb=os.path.join(processing_directory_path,"processing.gdb")

            if not os.path.exists(processing_directory_path):
                os.makedirs(processing_directory_path)
                arcpy.management.CreateFileGDB(processing_directory_path, "processing.gdb")

            elif not os.path.exists(processing_gdb):
                arcpy.management.CreateFileGDB(processing_directory_path, "processing.gdb")

            #set workspace and scratch to processing geodatabase
            arcpy.env.workspace=processing_gdb
            arcpy.env.scratchWorkspace=processing_gdb

            #set up input names while also checking they are in the correct projection
            if terrain_type.upper()=="AGRICULTURE":
                punched_dem_path=projection_fixer(location_dictionary["pElevFile"],current_WKID,processing_gdb)
                field_boundaries_path=projection_fixer(location_dictionary["fieldBoundaries"],current_WKID,processing_gdb)
                subcatchment_boundaries_path=projection_fixer(location_dictionary["pdCatch"],current_WKID,processing_gdb)
                channel_path=projection_fixer(location_dictionary["pdChnl"],current_WKID,processing_gdb)
                
                
                grid_order_threshold_value=optional_parameters[0]
                ag_land_threshold_value=optional_parameters[1]
                random_generator_seed_number=optional_parameters[2]
                subcatchment_ag_land_sorting_method=optional_parameters[3]
                
                

            elif terrain_type.upper()=="FOREST":
                punched_dem_path=projection_fixer(location_dictionary["pElevFile"],current_WKID,processing_gdb)         
                subcatchment_boundaries_path=projection_fixer(location_dictionary["pdCatch"],current_WKID,processing_gdb)
                channel_path=projection_fixer(location_dictionary["pdChnl"],current_WKID,processing_gdb)
                
                OSM_roads=projection_fixer(location_dictionary["roadsfc"],current_WKID,processing_gdb)
                OSM_railroads=projection_fixer(location_dictionary["rrsfc"],current_WKID,processing_gdb)
                OSM_waterways=projection_fixer(location_dictionary["waterwaysfc"],current_WKID,processing_gdb)
                OSM_water=projection_fixer(location_dictionary["waterfc"],current_WKID,processing_gdb)
                
                # USFS_roads=os.path.join(os.path.dirname(location_dictionary["roadsfc"]),"USFS_forest_roads")
                # USFS_trails=os.path.join(os.path.dirname(location_dictionary["roadsfc"]),"USFS_forest_trails")
                
                # state_abbreviation_string="MN"
                # state_FS_roads=os.path.join(os.path.dirname(location_dictionary["roadsfc"]),state_abbreviation_string+"_forest_roads")

                # temporary paths for roads until we get them worked into the update structure correctly
                temp_local_forest_roads_gdb=r"C:\Users\idep2\Desktop\tim_stuff\forest_roads_26915.gdb"
                
                USFS_roads=projection_fixer(os.path.join(temp_local_forest_roads_gdb,"USFS_forest_roads"),current_WKID,processing_gdb)
                USFS_trails=projection_fixer(os.path.join(temp_local_forest_roads_gdb,"USFS_forest_trails"),current_WKID,processing_gdb)
                
                state_abbreviation_string="MN"
                state_FS_roads=projection_fixer(os.path.join(temp_local_forest_roads_gdb, state_abbreviation_string+"_forest_roads"),current_WKID,processing_gdb)

                landfire_canopy_raster=location_dictionary["canopyCoverMap"]
                buffered_HUC_12_boundary=projection_fixer(location_dictionary["bufferedBoundaries"],current_WKID,processing_gdb)

                
                grid_order_threshold_value=optional_parameters[0]
                ag_land_threshold_value=optional_parameters[1]
                random_generator_seed_number=optional_parameters[2]
                subcatchment_ag_land_sorting_method=optional_parameters[3]

            
            else:
                print("Processing mode not properly defined")



            


            #%%
            with arcpy.EnvManager(scratchWorkspace=processing_gdb, workspace=processing_gdb):
                
                arcpy.env.snapRaster=punched_dem_path
                arcpy.env.cellSize=punched_dem_path

                #preprocess the HUC12 level raw inputs
                if terrain_type.upper()=="AGRICULTURE":
                    #%%
                    flow_direction_raster, flow_accumulation_raster, grid_order_raster, grid_order_threshold_raster = preprocessing_flow_rasters(input_DEM=punched_dem_path,grid_order_threshold=grid_order_threshold_value)
                    #%%
                    flow_end, ag_land, non_ag = preprocessing_make_flow_endpoints(field_boundaries_feature=field_boundaries_path,channels_feature=channel_path,flow_accumulation_raster=flow_accumulation_raster, flow_endpoints_raster_name="flow_ends", thresholded_grid_order=grid_order_threshold_raster, processing_geodatabase=processing_gdb)
               

                
                elif terrain_type.upper()=="FOREST":
                    forest_canopy_cover_path=preprocessing_cut_landfire_to_HUC(canopy_cover_raster=landfire_canopy_raster,HUC_12_DEM=punched_dem_path,HUC_12_boundary=buffered_HUC_12_boundary,processing_geodatabase=processing_gdb)
                    flow_direction_raster, flow_accumulation_raster, grid_order_raster, grid_order_threshold_raster = preprocessing_flow_rasters(input_DEM=punched_dem_path,grid_order_threshold=grid_order_threshold_value)
                    #%%
                    flow_end, forest_land, non_forest_land = preprocessing_make_forest_flow_endpoints(HUC_12_boundary=buffered_HUC_12_boundary, HUC_12_canopy_cover_raster=forest_canopy_cover_path,channels_feature=channel_path, OSM_roads_feature=OSM_roads,OSM_railroads_feature=OSM_railroads,OSM_water_feature=OSM_water,OSM_waterways_feature=OSM_waterways,USFS_roads_feature=USFS_roads,USFS_trails_feature=USFS_trails,state_forest_roads_feature=state_FS_roads, flow_accumulation_raster=flow_accumulation_raster, forest_flow_endpoints_raster_name="forest_flow_ends", thresholded_grid_order=grid_order_threshold_raster, processing_geodatabase=processing_gdb)
                    #%%
                

                try:
                    flow_endpoints_raster_name=f"flow_endpoints_{current_HUC_12}.tif"
                    flow_endpoints_raster_path=os.path.join(processing_directory_path,flow_endpoints_raster_name)
                    arcpy.MosaicToNewRaster_management(input_rasters=flow_end, output_location=processing_directory_path, raster_dataset_name_with_extension=flow_endpoints_raster_name, pixel_type="32_BIT_FLOAT",number_of_bands=1)
                    flow_ends_fixed=arcpy.sa.Con(in_conditional_raster=flow_endpoints_raster_path, in_true_raster_or_constant=1, where_clause="VALUE > 0")
                    flow_ends_fixed.save(flow_endpoints_raster_path)
                except:
                    print("Mosaic to new raster failed on flow endpoints mosaic dataset")

                