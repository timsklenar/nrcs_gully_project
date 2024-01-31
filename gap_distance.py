# -*- coding: utf-8 -*-
#用于计算每个指标，在每个阈值取值下的沟道最端点（预测的沟头位置点）和沟头点的距离 It is used to calculate the distance between the end point of the EG (predictedhead position point) and the EG head point under each threshold value for each indicator
# 分两种距离，高估和低估。高估情况需要做沟头汇水区，选汇水区内最长距离；There are two kinds of distance, overestimation and underestimation. The overestimated situation requires that the gully head catchment area be built and the longest distance within the catchment area be selected;
# 低估情况下求沟头到沟道端点最近距离 Find the nearest distance from the trench head to the end of the trench in case of underestimation
#思路：王春梅；王春梅搭MODELBUILDER框架，刘欣实现PYTHON代码下的循环。Chunmei had the idea of the calculation and tried an MODELBUILDER, and then Liu Xin (my master student) finished this Python script.
#2020.4.15

import arcpy
from arcpy import management
from arcpy import da
from arcpy.sa import *
import os
import platform
import datetime
import pandas
import numpy 
import os
import shutil
import sys
import json
import traceback

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

#check if the terminal is an interacive jupyter one and only save the first value of sys.argv if it is
run_info=sys.argv

for each in sys.argv:
    if "jupyter" in each:
        run_info=[sys.argv[0]]

start_time=datetime.datetime.strftime(datetime.datetime.now(), '%Y_%m_%d_%H_%M_%S')
print(f"{start_time}")

#changes the paths in brian's dict builder to local server paths
def LocationDictionaryFixer(year_string,WKID,HUC12,DEM_type_string,DEM_resolution_string,versioning_string,arguments_list):
    parameters = ["C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/pythonw.exe",
            "O:/DEP/Scripts/basics/cmd_FlowPath_v9.py",
            platform.node(),
            year_string,
            HUC12,
            WKID,
            DEM_type_string,
            DEM_resolution_string]

    if len(arguments_list) == 1:
        cleanup = False
        for i in parameters[2:]:
            arguments_list.append(i)
        
    else:
        cleanup = True
        for counter,current_parameter in enumerate(parameters[2:]):
            arguments_list[counter+1]=current_parameter
            

    ####            node = arguments_list[1]
    ACPFyear = arguments_list[2]
    huc12 = arguments_list[3]
    srOutCode = arguments_list[4]
    interpType = arguments_list[5]
    cellSize = int(arguments_list[6])
    # print(f"parameters are {parameters}")
    print(f"arguments are {arguments_list}")
    # print(f"srOutCode is {srOutCode}")

    nowYmd = datetime.datetime.strftime(datetime.datetime.now(), '%Y_%m_%d_%H_%M_%S')

    fixed_dictionary = {}
    for location_dictionary_key,location_dictionary_path in df.loadVariablesDict(platform.node(), ACPFyear, huc12, srOutCode, interpType, cellSize, nowYmd).items():
        if '\\\\EL3354-02\\O$' in location_dictionary_path:
            fixed_dictionary[location_dictionary_key]=location_dictionary_path.replace('\\\\EL3354-02\\O$', 'M:')
        elif "\\\\EL3354-02\\D$" in location_dictionary_path:
            fixed_dictionary[location_dictionary_key]=location_dictionary_path.replace('\\\\EL3354-02\\D$', 'D:')
        else:
            fixed_dictionary[location_dictionary_key]=location_dictionary_path

    fixed_versioned_dictionary = {}
    for location_dictionary_key,location_dictionary_path in df.loadVariablesDict(platform.node(), ACPFyear, huc12, srOutCode, interpType, cellSize, nowYmd, versioning_string).items():
        if '\\\\EL3354-02\\O$' in location_dictionary_path:
            fixed_versioned_dictionary[location_dictionary_key]=location_dictionary_path.replace('\\\\EL3354-02\\O$', 'M:')
        elif "\\\\EL3354-02\\D$" in location_dictionary_path:
            fixed_versioned_dictionary[location_dictionary_key]=location_dictionary_path.replace('\\\\EL3354-02\\D$', 'D:')
        else:
            fixed_versioned_dictionary[location_dictionary_key]=location_dictionary_path
    return(fixed_dictionary,fixed_versioned_dictionary)

def GapDistancePreprocessing(raw_EG_heads,EG_heads,projection_WKID,flow_accumulation,error_radius,snapped_EG_raster_name,snapped_EG_feature_name,flow_direction_raster,upstream_watersheds_raster,downstream_paths_raster,full_output_gdb,valid_prediction_area_raster,longest_upstream_flowpaths_feature,longest_downstream_flowpaths_feature):
    
    #project points into local NAD 1983 UTM zone 14/15 N 
    print("reprojecting...")
    arcpy.Project_management(raw_EG_heads, EG_heads, projection_WKID)    
    
    #add a field for an integer version of the flowpath ID so it can persist through being turned into a raster and back (rasters don'e accept text values)
    print("making integer...")
    arcpy.AddField_management(EG_heads,"unique_flowpath_number","LONG")
    flowpath_IDs=["Path_ID","unique_flowpath_number"]
    
    #calculate the raster ID integers
    with da.UpdateCursor(EG_heads,flowpath_IDs) as cursor:
        for row in cursor:
            path_ID=row[0]
            # print(path_ID)
            path_split=path_ID.split("_")
            # print(path_split)
            path_as_number= int(path_split[1])*10000+int(path_split[0])
            # print(path_as_number)
            row[1]=path_as_number

            cursor.updateRow(row)


    #snap the gully heads to a raster pour point, save the raster, convert back to newly centered point feature carrying the unique flowpath intEGer through to table on the new feature
    print("snapping...")
    snapped_gully_heads_raster=SnapPourPoint(EG_heads,flow_accumulation,error_radius,"unique_flowpath_number")
    # print(snapped_EG_raster_name)
    snapped_gully_heads_raster.save(snapped_EG_raster_name)
    arcpy.RasterToPoint_conversion(snapped_gully_heads_raster,snapped_EG_feature_name,"VALUE")

    #change the default grid_code field name back to the original field name
    management.AlterField(snapped_EG_feature_name,"grid_code","unique_flowpath_number")

    #calculate the upstream watershed and downstream flowpath areas for the gully head pour points
    print("calculating flow areas...")
    watersheds_upstream_of_gullies = Watershed(flow_direction_raster,snapped_gully_heads_raster, "VALUE")
    watersheds_upstream_of_gullies.save(upstream_watersheds_raster)

    #calculate flowpaths downstream of gully heads
    downstream_paths=arcpy.sa.OptimalPathAsRaster(in_destination_data=snapped_gully_heads_raster, in_distance_accumulation_raster=flow_accumulation, in_back_direction_raster=flow_direction_raster)
    downstream_paths.save(downstream_paths_raster)
    
    #not sure this is necessary but set the values on the flow cells to 1
    # upstream_con=Con(Raster(upstream_watersheds_raster)>=0,Raster(upstream_watersheds_raster))
    # downstream_con=Con(Raster(downstream_paths_raster)>=0,1)
    # valid_prediction_area_list=[upstream_con,downstream_con]

    #create a new raster with both up and downstream flow areas combined to be a valid TI index selection mask
    #print("mosaicing...)
    print("mosaicing...is currently not working")
    # arcpy.MosaicToNewRaster_management(input_rasters=upstream_watersheds, output_location= full_output_gdb, raster_dataset_name_with_extension=os.path.basename(upstream_watersheds_raster), number_of_bands=1)
    # arcpy.MosaicToNewRaster_management(input_rasters=downstream_flowpaths, output_location= full_output_gdb, raster_dataset_name_with_extension=os.path.basename(downstream_paths_raster), pixel_type="32_BIT_UNSIGNED", number_of_bands=1,mosaic_method="MAXIMUM")
    arcpy.MosaicToNewRaster_management(input_rasters=[upstream_watersheds_raster, downstream_paths_raster], output_location= full_output_gdb, raster_dataset_name_with_extension=valid_prediction_area_raster, pixel_type="32_BIT_UNSIGNED", number_of_bands=1)

    # upstream_watersheds=[]
    # downstream_flowpaths=[]

    #calculate downstream flowpaths for each unique gully head flowpath ID and set the values of the resulting raster to the unique ID number

    with arcpy.EnvManager(mask=downstream_paths_raster):
        with da.SearchCursor(snapped_EG_feature_name,"unique_flowpath_number") as unique_flowpath_number_list:
            for each_number in unique_flowpath_number_list:
                
                #take the unique flowpath ID and get only the matching cell from the gully head raster
                current_number=each_number[0]
                print(f"calculating flow areas for gully head {current_number}")
                current_number_selection_string=f"VALUE <> {current_number}"
                current_gully_head_raster=SetNull(snapped_gully_heads_raster,snapped_gully_heads_raster,current_number_selection_string)
                
                #calculate the dowstream flowpath, then set the values of the downstream flowpath to the unique gullyhead ID number and add the resulting raster to the mosaicing name list
                current_downstream_flowpath=OptimalPathAsRaster(in_destination_data=current_gully_head_raster, in_distance_accumulation_raster=flow_accumulation, in_back_direction_raster=flow_direction_raster)
                current_downstream_flowpath_raster_name=os.path.join(full_output_gdb, f"downstream_flowpath_{current_number}")
                current_downstream_flowpath_with_id_value=Con(current_downstream_flowpath,current_number)
                current_downstream_flowpath_with_id_value.save(current_downstream_flowpath_raster_name)
                # downstream_flowpaths.append(current_downstream_flowpath_raster_name)
    
    # print(downstream_flowpaths)    
    
    #cut each upstream watershed out of the watersheds raster for each unique gully head flowpath ID 

    """this needs reworked to recalculate uptream watersheds each time to deal with the case of gully heads from nested watersheds"""
    with arcpy.EnvManager(mask=upstream_watersheds_raster):
        with da.SearchCursor(snapped_EG_feature_name,"unique_flowpath_number") as unique_flowpath_number_list:
            for each_number in unique_flowpath_number_list:
                
                #take the unique flowpath ID and get only the matching cell from the gully head raster
                current_number=each_number[0]
                print(f"calculating upstream watershed for gully head {current_number}")
                
                #calculate the upstream watershed for the gully head cell, save it and add the file to the mosaicing raster name list
                current_upstream_watershed_raster_name=os.path.join(full_output_gdb, f"upstream_watershed_{current_number}")
                current_watershed_upstream_of_gully_selection_string=f"VALUE <> {current_number}"
                current_watershed_upstream_of_gully = SetNull(upstream_watersheds_raster,upstream_watersheds_raster,current_watershed_upstream_of_gully_selection_string)
                current_watershed_upstream_of_gully.save(current_upstream_watershed_raster_name)
                
    #mosaic the upstream watershed and downstream flowpath raster for each unique ID number and save it 
    with da.SearchCursor(snapped_EG_feature_name,"unique_flowpath_number") as unique_flowpath_number_list:
        for each_number in unique_flowpath_number_list:
            
            #take the unique flowpath ID and get only the matching cell from the gully head raster
            current_number=each_number[0]
            print(f"mosaicing flow areas for gully head {current_number}")

            current_upstream_watershed_raster_name=os.path.join(full_output_gdb, f"upstream_watershed_{current_number}")
            current_downstream_flowpath_raster_name=os.path.join(full_output_gdb, f"downstream_flowpath_{current_number}")
            total_flow_area_raster_name=f"valid_flow_cells_{current_number}"

            arcpy.MosaicToNewRaster_management(input_rasters=[current_downstream_flowpath_raster_name, current_upstream_watershed_raster_name], output_location= full_output_gdb, raster_dataset_name_with_extension=total_flow_area_raster_name, pixel_type="32_BIT_UNSIGNED", number_of_bands=1)

    #calculate longest flowpaths for each watershed and downstream flowpath    
    with arcpy.EnvManager(mask=upstream_watersheds_raster):
        print("calculating longest upstream flowpath...")
        longest_upstream_flowpaths_raster=FlowLength(flow_direction_raster,"UPSTREAM")
        longest_upstream_points_feature=ExtractValuesToPoints(snapped_EG_feature_name,longest_upstream_flowpaths_raster,longest_upstream_flowpaths_feature)
        
    with arcpy.EnvManager(mask=downstream_paths_raster):
        print("calculating longest downstream flowpath...")
        longest_downstream_flowpaths_raster=FlowLength(flow_direction_raster_location,"DOWNSTREAM")
        longest_downstream_points_feature=ExtractValuesToPoints(snapped_EG_feature_name,longest_downstream_flowpaths_raster,longest_downstream_flowpaths_feature)
        
    return

def TerrainIndexThresholdListMaker(minimum_value, maximum_value, step_size, decimal_precision):
    current_ti_list=numpy.arange(minimum_value,maximum_value+(step_size*0.1),step_size).tolist()
    fixed_ti_list=[]
    for each in current_ti_list:
        each = round(each,decimal_precision)
        fixed_ti_list.append(each)
        
    return(fixed_ti_list)

    
# def GridOderPointsGenerator(watershed_id,terrain_index_type,threshold_value, thresholded_terrain_index_raster, flow_direction,output_gdb):
#     stream_string=str(watershed_id) + "_" + str(terrain_index_type) + "_streams_" + str(threshold_value)
#     stream_feature_name=os.path.join(output_gdb,stream_string)
#     StreamToFeature(thresholded_terrain_index_raster,flow_direction,stream_feature_name)
    
#     dissolve_string=str(watershed_id) + "_" + str(terrain_index_type) + "_streams_dissolved_" + str(threshold_value)
#     management.Dissolve(stream_feature_name,dissolve_string)

#     point_string=str(watershed_id) + "_" +str(terrain_index_type) + "_points_" + str(threshold_value)
    
#     return()


#current huc12s by projection WKID
huc_projection_dict = {"26914": ["102400060304","102702050101","102702070207"], "26915": ["070801030408","070801050302","070802050807", "071000040910", "102300030509","102300031003","102300031403"]} 
#"102702060102" data is missing in the file structure for mnmx but exists for mean as of 7/24/23

# huc_projection_dict = {"26915": ["070801050302","071000040910", "102300031003","102300031403"],"26914": ["102400060304","102702070207","102702050101","102702060102"]} 
# huc_projection_dict = {"26915": ["070801050302"],"26914": ["102400060304"]}
#huc_projection_dict = {"26915": ["070801050302"]}
# huc_projection_dict = {"26914": ["102400060304"], "26915": ["070802050807","071000040910"]} # {"26915": ["070802050807","071000040910"]}#
# huc_projection_dict = {"26914": ["102400060304","102702050101","102702070207"], "26915": ["070801030408","070801050302","070802050807", "071000040910", "102300030509","102300031003","102300031403"]} 


#setting up dictionaries to iterate over for threshold calculations
#from chunmei's paper
chunmei_terrain_index_log_thresholds={"specific_contributing_area":[2.6],"stream_power_index": [1.5],"modified_stream_power_index": [0.25], "compound_topographic_index":[1.2]}
#chunmei_terrain_index_log_thresholds={"stream_power_index": [1.5]}

specific_area_min=1.6
specific_area_max=3.6
specific_area_step=0.1

# spi_min=120
# spi_max=180
# spi_step=5

# modified_spi_min=15
# modified_spi_max=35
# modified_spi_step=1

# cti_min=90
# cti_max=150
# cti_step=5

# specific_area_threshold_list=numpy.arange(specific_area_min,specific_area_max+specific_area_step,specific_area_step).tolist()
# spi_threshold_list=numpy.arange(spi_min,spi_max+1,spi_step).tolist()
# modified_spi_threshold_list=numpy.arange(modified_spi_min,modified_spi_max+1,modified_spi_step).tolist()
# cti_threshold_list=numpy.arange(cti_min,cti_max+1,cti_step).tolist()

#terrain_index_log_thresholds_dict = { "specific_contributing_area": TerrainIndexThresholdListMaker(1.6,3.6,0.1,2)}#, "stream_power_index": spi_threshold_list,"modified_stream_power_index": modified_spi_threshold_list, "compound_topographic_index": cti_threshold_list}    
#terrain_index_list = ["CA","AS","AS2","CTI"] #["CA_log","AS_log","AS2_log","CTI_log"] 
# print(terrain_index_log_thresholds_dict)
terrain_index_stats_dict={}
terrain_index_errors_dict={}
errors_counter=0

for projection_WKID, huc12s_list in huc_projection_dict.items():
    print(f"WKID = {projection_WKID}" )
    print(f"HUC12s with {projection_WKID} are {huc12s_list}")
    terrain_index_stats_dict[projection_WKID]={}
    terrain_index_errors_dict[projection_WKID]={}
    for huc12 in huc12s_list:
        terrain_index_stats_dict[projection_WKID][huc12]={}
        terrain_index_errors_dict[projection_WKID][huc12]={}
        print(f"current HUC12 is {huc12}")

        location_dictionary,versioned_location_dictionary=LocationDictionaryFixer("2022",projection_WKID,huc12,"mnmx18","2","nrcs_gully_test6",run_info) 

        # path names if defined using Brian's defined flowpath program
        flow_direction_raster_location = versioned_location_dictionary["ag_fd"]
        flow_accumulation_raster_location = versioned_location_dictionary["ag_fa"]
        planform_curvature_raster_location = versioned_location_dictionary["ag_plan_crv"]
        grid_order_raster_location = versioned_location_dictionary["GordRaster"]
        
        shifted_flow_direction_raster_location=flow_direction_raster_location.replace('.tif', '_shifted.tif')
        shifted_flow_accumulation_raster_location=flow_accumulation_raster_location.replace('.tif', '_shifted.tif')
        shifted_planform_curvature_raster_location=planform_curvature_raster_location.replace('.tif', '_shifted.tif')
        shifted_grid_order_raster_location=grid_order_raster_location.replace('.tif', '_shifted.tif')

        if os.path.exists(shifted_flow_direction_raster_location):
            # print("shifted flow direction exists, changed to shifted version")
            flow_direction_raster_location=shifted_flow_direction_raster_location
            # print(flow_direction_raster_location)

        if os.path.exists(shifted_flow_accumulation_raster_location):
            # print("shifted flow accumulation exists, changed to shifted version")
            flow_accumulation_raster_location=shifted_flow_accumulation_raster_location
            # print(flow_accumulation_raster_location)

        if os.path.exists(shifted_planform_curvature_raster_location):
            # print("shifted planform curvature exists, changed to shifted version")
            planform_curvature_raster_location=shifted_planform_curvature_raster_location
            # print(planform_curvature_raster_location)

        if os.path.exists(shifted_grid_order_raster_location):
            # print("shifted grid order exists, changed to shifted version")
            grid_order_raster_location=shifted_grid_order_raster_location
            # print(grid_order_raster_location)

        #check if output directories and geodatabases exist and make them if they don't
        
        # #mean18 version
        # output_gdb_directory="M:\\DEP_nrcs_gully_test\\Gully_Outputs\\testing_area"

        #minmax18 version
        output_gdb_directory="M:\\DEP_nrcs_gully_test6\\Gully_Outputs\\testing_area"
        
        output_gdb_name = huc12 + ".gdb"
        full_output_gdb_path= os.path.join(output_gdb_directory,output_gdb_name)
        if not os.path.exists(output_gdb_directory):
            print("output gdb directory doesn't exist. creating it...")
            os.makedirs(output_gdb_directory)
            management.CreateFileGDB(output_gdb_directory,output_gdb_name)
        elif not os.path.exists(full_output_gdb_path):
            print("output gdb directory exists, but the gdb doesn't exist. creating it...")
            management.CreateFileGDB(output_gdb_directory,output_gdb_name)

        observed_ephemeral_gully_heads_gdb_name = "M:\\DEP\\Man_Data_Other\\Gully_data\\NRCS\\NRCS_Gulley_" + huc12 + ".gdb"
        joined_name="GullyHeads_" + huc12 + "_Joined" #"""need the name maker dictionary version of this"""   
        normal_name="GullyHeads_" + huc12

        with arcpy.EnvManager(workspace=observed_ephemeral_gully_heads_gdb_name):
            if arcpy.Exists(joined_name):
                print("_joined name exists")
                observed_ephemeral_gully_heads_raw= os.path.join(observed_ephemeral_gully_heads_gdb_name,joined_name)
            else:
                print("_joined name doesn't exist")
                observed_ephemeral_gully_heads_raw= os.path.join(observed_ephemeral_gully_heads_gdb_name,normal_name)
            
        # print(observed_ephemeral_gully_heads_raw)

        #set up names and variables
        observed_ephemeral_gully_heads = full_output_gdb_path+"\\GullyHeads_" + huc12 + "_Joined" + "_" + projection_WKID + "_reprojection"
        desc_fd = da.Describe(flow_direction_raster_location)
        cell_size = desc_fd['meanCellHeight'] 
        # snap_radius=cell_size*2**0.5
        snap_radius=float(6)
        snapped_gully_heads_raster_name=observed_ephemeral_gully_heads + "_snapped_raster"
        snapped_gully_heads_feature_name=observed_ephemeral_gully_heads+"_snapped"
        upstream_watersheds_raster_name=os.path.join(full_output_gdb_path,"upstream_watersheds")
        downstream_paths_raster_name=os.path.join(full_output_gdb_path,"downstream_flowpaths")
        valid_prediction_area_raster_name="valid_prediction_area_raster"
        longest_upstream_flowpaths_feature_name=observed_ephemeral_gully_heads+"_longest_upstream_fp_length"
        longest_downstream_flowpaths_feature_name=observed_ephemeral_gully_heads+"_longest_downstream_fp_length"

        #these are the locations of the outputs from the data preperation algorithm which calculates all the terrain indices
        stream_power_index_location = versioned_location_dictionary["stream_power_raster"]
        modified_stream_power_index_location = versioned_location_dictionary["modified_stream_power_raster"]
        specific_contributing_area_location = versioned_location_dictionary["specific_contributing_area_raster"]
        compound_topographic_index_location = versioned_location_dictionary["compound_topo_index_raster"]

        stream_power_index_location_log = stream_power_index_location.replace('.tif', '_log.tif')
        modified_stream_power_index_location_log = modified_stream_power_index_location.replace('.tif', '_log.tif')
        specific_contributing_area_location_log = specific_contributing_area_location.replace('.tif', '_log.tif')
        compound_topographic_index_location_log = compound_topographic_index_location.replace('.tif', '_log.tif')

        HUC12_terrain_index_raster_locations_dict = {"stream_power_index_location":stream_power_index_location,
        "modified_stream_power_index_location":modified_stream_power_index_location, 
        "specific_contributing_area_location":specific_contributing_area_location, 
        "compound_topographic_index_location":compound_topographic_index_location,
        "stream_power_index_location_log":stream_power_index_location_log, 
        "modified_stream_power_index_location_log":modified_stream_power_index_location_log,
        "specific_contributing_area_location_log":specific_contributing_area_location_log,
        "compound_topographic_index_location_log":compound_topographic_index_location_log}

        #setting the environment variables so they line up
        arcpy.env.snapRaster=flow_direction_raster_location
        arcpy.env.cellSize=flow_direction_raster_location

        
        GapDistancePreprocessing(observed_ephemeral_gully_heads_raw,observed_ephemeral_gully_heads,projection_WKID,flow_accumulation_raster_location,snap_radius,snapped_gully_heads_raster_name,snapped_gully_heads_feature_name,flow_direction_raster_location,upstream_watersheds_raster_name,downstream_paths_raster_name,full_output_gdb_path,valid_prediction_area_raster_name,longest_upstream_flowpaths_feature_name,longest_downstream_flowpaths_feature_name)
       
                
        for terrain_index_name,current_threshold_values_list in chunmei_terrain_index_log_thresholds.items():#terrain_index_log_thresholds_dict.items():
            print("Starting terrain index stuff")
            print(f"Current terrain index is {terrain_index_name}")
            print(f"Threshold value list is {current_threshold_values_list}.")# Beginning preprocessing...")
            terrain_index_stats_dict[projection_WKID][huc12][terrain_index_name]={}
            terrain_index_errors_dict[projection_WKID][huc12][terrain_index_name]={}
            # print(HUC12_terrain_index_raster_locations_dict)
            current_terrain_log_raster_variable_name=f"{terrain_index_name}_location_log"
            # print(current_terrain_log_raster_variable_name)
            current_terrain_log_raster=HUC12_terrain_index_raster_locations_dict[current_terrain_log_raster_variable_name]
            # print(current_terrain_log_raster)
            current_valid_terrain_index_log_raster_name=current_terrain_log_raster.replace('.tif', '_valid_flow_cells.tif')
            current_valid_terrain_index_log_raster= Con(os.path.join(full_output_gdb_path, valid_prediction_area_raster_name),current_terrain_log_raster)
            current_valid_terrain_index_log_raster.save(current_valid_terrain_index_log_raster_name)

            for current_threshold_value in current_threshold_values_list:
                gully_stats_error_list=[]
                print(f"Current threshold value is {current_threshold_value}")
                #threshold with con
                current_threshold_value_selection_string=f"VALUE < {current_threshold_value*100}"
                current_stats_results=[]
                try:
                    valid_terrain_index_log_raster_thresholded=SetNull(current_valid_terrain_index_log_raster_name,current_valid_terrain_index_log_raster_name,current_threshold_value_selection_string)
                
                    with da.SearchCursor(snapped_gully_heads_feature_name,"unique_flowpath_number") as unique_flowpath_number_list: #,where_clause="unique_flowpath_number = 20454"
                        for each_number in unique_flowpath_number_list:
                            
                            #take the unique flowpath ID and get the valid flow area raster for that gully head
                            current_number=each_number[0]
                            print(f"finding threshold point for gully head {current_number} and calculating")
                            current_valid_flow_area_raster_name=os.path.join(full_output_gdb_path, f"valid_flow_cells_{current_number}")
                            
                            try:    
                            # with arcpy.EnvManager(mask=current_valid_flow_area_raster_name):
                                current_valid_flow_index_extraction=ExtractByMask(valid_terrain_index_log_raster_thresholded,current_valid_flow_area_raster_name)
                                current_minimum=int(current_valid_flow_index_extraction.minimum)
                                current_minimum_raster=Con(in_conditional_raster= current_valid_flow_index_extraction,in_true_raster_or_constant=current_valid_flow_index_extraction, where_clause= f"VALUE = {current_minimum}")
                                current_minimum_point_name=os.path.join(full_output_gdb_path,f"current_TI_point")
                                current_minimum_point=arcpy.conversion.RasterToPoint(current_minimum_raster,current_minimum_point_name,"VALUE")
                                # current_maximum=int(current_valid_flow_index_extraction.maximum)
                                print(f"lowest {terrain_index_name} value for gully {current_number} without going under is {current_minimum}")
                                current_gully_head_point_name=os.path.join(full_output_gdb_path,f"current_gully_point")
                                current_gully_selection_string='"unique_flowpath_number" = ' + str(current_number)
                                current_gully_head_point=arcpy.analysis.Select(snapped_gully_heads_feature_name,current_gully_head_point_name,current_gully_selection_string)

                                arcpy.analysis.Near(current_gully_head_point,current_minimum_point)

                                with da.SearchCursor(current_gully_head_point,["unique_flowpath_number","NEAR_DIST"]) as cursor:
                                    for row in cursor:
                                        current_stats_results.append(row)
                                        # print(current_stats_results)
                                        #print(f"current unique ID is {row[0]}, closest object ID is {row[1]} and closest distance is {row[2]}")
                            
                            except:
                                errors_counter+=1
                                path_errors = []  
                                for error_info_item in sys.exc_info():#[0]
                                    # print(error_info_item.__class__.__text_signature__)
                                    # if type(error_info_item) is class("traceback"):
                                    #     error_info_item=traceback.format_list(traceback.extract_tb(error_info_item))
                                    path_errors.append(str(error_info_item))
                                gully_stats_error_list.append([current_number,path_errors])
                                print(f"Error on terrain index {terrain_index_name} at threshold value {current_threshold_value} and unique ID {current_number}: {path_errors}")
                
                except:
                    errors_counter+=1
                    path_errors = []
                    for error_info_item in sys.exc_info():#[0]
                        path_errors.append(str(error_info_item))
                    gully_stats_error_list.append([path_errors])
                    print(f"Error on terrain index {terrain_index_name} thresholding at threshold value {current_threshold_value}: {path_errors}")

                terrain_index_stats_dict[projection_WKID][huc12][terrain_index_name][current_threshold_value]=current_stats_results
                terrain_index_errors_dict[projection_WKID][huc12][terrain_index_name][current_threshold_value]=gully_stats_error_list        
                        # print(f"highest {terrain_index_name} value for gully {current_number} is {current_maximum}")
                        # current_upstream_watershed_raster_name=os.path.join(full_output_gdb, f"upstream_watershed_{current_number}")
                        # current_downstream_flowpath_raster_name=os.path.join(full_output_gdb, f"downstream_flowpath_{current_number}")
                        
print(terrain_index_stats_dict)
print(terrain_index_errors_dict)
            
# create json object from dictionary
json_stats = json.dumps(terrain_index_stats_dict)
json_stats_file_path=os.path.join(output_gdb_directory, "testing_minmax_stats_6m_snap_results.json")#"102400060304_sca_10point_stats_results.json")
json_file = open(json_stats_file_path,"w")
json_file.write(json_stats)
json_file.close()


json_stats_errors = json.dumps(terrain_index_errors_dict)
json_stats_errors_file_path=os.path.join(output_gdb_directory, "testing_minmax_stats_6m_snap_errors.json")#"102400060304_sca_10point_stats_errors.json")
json_errors_file =open(json_stats_errors_file_path,"w")
json_errors_file.write(json_stats_errors)
json_errors_file.close()

end_time=datetime.datetime.strftime(datetime.datetime.now(), '%Y_%m_%d_%H_%M_%S')
print(f"started run at {start_time}")
print(f"ended run at {end_time}")

#make json path

# open file for writing, "w" 

# write json object to file

# close files




                


        # for terrain_index_type, threshold_values_list in terrain_index_log_thresholds_dict.items():
        #     for current_threshold_value in threshold_values_list:


        # # workbook=xlsxwriter.Workbook(out_folder+'/result.xlsx')
        # # worksheet=workbook.add_worksheet('result')
        # # worksheet.write_row('A1',['Threshold','max_value'])
        # threshold_list=[]
        # mean_list=[]

        # for TI in range(TI_start,TI_end+1,step):
        #     print TI

        #     gdb_name="result"+str(TI)+'.gdb'
        #     out_gdb = out_folder + '/' + gdb_name

        #     if os.path.exists(out_gdb):
        #         print gdb_name+"gdb exists"
        #         shutil.rmtree(out_gdb)

        #     arcpy.CreateFileGDB_management(out_folder,gdb_name)

        #     # Process: Snap Pour Point
        #     Gullyheadacc = out_gdb+'/Gullyhead'+str(TI) # provide a default value if unspecified

        #     StreamT_stream=out_gdb+'/StreamT_stream'+str(TI)
        #     stream_orderf=out_gdb+"/stream_orderf"+str(TI)
        #     Endpoint=out_gdb+"/Endpoint"+str(TI)
        #     Output_direction_raster=out_gdb+"/Output_direction_raster"+str(TI)
        #     outWatershed=out_gdb+"/watershed"+str(TI)
        #     # Local variables:
        #     outSnapPour=SnapPourPoint(EG_Heads, FlowAcc, "2", "OBJECTID")
        #     # outSnapPour.save(out_gdb+'/outSnapPour'+str(TI))
        #     # Process: Raster to Point (2)
        #     tempEnvironment0 = arcpy.env.snapRaster
        #     arcpy.env.snapRaster = FlowAcc
        #     tempEnvironment1 = arcpy.env.cellSize
        #     arcpy.env.cellSize = FlowAcc
        #     arcpy.RasterToPoint_conversion(outSnapPour, Gullyheadacc)
        #     arcpy.env.snapRaster = tempEnvironment0
        #     arcpy.env.cellSize = tempEnvironment1

        #     # outCon=Con(stream, "1", stream,  "Value >= 234")
        #     # arcpy.gp.Con_sa(TOPO_INDEX, Input_true_raster_or_constant_value__2_, stream, "", "\"Value\" >=" + str(TI))
        #     outGreater=GreaterThanEqual(TOPO_INDEX, TI)
        #     outSetNull = SetNull(outGreater, outGreater, "VALUE=0")
        #     outSetNull.save(out_gdb+'/outSetNULL'+str(TI))
        #     # Process: Stream to Feature (3)
        #     StreamToFeature(outSetNull, FlowDir, StreamT_stream, "SIMPLIFY")

        #     # Process: Dissolve (3)
        #     arcpy.Dissolve_management(StreamT_stream, stream_orderf, "", "", "MULTI_PART", "DISSOLVE_LINES")

        #     # Process: Feature Vertices To Points (3)
        #     arcpy.FeatureVerticesToPoints_management(stream_orderf, Endpoint, "START")

        #     # Process: Near (2)
        #     tempEnvironment0 = arcpy.env.cellSize
        #     arcpy.env.cellSize = FlowAcc
        #     arcpy.Near_analysis(Gullyheadacc, Endpoint, "300 Meters", "NO_LOCATION", "NO_ANGLE")#用原有沟头？？
        #     arcpy.env.cellSize = tempEnvironment0

        #     # Process: Watershed (2)
        #     tempEnvironment0 = arcpy.env.snapRaster
        #     arcpy.env.snapRaster = FlowAcc
        #     outWatershed1=Watershed(FlowDir, outSnapPour, "VALUE")
        #     outWatershed1.save(outWatershed)
        #     arcpy.env.snapRaster = tempEnvironment0

        #     # Process: Euclidean Distance (2)
        #     tempEnvironment0 = arcpy.env.snapRaster
        #     arcpy.env.snapRaster = FlowAcc
        #     tempEnvironment1 = arcpy.env.extent
        #     arcpy.env.extent = EG_Heads
        #     outEucDistance=EucDistance(Gullyheadacc,  "300", TOPO_INDEX, Output_direction_raster)#原有沟头？？
        #     outEucDistance.save(out_gdb+'/EucDistance'+str(TI))
        #     arcpy.env.snapRaster = tempEnvironment0
        #     arcpy.env.extent = tempEnvironment1

        #     # Process: Extract by Mask
        #     outExtractByMask=ExtractByMask(outEucDistance, outSetNull)
        #     outExtractByMask.save(out_gdb+'/ExtractByMask'+str(TI))
        #     # Process: Zonal Statistics (2)
        #     outZonal=ZonalStatistics(outWatershed, "VALUE", outExtractByMask, "MAXIMUM", "DATA")

        #     # Process: Extract Multi Values to Points (2)
        #     ExtractMultiValuesToPoints(Gullyheadacc, outZonal, "NONE")


        #     arcpy.AddField_management(Gullyheadacc,'Dis_merge','FLOAT')
        #     near_dist=[]
        #     zonalst_list=[]
        #     merge_list=[]

        #     with arcpy.da.UpdateCursor(Gullyheadacc,["NEAR_DIST","ZonalSt_Wate1","Dis_merge"]) as cursor:
        #         for row in cursor:
        #             if row[1]==None:
        #                 merge_list.append(row[0])
        #                 row[2] = row[0]
        #             else:
        #                 merge_list.append(row[1])
        #                 row[2] = row[1]
        #             cursor.updateRow((row))

        #     mean_list.append(np.mean(merge_list))
        #     threshold_list.append(TI)

        # # worksheet.write_column('A2',threshold_list)
        # # worksheet.write_column('B2',mean_list)
        # # workbook.close()

