#%%
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


run_info=sys.argv

for each in sys.argv:
    if "jupyter" in each:
        run_info=[sys.argv[0]]


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



#%%
start_time=datetime.datetime.now()
# start_time_string=datetime.datetime.strftime(start_time, '%Y_%m_%d_%H_%M_%S')
print(f"{start_time}")

#data parameters for the run
year="2022"
resolution="2"
dem_method="mnmx18"
versioning_string="nrcs_gully_test6"

#dictionary of the HUCs to run
# huc_projection_dict = {"26914": ["102400060304","102702050101","102702070207"], "26915": ["070801030408","070801050302","070802050807", "071000040910", "102300030509","102300031003","102300031403"]} 
huc_projection_dict = {"26914": ["102400060304"]}

for projection_WKID, huc12s_list in huc_projection_dict.items():
    for huc12 in huc12s_list:
        location_dictionary,versioned_location_dictionary=LocationDictionaryFixer(year,projection_WKID,huc12,dem_method,resolution,versioning_string,run_info)
        # arcpy.env.parallelProcessingFactor = "75%"
        c_dem=location_dictionary["cElevFile"]
        #%%
        c_dem_basename=os.path.basename(os.path.splitext(c_dem)[0])
        print("starting fill")
        test_fill=Fill(c_dem)
        print("fill save")
        current_gdb_path=os.path.join(versioned_location_dictionary["depBase"],fr"Gully_Outputs\testing_area\{huc12}.gdb")
        filled_c_dem_path=os.path.join(current_gdb_path,f"filled_{c_dem_basename}")
        #%%
        test_fill.save(filled_c_dem_path)
        print("flow direction")
        fd=FlowDirection(test_fill)
        print("flow acc")
        fa=FlowAccumulation(fd)
        print("fa save")
        filled_c_fa_path=os.path.join(current_gdb_path,f"fa_filled_{c_dem_basename}")
        fa.save(filled_c_fa_path)
        unfilled_fa=versioned_location_dictionary["ag_fa"]
        print("con statement")
        headwaters_gully_mask_raster=Con(unfilled_fa>=fa,1)
        print("con save")
        valid_headwaters_raster_path=os.path.join(current_gdb_path,f"valid_headwaters_{c_dem_basename}")
        headwaters_gully_mask_raster.save(valid_headwaters_raster_path)

        p_dem=location_dictionary["pElevFile"]
        p_dem_basename=os.path.basename(os.path.splitext(p_dem)[0])
        print("starting fill")
        test_fill=Fill(p_dem)
        print("fill save")
        current_gdb_path=os.path.join(versioned_location_dictionary["depBase"],fr"Gully_Outputs\testing_area\{huc12}.gdb")
        filled_p_dem_path=os.path.join(current_gdb_path,f"filled_{p_dem_basename}")
        test_fill.save(filled_p_dem_path)
        print("flow direction")
        fd=FlowDirection(test_fill)
        print("flow acc")
        fa=FlowAccumulation(fd)
        print("fa save")
        filled_p_fa_path=os.path.join(current_gdb_path,f"fa_filled_{p_dem_basename}")
        fa.save(filled_p_fa_path)
        unfilled_fa=versioned_location_dictionary["ag_fa"]
        print("con statement")
        headwaters_gully_mask_raster=Con(unfilled_fa>=fa,1)
        print("con save")
        valid_headwaters_raster_path=os.path.join(current_gdb_path,f"valid_headwaters_{p_dem_basename}")
        headwaters_gully_mask_raster.save(valid_headwaters_raster_path)

        # p_dem=r"M:\DEP\LiDAR_Current\elev_PLib_mnmx18\07080205\ep2m070802050807.tif"
        # test_fill=Fill(p_dem)
        # test_fill.save(r"M:\DEP_nrcs_gully_test6\Gully_Outputs\testing_area\070802050807.gdb\filled_p_dem")
        # fd=FlowDirection(test_fill)
        # fa=FlowAccumulation(fd)
        # fa.save(r"M:\DEP_nrcs_gully_test6\Gully_Outputs\testing_area\070802050807.gdb\filled_p_fa_2m_070802050807")
        # unfilled_fa=r"M:\DEP_nrcs_gully_test6\Gully_Outputs\AG_FA_mnmx18\07080205\fa_2m_070802050807.tif"
        # headwaters_gully_mask_raster=Con(unfilled_fa>=fa,1)
        # headwaters_gully_mask_raster.save(r"M:\DEP_nrcs_gully_test6\Gully_Outputs\testing_area\070802050807.gdb\valid_headwaters_gullies_p_dem")

end_time=datetime.datetime.now()
# end_time_string=datetime.datetime.strftime(end_time, '%Y_%m_%d_%H_%M_%S')
print(f"{end_time}")
# %%
run_time=end_time-start_time
# run_time_string=datetime.datetime.strftime(run_time, '%Y_%m_%d_%H_%M_%S')
print(f"{run_time}")
# %%
