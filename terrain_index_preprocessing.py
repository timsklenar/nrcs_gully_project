# -*- coding: utf-8 -*-
#输入全部流域文件夹，输出结果全部保存输入文件夹下对应流域的数据库中。Enter all watershed folders, and save the output results in the database of the corresponding watershed under the input folder. (That means this indexes calculation script can loop from watershed to watershed)
#stream_power_name、modified_stream_power_name、compound_topo_index_name、specific_contributing_area_name分别为输出四个地形指数文件名 The names for the output terrain indexes
#输入文件夹对应流域数据库中要确保存在Flowacc(汇流累积）、FlowDir(流向）、Curvatu_Plan（平面曲率）Make sure you have Flowacc FlowDir Curvatu_Plan in your geodatabase for each watershed  (calculated using Arctoolbox)
#对应流域文件夹中要确保存在sd8.tif(坡度） Make sure you have slope gradient (calculated using TAUDEM) in your foder for each watershed 
import arcpy
from arcpy.sa import *
import os
from os.path import join as opj
netid = os.getlogin()
# connect to netid
import sys
if netid == 'bkgelder':
    sys.path.append(opj('C:\\Users', netid, 'Box\\Data_Sharing\\Scripts\\basics'))
else:
    sys.path.append(opj('C:\\Users', netid, 'Box\\Scripts\\basics'))
import dem_functions2 as df
import platform
import datetime

arcpy.CheckOutExtension("spatial")
arcpy.env.overwriteOutput = True

##in_folder='F:/wang_teacher/GRID_ORDER/wsds2'
##
##in_slp_d8_name = 'sd8.tif'
##in_FlowDir_name = 'FlowDir'
##in_FlowAcc_name = 'FlowAcc'
##in_cur_plan_name= 'Curvatu_Plan'

# stream_power_name = "stream_power_LOG"
# modified_stream_power_name = "modified_stream_power_LOG"
# compound_topo_index_name = "compound_topo_index_LOG"
# specific_contributing_area_name = "specific_contributing_area_LOG"
# flowlength_up_name = "Flowlength_up1"
#EG_heads_vector_name = "EG_heads_vector"
#out_folder="F:/wang_teacher/proc_result2"
#gdb_name='INDEX.gdb'


# def read_raster(raster):
#     inras = arcpy.Raster(raster)
#     lowerLeft = arcpy.Point(inras.extent.XMin, inras.extent.YMin)
#
#     prj = arcpy.Describe(raster).spatialReference
#
#     ras_array=arcpy.RasterToNumPyArray(raster, nodata_to_value=-9999)
#
#     cell_size = float(str(arcpy.GetRasterProperties_management(raster, "CELLSIZEX")))
#
#     return ras_array, lowerLeft, cell_size, prj

#对栅格文件求log10后取整 have log10 and int calculation (actually also have times 100)
def log_int(in_tif):
    a=Log10(in_tif)
    b=Times(a, 100)
    c=Int(b)
    return c


#check the bottom left extent values to see if there is a difference between 2 rasters, return a list with float values for the differences. 
# values returned are first minus second and can be used to shift the second to match the first 
def extent_checker(first_input_raster, second_input_raster):
    first_raster_left_extent_result = arcpy.GetRasterProperties_management(first_input_raster,"LEFT")
    first_raster_left_extent = first_raster_left_extent_result.getOutput(0)
    
    first_raster_bottom_extent_result = arcpy.GetRasterProperties_management(first_input_raster,"BOTTOM")
    first_raster_bottom_extent = first_raster_bottom_extent_result.getOutput(0)
    
    second_raster_left_extent_result = arcpy.GetRasterProperties_management(second_input_raster,"LEFT")
    second_raster_left_extent = second_raster_left_extent_result.getOutput(0)
    
    second_raster_bottom_extent_result = arcpy.GetRasterProperties_management(second_input_raster,"BOTTOM")
    second_raster_bottom_extent = second_raster_bottom_extent_result.getOutput(0)

    left_shift_value = float(first_raster_left_extent) - float(second_raster_left_extent)
    bottom_shift_value = float(first_raster_bottom_extent) - float(second_raster_bottom_extent)
    return [left_shift_value,bottom_shift_value]


# def caculate_index(out_folder,sd8,Curvatu_Plan,FlowAcc,Snap_EG_HEAD,Flowdir,gdb_name):
def calculate_index(sd8,Curvatu_Plan,FlowAcc,Flowdir, stream_power_raster, modified_stream_power_raster, specific_contributing_area_raster, compound_topo_index_raster, cell_size, stream_power_raster_log, modified_stream_power_raster_log, specific_contributing_area_raster_log, compound_topo_index_raster_log):
    ##创建数据库 create an database
    #arcpy.CreateFileGDB_management(out_folder, gdb_name)
    # gdb_dir=out_folder+'/'+gdb_name

    # stream_power_dir = gdb_dir + "/" +stream_power_name
    # modified_stream_power_dir = gdb_dir + "/" + modified_stream_power_name
    # compound_topo_index_dir = gdb_dir + "/" + compound_topo_index_name
    # specific_contributing_area_dir = gdb_dir + "/" + specific_contributing_area_name
    # #EG_heads_vector_dir = gdb_dir + "/" + EG_heads_vector_name
    # flowlength_up_dir = gdb_dir + "/" + flowlength_up_name

    # Process: Flow Length and vector heads
    # flowlength=FlowLength(Flowdir, "UPSTREAM", "")
    # flowlength.save(flowlength_up_dir)
    #arcpy.RasterToPoint_conversion(Snap_EG_HEAD, EG_heads_vector_dir, "Value")

    for i in [stream_power_raster, modified_stream_power_raster, specific_contributing_area_raster, compound_topo_index_raster]:
        dir_to_make = os.path.dirname(i)
        if not os.path.isdir(dir_to_make):
            os.makedirs(dir_to_make)

    #分别计算坡度、汇流累积量、曲率 Preparing from sd8, flowacc and curvature
    # make sure 0 is not an error when doing log calculations
    slope_percent=Plus(sd8, 0.0001)
    slope_fraction=Divide (slope_percent, 100)

    accumulated_flow_area=Times(FlowAcc , cell_size**2)#9)
    avg_cell_distance = (cell_size + cell_size * 2**0.5)/2
    specific_contributing_area=Divide(accumulated_flow_area, avg_cell_distance)#3.621)
    # adjust to make sure 0 is not an error when doing log calculations
    specific_contributing_area_adjusted=Plus(specific_contributing_area, 0.0001)
    specific_contributing_area_adjusted.save(specific_contributing_area_raster)
    specific_contributing_area_log10=log_int(specific_contributing_area_adjusted)
    specific_contributing_area_log10.save(specific_contributing_area_raster_log)

    negative_planform_curvature=Times(Curvatu_Plan,-1)
    negative_planform_curvature_adjusted=Divide (negative_planform_curvature,100)

    #arcpy.gp.RasterCalculator_sa(Int(Log10((FlowAcc * 9 / 3.621 + 0.0001) * (dem + 0.0001) / 100) * 100), SA)
    stream_power_index=Times(slope_fraction, specific_contributing_area_adjusted)
    stream_power_index.save(stream_power_raster)
    stream_power_index_log10=log_int(stream_power_index)
    stream_power_index_log10.save(stream_power_raster_log)#dir)

    #arcpy.gp.RasterCalculator_sa("Int(Log10((\"%FlowAcc%\"  * 9 / 3.621 + 0.0001) * (\"%ep3m070600040608sd8.tif%\" + 0.0001)/ 100 *(\"%ep3m070600040608sd8.tif%\" + 0.0001)/ 100) * 100) ", AS2)
    stream_power_index_variant=Times(stream_power_index,slope_fraction)
    stream_power_index_variant.save(modified_stream_power_raster)
    stream_power_index_variant_log10=log_int(stream_power_index_variant)
    stream_power_index_variant_log10.save(modified_stream_power_raster_log)#dir)

    #arcpy.gp.RasterCalculator_sa("Int(Log10((\"%FlowAcc%\" * 9  + 0.0001)* (\"%ep3m070600040608sd8.tif%\" + 0.0001) / 100 * ( - 1) * \"%Curvatu_Plan%\" / 100) * 100)", CTI)
    compound_topographic_index=Times(stream_power_index,negative_planform_curvature_adjusted)
    compound_topographic_index.save(compound_topo_index_raster)
    compound_topographic_index_log10=log_int(compound_topographic_index)
    compound_topographic_index_log10.save(compound_topo_index_raster_log)#dir)

    #arcpy.gp.RasterCalculator_sa("Int(Log10(\"%FlowAcc%\"  * 9 / 3.621 + 0.0001) * 100)", CA)
    """this doesn't do anything, it's just plan curve times specific area when specific area is what we wanted saved out as CA. It's a confusion because of the 
    poor naming scheme and is meaningless"""
    # CA=Times(A,C)
    # specific_contributing_area_log10=log_int(CA)
    # specific_contributing_area_log10.save(specific_contributing_area_raster)#dir)

#按流域文件夹进行处理 loop in watersheds
##def deal(in_folder):
##    # wsds=os.listdir(in_folder)
##    # for wsd in wsds: 


huc_projection_dict = {"26914": ["102400060304","102702050101","102702070207"], "26915": ["070801030408","070801050302","070802050807", "071000040910", "102300030509","102300031003","102300031403"]} 
#"102702060102" data is missing in the file structure for mnmx but exists for mean as of 7/24/23

if __name__ == '__main__':
##    deal(in_folder)
    for projection_WKID, huc12s_list in huc_projection_dict.items():
        print(f"WKID = {projection_WKID}" )
        print(f"HUC12s with {projection_WKID} are {huc12s_list}")
        huc_list = huc12s_list
        
        for huc12 in huc_list:
            print(huc12)
            
            ## mean18 parameters
            # if len(sys.argv) == 1:
            #     cleanup = False
            #     parameters = ["C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/pythonw.exe",
            # "O:/DEP/Scripts/basics/cmd_FlowPath_v9.py",
            # platform.node(),
            # "2021",
            # "070801050302",
            # "26914", #"26915", #
            # "mean18",
            # "2"]

            # minmax18 parameters
            if len(sys.argv) == 1:
                cleanup = False
                parameters = ["C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/pythonw.exe",
            "O:/DEP/Scripts/basics/cmd_FlowPath_v9.py",
            platform.node(),
            "2021",
            "070801050302",
            projection_WKID,
            "mnmx18",
            "2"]

                for i in parameters[2:]:
                    sys.argv.append(i)

            else:
                cleanup = True
    ####            node = sys.argv[1]
            ACPFyear = sys.argv[2]
    ####            huc12 = sys.argv[3]
            srOutCode = sys.argv[4]
            interpType = sys.argv[5]
            cellSize = int(sys.argv[6])

            nowYmd = datetime.datetime.strftime(datetime.datetime.now(), '%Y_%m_%d_%H_%M_%S')

            locDict = df.loadVariablesDict(platform.node(), ACPFyear, huc12, srOutCode, interpType, cellSize, nowYmd)

            # #mean18 version
            # locDictVer = df.loadVariablesDict(platform.node(), ACPFyear, huc12, srOutCode, interpType, cellSize, nowYmd, 'nrcs_gully_test')

            #minmax18 version
            locDictVer = df.loadVariablesDict(platform.node(), ACPFyear, huc12, srOutCode, interpType, cellSize, nowYmd, 'nrcs_gully_test6')

            # Brian's defined flowpath program creates this stuff
            slp_d8 = locDictVer["tau_slope"].replace('\\\\EL3354-02\\O$', 'M:')
            FlowDir = locDictVer["ag_fd"].replace('\\\\EL3354-02\\O$', 'M:')
            FlowAcc = locDictVer["ag_fa"].replace('\\\\EL3354-02\\O$', 'M:')
            Curvatu_Plan = locDictVer["ag_plan_crv"].replace('\\\\EL3354-02\\O$', 'M:')
            grid_order_raster= locDictVer["GordRaster"].replace('\\\\EL3354-02\\O$', 'M:')
            
            # These are new things to create to calculate indices        
            stream_power_location = locDictVer["stream_power_raster"].replace('\\\\EL3354-02\\O$', 'M:')
            modified_stream_power_location = locDictVer["modified_stream_power_raster"].replace('\\\\EL3354-02\\O$', 'M:')
            specific_contributing_area_location = locDictVer["specific_contributing_area_raster"].replace('\\\\EL3354-02\\O$', 'M:')
            compound_topo_index_location = locDictVer["compound_topo_index_raster"].replace('\\\\EL3354-02\\O$', 'M:')

            stream_power_location_log = stream_power_location.replace('.tif', '_log.tif')
            modified_stream_power_location_log = modified_stream_power_location.replace('.tif', '_log.tif')
            specific_contributing_area_location_log = specific_contributing_area_location.replace('.tif', '_log.tif')
            compound_topo_index_location_log = compound_topo_index_location.replace('.tif', '_log.tif')
            
            #setting up pararmeters to make sure the input files align to the extent of the d8 slope files
            desc_fd = arcpy.da.Describe(slp_d8)
            cell_size = desc_fd['meanCellHeight']                    
            arcpy.env.snapRaster=slp_d8
            arcpy.env.cellSize=slp_d8

            #checking the extents and aligning if neccessary
            flow_direction_shift_values=extent_checker(slp_d8,FlowDir)
            flow_accumulation_shift_values=extent_checker(slp_d8,FlowAcc)
            grid_order_shift_values=extent_checker(slp_d8,grid_order_raster)
            planform_curvature_shift_values=extent_checker(slp_d8,Curvatu_Plan)

            if flow_direction_shift_values[0] != 0 or flow_direction_shift_values[1] != 0:
                shifted_flow_direction_location=FlowDir.replace('.tif', '_shifted.tif')
                arcpy.management.Shift(FlowDir, shifted_flow_direction_location, flow_direction_shift_values[0],flow_direction_shift_values[1], slp_d8)
                FlowDir=shifted_flow_direction_location

            if flow_accumulation_shift_values[0] != 0 or flow_accumulation_shift_values[1] != 0:
                shifted_flow_accumulation_location=FlowAcc.replace('.tif', '_shifted.tif')
                arcpy.management.Shift(FlowAcc, shifted_flow_accumulation_location, flow_accumulation_shift_values[0],flow_accumulation_shift_values[1], slp_d8)
                FlowAcc=shifted_flow_accumulation_location

            if grid_order_shift_values[0] != 0 or grid_order_shift_values[1] != 0:
                shifted_grid_order_location=grid_order_raster.replace('.tif', '_shifted.tif')
                arcpy.management.Shift(grid_order_raster, shifted_grid_order_location, grid_order_shift_values[0],grid_order_shift_values[1], slp_d8)
                grid_order_raster=shifted_grid_order_location

            if planform_curvature_shift_values[0] != 0 or planform_curvature_shift_values[1] != 0:
                shifted_planform_curvature_location=Curvatu_Plan.replace('.tif', '_shifted.tif')
                arcpy.management.Shift(Curvatu_Plan, shifted_planform_curvature_location, planform_curvature_shift_values[0],planform_curvature_shift_values[1], slp_d8)
                Curvatu_Plan=shifted_planform_curvature_location

            calculate_index(slp_d8, Curvatu_Plan, FlowAcc, FlowDir, stream_power_location, modified_stream_power_location, specific_contributing_area_location, compound_topo_index_location, cell_size, stream_power_location_log, modified_stream_power_location_log, specific_contributing_area_location_log, compound_topo_index_location_log)

            # this should get you the HUC12 point feature class of gully heads for the specified HUC12
            # defined_ends = locDict['nrcs_gullies'].replace('\\\\EL3354-02\\O$', 'M:')
