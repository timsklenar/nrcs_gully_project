#%%
import arcpy
from arcpy import management
from arcpy import da
from arcpy.sa import *
import os
import platform
import datetime
import pandas as pd
import numpy as np
import os
import shutil
import sys
import json
from itertools import repeat
given_value ='Hello! '
new_list=[]
new_list.extend(repeat(given_value,5))
#%%
#found these handy functions to take a featureclass, input fields list and sql query to create a pandas dataframe 
def arcgis_table_to_df(in_fc, input_fields=None, query=""):
    """Function will convert an arcgis table into a pandas dataframe with an object ID index, and the selected
    input fields using an arcpy.da.SearchCursor.
    :param - in_fc - input feature class or table to convert
    :param - input_fields - fields to input to a da search cursor for retrieval
    :param - query - sql query to grab appropriate values
    :returns - pandas.DataFrame"""
    OIDFieldName = arcpy.Describe(in_fc).OIDFieldName
    if input_fields:
        final_fields = [OIDFieldName] + input_fields
    else:
        final_fields = [field.name for field in arcpy.ListFields(in_fc)]
    data = [row for row in arcpy.da.SearchCursor(in_fc,final_fields,where_clause=query)]
    fc_dataframe = pd.DataFrame(data,columns=final_fields)
    fc_dataframe = fc_dataframe.set_index(OIDFieldName,drop=True)
    return fc_dataframe

def arcgis_table_to_dataframe(in_fc, input_fields, query="", skip_nulls=False, null_values=None):
    """Function will convert an arcgis table into a pandas dataframe with an object ID index, and the selected
    input fields. Uses TableToNumPyArray to get initial data.
    :param - in_fc - input feature class or table to convert
    :param - input_fields - fields to input into a da numpy converter function
    :param - query - sql like query to filter out records returned
    :param - skip_nulls - skip rows with null values
    :param - null_values - values to replace null values with.
    :returns - pandas dataframe"""
    OIDFieldName = arcpy.Describe(in_fc).OIDFieldName
    if input_fields:
        final_fields = [OIDFieldName] + input_fields
    else:
        final_fields = [field.name for field in arcpy.ListFields(in_fc)]
    np_array = arcpy.da.TableToNumPyArray(in_fc, final_fields, query, skip_nulls, null_values)
    object_id_index = np_array[OIDFieldName]
    fc_dataframe = pd.DataFrame(np_array, index=object_id_index, columns=input_fields)
    return fc_dataframe

#funtion to create a pandas dataframe and a calculate the NSE for longest upstream flowpath for each HUC12/terrain index type/terrain index threshold combination
def GullyStatsDataFrameConverter(json_data,WKID, HUC12, terrain_index_name, terrain_index_threshold_string):
    sliced_stats_dataframe=pd.DataFrame(json_data[WKID][HUC12][terrain_index_name][terrain_index_threshold_string], columns = ['unique_flowpath_number', 'gap_distance'])  
    # print(sliced_stats_dataframe)
    longest_upstream_paths_feature=f"M:\\DEP_nrcs_gully_test6\\Gully_Outputs\\testing_area\\{HUC12}.gdb\\GullyHeads_{HUC12}_Joined_{WKID}_reprojection_longest_upstream_fp_length"
    upstream_fields_list=[field.name for field in arcpy.ListFields(longest_upstream_paths_feature)]
    if "grid_code" in upstream_fields_list:
        longest_upstream_paths_df=arcgis_table_to_df(longest_upstream_paths_feature,["grid_code","RASTERVALU"])
        longest_upstream_paths_df = longest_upstream_paths_df.rename(columns={"grid_code": "unique_flowpath_number","RASTERVALU":"longest_upstream_path"})
    else:
        longest_upstream_paths_df=arcgis_table_to_df(longest_upstream_paths_feature,["unique_flowpath_number","RASTERVALU"])
        longest_upstream_paths_df = longest_upstream_paths_df.rename(columns={"RASTERVALU":"longest_upstream_path"})

    
    longest_downstream_paths_feature=f"M:\\DEP_nrcs_gully_test6\\Gully_Outputs\\testing_area\\{HUC12}.gdb\\GullyHeads_{HUC12}_Joined_{WKID}_reprojection_longest_downstream_fp_length"
    downstream_fields_list=[field.name for field in arcpy.ListFields(longest_downstream_paths_feature)]
    if "grid_code" in downstream_fields_list:
        longest_downstream_paths_df=arcgis_table_to_df(longest_downstream_paths_feature,["grid_code","RASTERVALU"])
        longest_downstream_paths_df = longest_downstream_paths_df.rename(columns={"grid_code":"unique_flowpath_number","RASTERVALU":"longest_downstream_path"})
    else:    
        longest_downstream_paths_df=arcgis_table_to_df(longest_downstream_paths_feature,["unique_flowpath_number","RASTERVALU"])
        longest_downstream_paths_df = longest_downstream_paths_df.rename(columns={"RASTERVALU":"longest_downstream_path"})
    
    # print(longest_upstream_paths_df)
    # print(longest_downstream_paths_df)
    
    full_stats_df=pd.merge(sliced_stats_dataframe,longest_upstream_paths_df,on="unique_flowpath_number",how="outer")
    
    full_stats_df=pd.merge(full_stats_df,longest_downstream_paths_df,on="unique_flowpath_number",how="outer")
    # print(full_stats_df)
    full_stats_df["failure_value"]=np.where(full_stats_df["longest_downstream_path"] < full_stats_df["longest_upstream_path"],full_stats_df["longest_upstream_path"],full_stats_df["longest_downstream_path"])
    # print(full_stats_df)
    upstream_flowpath_mean=np.mean(full_stats_df["longest_upstream_path"])
    upstream_flowpath_mean_column=[]
    upstream_flowpath_mean_column.extend(repeat(upstream_flowpath_mean,len(full_stats_df)))
    print(f"upstream_flowpath_mean is {upstream_flowpath_mean}")
    full_stats_df["HUC12_upstream_flowpath_mean"]=upstream_flowpath_mean_column
    full_stats_df["population_variance"]=(full_stats_df["longest_upstream_path"]-upstream_flowpath_mean)**2
    # print(full_stats_df)
    full_stats_df["error_variance"]=np.where(full_stats_df["gap_distance"].isnull(), full_stats_df["failure_value"]**2,full_stats_df["gap_distance"]**2)
    
    NSE_denominator=np.sum(full_stats_df["population_variance"])
    NSE_nummerator=np.sum(full_stats_df["error_variance"])
    NSE_value=float(1)-(NSE_nummerator/NSE_denominator)
    # print(full_stats_df.dtypes)
    print(full_stats_df)
    print(NSE_nummerator)
    print(NSE_denominator)
    print(NSE_value)
    return full_stats_df,NSE_value
#%%
#setting up the dictionaries used to iterate throught the data (projection WKIDs as keys with HUC12 lists as values, second is terrain index names as keys and threshold lists as values)
huc_projection_dict = {"26914": ["102400060304","102702050101","102702070207"], "26915": ["070801030408","070801050302","070802050807", "071000040910", "102300030509","102300031003","102300031403"]} 
# {"26915": ["070801050302","071000040910", "102300031003","102300031403"],"26914": ["102400060304","102702070207","102702050101","102702060102"]} 
chunmei_terrain_index_log_thresholds={"specific_contributing_area":[2.6],"stream_power_index": [1.5],"modified_stream_power_index": [0.25], "compound_topographic_index":[1.2]}
json_file_path="M:\\DEP_nrcs_gully_test6\\Gully_Outputs\\testing_area\\testing_minmax_stats_6m_snap_results.json"

#loading in the stats jason
# huc_projection_dict = {"26915": ["070801050302"]}
# chunmei_terrain_index_log_thresholds={"stream_power_index": [1.5]}
# json_file_path="M:\\DEP_nrcs_gully_test\\Gully_Outputs\\testing_area\\full_stats_results_070801050302.json"

# huc_projection_dict = {"26914": ["102400060304"]}
# chunmei_terrain_index_log_thresholds={"specific_contributing_area":[2.6],"stream_power_index": [1.5]}
# json_file_path="M:\\DEP_nrcs_gully_test\\Gully_Outputs\\testing_area\\full_stats_results_102400060304_sca_spi.json"

# huc_projection_dict = {"26914": ["102400060304"]}
# specific_area_min=1.6
# specific_area_max=3.6
# specific_area_step=0.1
# specific_area_samples_number=int(((specific_area_max-specific_area_min)/specific_area_step)+1)
# specific_area_threshold_list=np.linspace(specific_area_min,specific_area_max,num=specific_area_samples_number).tolist()
# print(specific_area_threshold_list)
# specific_area_threshold_list=list(np.around(np.array(specific_area_threshold_list),2))
# print(specific_area_threshold_list)
# terrain_index_log_thresholds_dict = { "specific_contributing_area": specific_area_threshold_list}
# print(terrain_index_log_thresholds_dict)
# json_file_path="M:\\DEP_nrcs_gully_test6\\Gully_Outputs\\testing_area\\102400060304_sca_10point_stats_results.json"

json_stats_file = open(json_file_path)
stats_data = json.load(json_stats_file)
 
#demo_df,demo_NSE=GullyStatsDataFrameConverter(data,"26915","070801050302","specific_contributing_area","2.6")

NSE_list=[]
csv_paths_list=[]
#%%
for projection_WKID, huc12s_list in huc_projection_dict.items():
    for huc12 in huc12s_list: 
        for terrain_index_name,current_threshold_values_list in chunmei_terrain_index_log_thresholds.items():#terrain_index_log_thresholds_dict.items():#
            for current_threshold_value in current_threshold_values_list:
                current_threshold_value_string=str(current_threshold_value)
                current_df,current_NSE=GullyStatsDataFrameConverter(stats_data,projection_WKID,huc12,terrain_index_name,current_threshold_value_string)
                current_csv_path=f"M:\\DEP_nrcs_gully_test6\\Gully_Outputs\\testing_area\\stats_csvs\\{huc12}_{terrain_index_name}_{current_threshold_value}.csv"
                projection_WKID_column=[]
                huc12_column=[]
                terrain_index_name_column=[]
                current_threshold_value_column=[]
                projection_WKID_column.extend(repeat(projection_WKID,len(current_df)))
                huc12_column.extend(repeat(huc12,len(current_df)))
                terrain_index_name_column.extend(repeat(terrain_index_name,len(current_df)))
                current_threshold_value_column.extend(repeat(current_threshold_value_string,len(current_df)))
                
                current_df["WKID"]=projection_WKID_column
                current_df["HUC12"]=huc12_column
                current_df["terrain_index"]=terrain_index_name_column
                current_df["threshold_value"]=current_threshold_value_column
                
                current_df.to_csv(current_csv_path,index=False)
                NSE_list.append([huc12,terrain_index_name,current_threshold_value,current_NSE])
                del(current_df)
                del(current_NSE)

NSE_df=pd.DataFrame(NSE_list, columns=["HUC12","Terrain Index type", "Terrain Index Threshold", "NSE"])
# NSE_file_name="M:\\DEP_nrcs_gully_test\\Gully_Outputs\\testing_area\\stats_csvs\\NSE_values.csv"
NSE_file_name="M:\\DEP_nrcs_gully_test6\\Gully_Outputs\\testing_area\\stats_csvs\\testing_minmax_stats_6m_snap_results.csv" #102400060304_sca_10point_NSE_values.csv"
NSE_df.to_csv(NSE_file_name,index=False)


current_df=pd.DataFrame(columns = ["WKID", "HUC12", "terrain_index","threshold_value",'unique_flowpath_number', 'gap_distance',"longest_upstream_path","longest_downstream_path","failure_value","HUC12_upstream_flowpath_mean","population_variance","error_variance"])
for projection_WKID, huc12s_list in huc_projection_dict.items():
    for huc12 in huc12s_list: 
        for terrain_index_name,current_threshold_values_list in chunmei_terrain_index_log_thresholds.items():#terrain_index_log_thresholds_dict.items():#
            for current_threshold_value in current_threshold_values_list:
                current_threshold_value_string=str(current_threshold_value)
                current_csv_path=f"M:\\DEP_nrcs_gully_test6\\Gully_Outputs\\testing_area\\stats_csvs\\{huc12}_{terrain_index_name}_{current_threshold_value}.csv"
                this_df=pd.read_csv(current_csv_path)
                current_df=pd.concat([current_df,this_df],ignore_index=True)#.append(this_df,ignore_index = True)
                csv_paths_list.append(current_csv_path)

grouped_mean_fps=current_df.groupby(by=["terrain_index","threshold_value"])["longest_upstream_path"].mean()
grouped_mean_fps.name="longest_upstream_path_mean"
# current_df["mean_upstream_flowpath_length"] = current_df.groupby(by=["terrain_index","threshold_value"])["longest_upstream_path"].mean()            
current_df=current_df.merge(grouped_mean_fps,on=["terrain_index","threshold_value"],how="outer")

current_df["population_variance"]=(current_df["longest_upstream_path"]-current_df["longest_upstream_path_mean"])**2
    # print(full_stats_df)
current_df["error_variance"]=np.where(current_df["gap_distance"].isnull(), current_df["failure_value"]**2,current_df["gap_distance"]**2)

NSE_denominator=current_df.groupby(by=["terrain_index","threshold_value"])["population_variance"].sum()
NSE_nummerator=current_df.groupby(by=["terrain_index","threshold_value"])["error_variance"].sum()
NSE_df=pd.merge(NSE_nummerator,NSE_denominator,on=["terrain_index","threshold_value"],how="outer")

NSE_df["NSE_value"]=1-(NSE_df["error_variance"]/NSE_df["population_variance"])
# NSE_value=float(1)-(NSE_nummerator/NSE_denominator)


current_df.to_csv("M:\\DEP_nrcs_gully_test6\\Gully_Outputs\\testing_area\\stats_csvs\\full_domain_data.csv",index=False)
NSE_df.to_csv("M:\\DEP_nrcs_gully_test6\\Gully_Outputs\\testing_area\\stats_csvs\\full_domain_NSE.csv")
#%%

json_errors_file_path="M:\\DEP_nrcs_gully_test6\\Gully_Outputs\\testing_area\\testing_minmax_stats_6m_snap_errors.json"
#102400060304_sca_10point_stats_errors.json#full_stats_errors.json"
json_errors_file = open(json_errors_file_path)
errors_data = json.load(json_errors_file)

for projection_WKID, huc12s_list in huc_projection_dict.items():
    for huc12 in huc12s_list: 
        for terrain_index_name,current_threshold_values_list in chunmei_terrain_index_log_thresholds.items():#terrain_index_log_thresholds_dict.items():# 
            for current_threshold_value in current_threshold_values_list:
                print(f"errors in huc12 {huc12} with {terrain_index_name} set to {current_threshold_value}:")
                for each_error in errors_data[projection_WKID][huc12][terrain_index_name][str(current_threshold_value)]:
                    print(f"gully {each_error[0]}: {each_error[1][1]}")

# Iterating through the
# list
# for i in data:
#     print(i)
#     for j in data[i]:
#         print(j)
#         for k in data[i][j]:
#             print(k)
#             for l in data[i][j][k]:
#                 print(l)
#                 for m in data[i][j][k][l]:
#                     print(m)

# stats_dataframe=pd.json_normalize(data)
# print(stats_dataframe)



json_stats_file.close()