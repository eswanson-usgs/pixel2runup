#rename wrongly named forecast images in S3
##### IMPORT #####
import os
import time
#will need fs3 package to use s3 in fsspec
import fsspec 
import imageio
import calendar
import datetime
import csv
import concurrent.futures


### FUNCTIONS ###


### MAIN ###
print("start:", datetime.datetime.now())

file_system = fsspec.filesystem('s3', profile='coastcam')
source_folder = "s3://cmgp-coastcam/cameras/madeira_beach/cx/2023/069_Mar.10/"
image_list = file_system.glob(source_folder+'/*')

for image in image_list:
    
    filename = image.split('/')[-1]
    image_type = filename.split('.')[8]

    if image_type == 'snap_forecast':
        destination = image.replace(image_type, 'snapForecast')
        file_system.copy(image, destination)
        file_system.delete(image, destination)
    
    elif image_type == 'timex_forecast':
        destination = image.replace(image_type, 'timexForecast')
        file_system.copy(image, destination)
        file_system.delete(image, destination)
    
    elif image_type == 'bright_forecast':
        destination = image.replace(image_type, 'brightForecast')
        file_system.copy(image, destination)
        file_system.delete(image, destination)


print("end:", datetime.datetime.now())
