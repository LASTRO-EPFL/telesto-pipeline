import os
import shutil

# Folder creation
path = os.path.dirname(os.path.abspath('__Sorting.py__'))
if not os.path.exists(path+'\images'):  #
    os.mkdir(path+'\images')
if not os.path.exists(path+'\Flat_field'):
    os.mkdir(path+'\Flat_field')
if not os.path.exists(path+'\Dark_frame'):
    os.mkdir(path+'\Dark_frame')
if not os.path.exists(path+'\Illumination'):
    os.mkdir(path+'\Illumination')

# Moving files into correct folder
files=os.listdir('.')
print(files)
for file in files:
    if "fit" in file:
        if "FlatField" in file:
            shutil.copyfile(file,path+'\Flat_field\ '+file)
        elif "Dark" in file:
            shutil.copyfile(file,path+'\Dark_frame\ '+file)
        elif "illumination" in file:
            shutil.copyfile(file, path + '\Illumination\ ' + file)
        else:
            shutil.copyfile(file, path + '\images\ ' + file)
        os.remove(file)
