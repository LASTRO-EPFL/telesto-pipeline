import os
import shutil
import astropy
from astropy.io import fits

# Folder creation
path = os.path.dirname(os.path.abspath('__Test_shuttle_tri.py__'))
if not os.path.exists(path+'\Deleted'):  #
    os.mkdir(path+'\Deleted')

files=os.listdir('.')
print(files)
for file in files:
    if "fit" in file:
        somme1=0
        somme2=0
        compteur1=0
        compteur2=0
        test_file = fits.getdata(file)
        numrow = len(test_file)
        numcol=len(test_file[0])
        print(file)
        bool=1
        #Checking if the shuttle is half open
        for i in range(int(numrow/2-numrow/15),int(numrow/2+numrow/15)):
            for j in range(int(numcol/2-numcol/15),int(numcol/2+numcol/15)): #Centered loop to check if the shuttle is half open
                compteur1=compteur1+1
                somme1=somme1+test_file[i,j]
        for i in range(0,int(numrow/8)):
            for j in range(0,int(numcol/8)): #loop on the edge of the image
                compteur2=compteur2+1
                somme2=somme2+test_file[i,j]
        mean_center=somme1/compteur1
        mean_outside=somme2/compteur2
        if mean_center>3200 and mean_center<4000 and mean_outside<2000:
            shutil.copyfile(file, path + '\Deleted\ ' + file)
            os.remove(file)
            bool=0
            print('1')
        #checking if the shuttle was closed
        if bool==1:
            if "Dark" not in file :
                somme=0
                compteur=0
                bool2=1
                for i in range(0, numrow):
                    for j in range(0, numcol):
                        if test_file[i,j]>3000:
                            bool2=0
                            break
                        else:
                            continue
                if bool2==1:
                    shutil.copyfile(file, path + '\Deleted\ ' + file)
                    os.remove(file)
                    print('2')

    '''
                        compteur = compteur+1
                        somme=somme+test_file[i,j]
                mean=somme/compteur
                if mean < 1200:
                    shutil.copyfile(file, path + '\Deleted\ ' + file)
                    os.remove(file)
                    print('2')
                    '''
