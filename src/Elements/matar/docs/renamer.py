import os 
import sys

def rename(path):
    files = os.listdir(path)
    for file in files:
        if(file in [".", "./", "..", "../"]):
            continue
        if( os.path.isdir(file)):
            rename(path+"/"+file)
        if(".html" not in file):
            continue
        test = path + "/" + file 

        file = open(test, "r")    
        data = file.read()
        file.close()

        file = open(test, "w")
        newdata = data.replace("_static", "static")
        file.write(newdata)
        file.close()




rename(".")
