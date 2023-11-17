import os
unique_files=[]
large_files=[]
os.system("git status > change_in_git")
os.system("git add -n . > git_to_be_added")
with open("change_in_git") as f:
    lines=f.readlines()
    for line in lines:
        filename=line.split("\t")[-1].replace("\n","")
        #print(filename)
        if os.path.isfile(filename):
            megabyte=os.stat(filename).st_size / (1024 * 1024)
            #print(megabyte)
            if megabyte >20:
                raise ValueError(f"{filename}larger than 20" )
with open("git_to_be_added") as f:
    lines=f.readlines()
    for line in lines:
        filename=line.replace("add '","").replace("'","").replace("\n","")
        #print(filename)
        if os.path.isfile(filename):
            #print(1)
            megabyte=os.stat(filename).st_size / (1024 * 1024)
            print(filename,":",megabyte,"Mb")
            if megabyte >20:
                raise ValueError(f"{filename}larger than 20" )
                large_files.append(filename)
if len(large_files)>0:
    raise ValueError(large_files)
 
print("check done, all is not large file")
