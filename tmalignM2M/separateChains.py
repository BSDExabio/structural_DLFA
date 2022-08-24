import os,sys
from glob import glob

inpdir=os.path.abspath(sys.argv[1])+"/"
outputdir=os.path.abspath(sys.argv[2])+"/"

if not os.path.isdir(outputdir): os.makedirs(outputdir)

perl_script="/data/farhan/DeepComplex_runs/pdb2seq.pl "

file_list=glob(inpdir+"*.pdb")

print (os.path.basename(file_list[0])[0:5])
for pdb in file_list:
    name=os.path.basename(pdb)[0:5]
    os.system("cp "+pdb+" "+outputdir+name+".pdb")
    print (name)
#print (len(file_list))
#with open ()