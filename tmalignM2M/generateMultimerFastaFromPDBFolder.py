import os,sys
from glob import glob

inpdir=os.path.abspath(sys.argv[1])+"/"
outputdir=os.path.abspath(sys.argv[2])+"/"

if not os.path.isdir(outputdir): os.makedirs(outputdir)

perl_script="/data/farhan/DeepComplex_runs/pdb2seq.pl "

folder_list=glob(inpdir+"*")

for pdb_folder in folder_list:
    pdb_list=glob(pdb_folder+"/*.pdb")
    name=os.path.basename(pdb_folder)
    if not os.path.isdir(outputdir+name+"/"):os.makedirs(outputdir+name+"/")
    #print (pdb_list)
    for pdb_name in pdb_list:
        chained_name=os.path.basename(pdb_name).replace(".pdb",".fasta")
        os.system("echo '>"+chained_name.replace(".fasta","")+"' > "+outputdir+name+"/"+chained_name)
        os.system("perl "+perl_script+pdb_name+" >> "+outputdir+name+"/"+chained_name)