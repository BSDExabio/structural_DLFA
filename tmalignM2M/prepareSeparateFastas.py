import os, sys
from glob import glob

folder="/data/farhan/bench/input_fasta/"
outputfolder="/data/farhan/bench/input_fasta_separated/"

if not os.path.isdir(outputfolder): os.makedirs(outputfolder)

fasta_list=glob(folder+"*.fasta")

#fasta=fasta_list[3]
#print (fasta)

for fasta in fasta_list:
    name=os.path.basename(fasta).split(".")[0]
    #print (name)
    with open (fasta) as f:
        for line in f:
            #print (line)
            if line.startswith(">"):
                chains=line.strip().split("|")[1].split("[")[0].replace("Chains","").replace("Chain","").strip().split(",")
                sequence=f.readline()
                for chain in chains:
                    tag=">"+name+"_"+chain.strip()+"\n"
                    outdir=outputfolder+name+"/"
                    if not os.path.isdir(outdir): os.makedirs(outdir)
                    filename=outdir+name+"_"+chain.strip()+".fasta"
                    with open (filename,"w") as f_out:
                        f_out.write(tag)
                        f_out.write(sequence)