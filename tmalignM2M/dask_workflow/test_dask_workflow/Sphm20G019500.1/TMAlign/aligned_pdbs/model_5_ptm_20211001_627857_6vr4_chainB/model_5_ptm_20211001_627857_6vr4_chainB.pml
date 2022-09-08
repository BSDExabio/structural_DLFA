#!/usr/bin/env pymol
load /gpfs/alpine/bif135/proj-shared/rbd_work/dask_testing/tmalign_andes_workflow/test8/Sphm20G019500.1/TMAlign/aligned_pdbs/model_5_ptm_20211001_627857_6vr4_chainB/model_5_ptm_20211001_627857_6vr4_chainB, format=pdb
hide all
show stick
color blue, chain A
color red, chain B
set ray_shadow, 0
set stick_radius, 0.3
set sphere_scale, 0.25
show stick, not polymer
show sphere, not polymer
bg_color white
set transparency=0.2
zoom polymer

