#!/usr/bin/env pymol
load /data/farhan/Collaboration/structural_DLFA/tmalignM2M/cif_test_out/aligned_pdbs/1bxw_1bxw_cif/1bxw_1bxw_cif.pdb, format=pdb
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

