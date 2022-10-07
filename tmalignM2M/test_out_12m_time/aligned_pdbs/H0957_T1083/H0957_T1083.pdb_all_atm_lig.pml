#!/usr/bin/env pymol
load /data/farhan/Collaboration/structural_DLFA/tmalignM2M/test_out_12m_time/aligned_pdbs/H0957_T1083/H0957_T1083.pdb_all_atm_lig, format=pdb
hide all
show cartoon
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

