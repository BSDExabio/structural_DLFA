#!/usr/bin/env bash
#
# Ingest SAdLSA data into the given database.
#
# SAdLSA output is generally stored in a set of directory hierarchies like this:
#
# {PROTEIN}/{SADLSA_RUN}/data_files
#
# E.g., WP_164928139.1/sadlsa_pdb70_210310/ will contain the SAdLSA score
# and alignment files for the sadlsa_pdb70_210310 SAdLSA run for protein
# WP_164928139.1.  The corresponding data files in that directory will be
# WP_164928139.1_aln.dat.gz and WP_164928139.1_sco.dat.gz.
