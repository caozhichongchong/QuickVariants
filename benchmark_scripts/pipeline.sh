#!/bin/bash
source ~/.bashrc
python SNP_model.py -i your_input_genome_dir
python SNPfilter.py
python Indelfilter.py
python SNP_model_compare.py

