#!bin/bash

cmsDriver.py SingleMuPt1000_pythia8_cfi.py --step GEN,SIM,RECO --fileout file:Muons_PT_10_1500.root --mc --datatier RECOSIM --conditions auto:mc --python_filename make_FlatMuonPTMC.py --no_exec -n 100000
