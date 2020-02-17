# energyScale


order of processing for LOFAR event analysis:


1) lofar_data/getEnergyEvents_final.py

  -looks in folders production_analysis_oct2018 for most recent iteration of reco_.dat
  -reads analysis info, finds sim runnr with closest xmax, calculates fluence, puts info in event files at 
    /vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/energyScale/radio
    
    
2) lofar_data/compileData_final.py
  - looks in /vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/energyScale/radio
  - collects all event info and writes 'compiled_sims_all_final.dat' file
  
3) compare_final.py
  -  collects fit info and 'uncertainty_info.dat' for final analysis



