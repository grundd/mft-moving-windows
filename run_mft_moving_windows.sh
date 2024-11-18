# chmod +x run_mft_aqc_mw.sh
#!/bin/bash

run=559544
pass=online
path="qc/MFT/MO/MFTClusterTask/mw/" #qc_async/MFT/MO/Tracks/mw/
hname=mClustersROFSize
title="#Clusters/ROF"
option_hist="rebinROF logy logx"
option_plot="hist"
str_sor="STF"
str_eor="EOR"
plot_next=true
rewrite_root=false
aggr_histo=0

root -b -l <<EOF
.L mft_moving_windows.cxx
mft_moving_windows($run,"$pass","$path","$hname","$title","$option_hist","$option_plot","$str_sor","$str_eor",$plot_next,$rewrite_root,$aggr_histo)
.q
EOF