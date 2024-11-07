# chmod +x run_mft_aqc_mw.sh
#!/bin/bash

run=559596
pass=online
path="qc/MFT/MO/MFTClusterTask/mw/" #qc_async/MFT/MO/Tracks/mw/
hname=mClustersROFSize
title="#Clusters/ROF"
option_hist="rebinROF logy logx"
option_plot="hist"
str_sor="STF"
str_eor="EOR"
rewrite_root=false

root -b -l <<EOF
.L mft_moving_windows.cxx
mft_moving_windows($run,"$pass","$path","$hname","$title","$option_hist","$option_plot","$str_sor","$str_eor",$rewrite_root)
.q
EOF