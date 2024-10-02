# chmod +x run_mft_aqc_mw.sh
#!/bin/bash

run=553588
pass=apass1
hname=mMFTTrackPhi

root 'mft_moving_windows.cxx('$run',"'$pass'","'$hname'")'