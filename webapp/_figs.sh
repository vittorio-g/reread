cd /c/Users/vitto/Desktop/ReReReRe/webapp
L=/c/Users/vitto/AppData/Local/Temp/claude/C--Users-vitto-Desktop-ReReReRe/8f426989-5c28-4009-a838-215f9a1e2ca0/scratchpad/figs.log
: > "$L"
echo "== eval_sens2_ab (sensitivity + multimetric data) ==" >> "$L"
node eval_sens2_ab.js >> "$L" 2>&1
echo "== gen_fig_calib (calibration data) ==" >> "$L"
node gen_fig_calib.js >> "$L" 2>&1
echo "== plots ==" >> "$L"
python plot_sens2.py >> "$L" 2>&1
python plot_multimetric2.py >> "$L" 2>&1
python plot_fig_calib.py >> "$L" 2>&1
python plot_robust.py >> "$L" 2>&1
echo "FIGS_DONE" >> "$L"
