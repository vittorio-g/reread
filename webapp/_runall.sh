cd /c/Users/vitto/Desktop/ReReReRe/webapp
L=/c/Users/vitto/AppData/Local/Temp/claude/C--Users-vitto-Desktop-ReReReRe/8f426989-5c28-4009-a838-215f9a1e2ca0/scratchpad/reruns.log
: > "$L"
until [ -f clean_envelope_ab.csv ]; do sleep 15; done
echo "=== STUDY 1 (gateless vs gated, same resamples) ===" >> "$L"
node study1_gateless.js >> "$L" 2>&1
echo "" >> "$L"; echo "=== STUDY 2 (external, deployable) ===" >> "$L"
node ext_deployable.js >> "$L" 2>&1
echo "" >> "$L"; echo "=== STUDY 3 (simulation) ===" >> "$L"
node sim_shipped.js >> "$L" 2>&1
echo "" >> "$L"; echo "=== DEMO ===" >> "$L"
node extract_demo_study1.js >> "$L" 2>&1
node check_demo_study1.js >> "$L" 2>&1
echo "ALL_RERUNS_DONE" >> "$L"
