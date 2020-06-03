set -x
source_path="${BASH_SOURCE[0]}"
source_dir=$(dirname "${source_path}")

for i in {1..50}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=$((0.01 * i)) -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolge_1_1000"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolge_1_1000"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolge_1_1000"