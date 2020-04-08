#!/bin/bash
set -x
source_path="${BASH_SOURCE[0]}"
source_dir=$(dirname "${source_path}")

rm "${source_dir}/classroom/scene_"*
rm "${source_dir}/classroom/timings"
rm "${source_dir}/classroom/samplerates"

#MDLC with matrix sep
for i in {1..5}
do
	mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.025 -Dcompcstrat=mdlc -Dcps=500 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -iet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/scene_amr_500"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_amr_500"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_amr_500"

for i in {1..5}
do
	mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.025 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -iet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/scene_amr_1000"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_amr_1000"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_amr_1000"

for i in {1..5}
do
	mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.025 -Dcompcstrat=mdlc -Dcps=2000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -iet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/scene_amr_2000"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_amr_2000"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_amr_2000"

for i in {1..5}
do
	mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.025 -Dcompcstrat=mdlc -Dcps=4000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -iet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/scene_amr_4000"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_amr_4000"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_amr_4000"

for i in {1..5}
do
	mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.025 -Dcompcstrat=mdlc -Dcps=8000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -iet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/scene_amr_8000"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_amr_8000"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_amr_8000"


#MDLC
for i in {1..5}
do
	mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=mdlc -Dcps=250 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -iet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/scene_mdlc_250"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_mdlc_250"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_mdlc_250"

for i in {1..5}
do
	mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=mdlc -Dcps=500 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -iet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/scene_mdlc_500"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_mdlc_500"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_mdlc_500"

for i in {1..5}
do
	mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -iet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/scene_mdlc_1000"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_mdlc_1000"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_mdlc_1000"

for i in {1..5}
do
	mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=mdlc -Dcps=2000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -iet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/scene_mdlc_2000"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_mdlc_2000"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_mdlc_2000"

for i in {1..5}
do
	mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=mdlc -Dcps=4000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -iet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/scene_mdlc_4000"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_mdlc_4000"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_mdlc_4000"


#LS
for i in {1..5}
do
	mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=ls -Dcps=250 -Dsps=10 -Dspp=4 -Dslice=5000 -Dis=1 -Dbv=1 -Dvp=0.025 -iet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/scene_ls_250"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_ls_250"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_ls_250"

for i in {1..5}
do
	mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=ls -Dcps=500 -Dsps=10 -Dspp=4 -Dslice=5000 -Dis=1 -Dbv=1 -Dvp=0.025 -iet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/scene_ls_500"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_ls_500"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_ls_500"

for i in {1..5}
do
	mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=ls -Dcps=1000 -Dsps=10 -Dspp=4 -Dslice=5000 -Dis=1 -Dbv=1 -Dvp=0.025 -iet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/scene_ls_1000"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_ls_1000"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_ls_1000"

for i in {1..5}
do
	mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=ls -Dcps=2000 -Dsps=10 -Dspp=4 -Dslice=5000 -Dis=1 -Dbv=1 -Dvp=0.025 -iet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/scene_ls_2000"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_ls_2000"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_ls_2000"

for i in {1..5}
do
	mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=ls -Dcps=4000 -Dsps=10 -Dspp=4 -Dslice=5000 -Dis=1 -Dbv=1 -Dvp=0.025 -iet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/scene_ls_4000"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_ls_4000"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_ls_4000"