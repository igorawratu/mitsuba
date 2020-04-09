#!/bin/bash
set -x
source_path="${BASH_SOURCE[0]}"
source_dir=$(dirname "${source_path}")

rm "${source_dir}/bathroom/bathroom/scene_"*
rm "${source_dir}/bathroom/bathroom/timings"
rm "${source_dir}/bathroom/bathroom/samplerates"

#MDLC with matrix sep
for i in {1..5}
do
	mitsuba "${source_dir}/bathroom/bathroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.025 -Dcompcstrat=mdlc -Dcps=500 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -Diet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/bathroom/bathroom/scene.exr" "${source_dir}/bathroom/bathroom/groundtruth.exr" "${source_dir}/bathroom/bathroom/scene_amr_500"
done
mv "${source_dir}/bathroom/bathroom/timings" "${source_dir}/bathroom/bathroom/timings_amr_500"
mv "${source_dir}/bathroom/bathroom/samplerates" "${source_dir}/bathroom/bathroom/samplerates_amr_500"

for i in {1..5}
do
	mitsuba "${source_dir}/bathroom/bathroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.025 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -Diet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/bathroom/bathroom/scene.exr" "${source_dir}/bathroom/bathroom/groundtruth.exr" "${source_dir}/bathroom/bathroom/scene_amr_1000"
done
mv "${source_dir}/bathroom/bathroom/timings" "${source_dir}/bathroom/bathroom/timings_amr_1000"
mv "${source_dir}/bathroom/bathroom/samplerates" "${source_dir}/bathroom/bathroom/samplerates_amr_1000"

for i in {1..5}
do
	mitsuba "${source_dir}/bathroom/bathroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.025 -Dcompcstrat=mdlc -Dcps=2000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -Diet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/bathroom/bathroom/scene.exr" "${source_dir}/bathroom/bathroom/groundtruth.exr" "${source_dir}/bathroom/bathroom/scene_amr_2000"
done
mv "${source_dir}/bathroom/bathroom/timings" "${source_dir}/bathroom/bathroom/timings_amr_2000"
mv "${source_dir}/bathroom/bathroom/samplerates" "${source_dir}/bathroom/bathroom/samplerates_amr_2000"

for i in {1..5}
do
	mitsuba "${source_dir}/bathroom/bathroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.025 -Dcompcstrat=mdlc -Dcps=4000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -Diet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/bathroom/bathroom/scene.exr" "${source_dir}/bathroom/bathroom/groundtruth.exr" "${source_dir}/bathroom/bathroom/scene_amr_4000"
done
mv "${source_dir}/bathroom/bathroom/timings" "${source_dir}/bathroom/bathroom/timings_amr_4000"
mv "${source_dir}/bathroom/bathroom/samplerates" "${source_dir}/bathroom/bathroom/samplerates_amr_4000"

for i in {1..5}
do
	mitsuba "${source_dir}/bathroom/bathroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.025 -Dcompcstrat=mdlc -Dcps=8000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -Diet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/bathroom/bathroom/scene.exr" "${source_dir}/bathroom/bathroom/groundtruth.exr" "${source_dir}/bathroom/bathroom/scene_amr_8000"
done
mv "${source_dir}/bathroom/bathroom/timings" "${source_dir}/bathroom/bathroom/timings_amr_8000"
mv "${source_dir}/bathroom/bathroom/samplerates" "${source_dir}/bathroom/bathroom/samplerates_amr_8000"


#MDLC
for i in {1..5}
do
	mitsuba "${source_dir}/bathroom/bathroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=mdlc -Dcps=250 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -Diet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/bathroom/bathroom/scene.exr" "${source_dir}/bathroom/bathroom/groundtruth.exr" "${source_dir}/bathroom/bathroom/scene_mdlc_250"
done
mv "${source_dir}/bathroom/bathroom/timings" "${source_dir}/bathroom/bathroom/timings_mdlc_250"
mv "${source_dir}/bathroom/bathroom/samplerates" "${source_dir}/bathroom/bathroom/samplerates_mdlc_250"

for i in {1..5}
do
	mitsuba "${source_dir}/bathroom/bathroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=mdlc -Dcps=500 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -Diet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/bathroom/bathroom/scene.exr" "${source_dir}/bathroom/bathroom/groundtruth.exr" "${source_dir}/bathroom/bathroom/scene_mdlc_500"
done
mv "${source_dir}/bathroom/bathroom/timings" "${source_dir}/bathroom/bathroom/timings_mdlc_500"
mv "${source_dir}/bathroom/bathroom/samplerates" "${source_dir}/bathroom/bathroom/samplerates_mdlc_500"

for i in {1..5}
do
	mitsuba "${source_dir}/bathroom/bathroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -Diet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/bathroom/bathroom/scene.exr" "${source_dir}/bathroom/bathroom/groundtruth.exr" "${source_dir}/bathroom/bathroom/scene_mdlc_1000"
done
mv "${source_dir}/bathroom/bathroom/timings" "${source_dir}/bathroom/bathroom/timings_mdlc_1000"
mv "${source_dir}/bathroom/bathroom/samplerates" "${source_dir}/bathroom/bathroom/samplerates_mdlc_1000"

for i in {1..5}
do
	mitsuba "${source_dir}/bathroom/bathroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=mdlc -Dcps=2000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -Diet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/bathroom/bathroom/scene.exr" "${source_dir}/bathroom/bathroom/groundtruth.exr" "${source_dir}/bathroom/bathroom/scene_mdlc_2000"
done
mv "${source_dir}/bathroom/bathroom/timings" "${source_dir}/bathroom/bathroom/timings_mdlc_2000"
mv "${source_dir}/bathroom/bathroom/samplerates" "${source_dir}/bathroom/bathroom/samplerates_mdlc_2000"

for i in {1..5}
do
	mitsuba "${source_dir}/bathroom/bathroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=mdlc -Dcps=4000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.025 -Diet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/bathroom/bathroom/scene.exr" "${source_dir}/bathroom/bathroom/groundtruth.exr" "${source_dir}/bathroom/bathroom/scene_mdlc_4000"
done
mv "${source_dir}/bathroom/bathroom/timings" "${source_dir}/bathroom/bathroom/timings_mdlc_4000"
mv "${source_dir}/bathroom/bathroom/samplerates" "${source_dir}/bathroom/bathroom/samplerates_mdlc_4000"


#LS
for i in {1..5}
do
	mitsuba "${source_dir}/bathroom/bathroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=ls -Dcps=250 -Dsps=5 -Dspp=4 -Dslice=5000 -Dis=1 -Dbv=1 -Dvp=0.025 -Diet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/bathroom/bathroom/scene.exr" "${source_dir}/bathroom/bathroom/groundtruth.exr" "${source_dir}/bathroom/bathroom/scene_ls_250"
done
mv "${source_dir}/bathroom/bathroom/timings" "${source_dir}/bathroom/bathroom/timings_ls_250"
mv "${source_dir}/bathroom/bathroom/samplerates" "${source_dir}/bathroom/bathroom/samplerates_ls_250"

for i in {1..5}
do
	mitsuba "${source_dir}/bathroom/bathroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=ls -Dcps=500 -Dsps=5 -Dspp=4 -Dslice=5000 -Dis=1 -Dbv=1 -Dvp=0.025 -Diet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/bathroom/bathroom/scene.exr" "${source_dir}/bathroom/bathroom/groundtruth.exr" "${source_dir}/bathroom/bathroom/scene_ls_500"
done
mv "${source_dir}/bathroom/bathroom/timings" "${source_dir}/bathroom/bathroom/timings_ls_500"
mv "${source_dir}/bathroom/bathroom/samplerates" "${source_dir}/bathroom/bathroom/samplerates_ls_500"

for i in {1..5}
do
	mitsuba "${source_dir}/bathroom/bathroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=ls -Dcps=1000 -Dsps=5 -Dspp=4 -Dslice=5000 -Dis=1 -Dbv=1 -Dvp=0.025 -Diet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/bathroom/bathroom/scene.exr" "${source_dir}/bathroom/bathroom/groundtruth.exr" "${source_dir}/bathroom/bathroom/scene_ls_1000"
done
mv "${source_dir}/bathroom/bathroom/timings" "${source_dir}/bathroom/bathroom/timings_ls_1000"
mv "${source_dir}/bathroom/bathroom/samplerates" "${source_dir}/bathroom/bathroom/samplerates_ls_1000"

for i in {1..5}
do
	mitsuba "${source_dir}/bathroom/bathroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=ls -Dcps=2000 -Dsps=5 -Dspp=4 -Dslice=5000 -Dis=1 -Dbv=1 -Dvp=0.025 -Diet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/bathroom/bathroom/scene.exr" "${source_dir}/bathroom/bathroom/groundtruth.exr" "${source_dir}/bathroom/bathroom/scene_ls_2000"
done
mv "${source_dir}/bathroom/bathroom/timings" "${source_dir}/bathroom/bathroom/timings_ls_2000"
mv "${source_dir}/bathroom/bathroom/samplerates" "${source_dir}/bathroom/bathroom/samplerates_ls_2000"

for i in {1..5}
do
	mitsuba "${source_dir}/bathroom/bathroom/scene.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=ls -Dcps=4000 -Dsps=5 -Dspp=4 -Dslice=5000 -Dis=1 -Dbv=1 -Dvp=0.025 -Diet=0
	python3 "${source_dir}/../compare_err.py" "${source_dir}/bathroom/bathroom/scene.exr" "${source_dir}/bathroom/bathroom/groundtruth.exr" "${source_dir}/bathroom/bathroom/scene_ls_4000"
done
mv "${source_dir}/bathroom/bathroom/timings" "${source_dir}/bathroom/bathroom/timings_ls_4000"
mv "${source_dir}/bathroom/bathroom/samplerates" "${source_dir}/bathroom/bathroom/samplerates_ls_4000"