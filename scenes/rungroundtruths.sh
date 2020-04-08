#!/bin/bash
source_path="${BASH_SOURCE[0]}"
source_dir=$(dirname "${source_path}")

mitsuba "${source_dir}/cbox/cbox.xml" -Dcstrat=0 -Dsamplerate=0 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -o "${source_dir}/cbox/groundtruth.exr"
#mitsuba "${source_dir}/living_room_3/living-room-3/scene.xml" -Dcstrat=0 -Dsamplerate=0 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -o "${source_dir}/living_room_3/living-room-3/groundtruth.exr"
#mitsuba "${source_dir}/san_miguel/san-miguel.xml" -Dcstrat=0 -Dsamplerate=0 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -o "${source_dir}/san_miguel/groundtruth.exr"
#mitsuba "${source_dir}/staircase/staircase/scene.xml" -Dcstrat=0 -Dsamplerate=0 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -o "${source_dir}/staircase/staircase/groundtruth.exr"
#mitsuba "${source_dir}/breakfast-room/breakfast-room/scene.xml" -Dcstrat=0 -Dsamplerate=0 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -o "${source_dir}/breakfast-room/breakfast-room/groundtruth.exr"
mitsuba "${source_dir}/hairball/hairball.xml" -Dcstrat=0 -Dsamplerate=0 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -o "${source_dir}/hairball/groundtruth.exr"
#mitsuba "${source_dir}/sponza/sponza.xml" -Dcstrat=0 -Dsamplerate=0 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -o "${source_dir}/sponza/groundtruth.exr"
#mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=0 -Dsamplerate=0 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -o "${source_dir}/classroom/groundtruth.exr"
#mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=0 -Dsamplerate=0 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -o "${source_dir}/kitchen/groundtruth.exr"
#mitsuba "${source_dir}/bathroom/bathroom/scene.xml" -Dcstrat=0 -Dsamplerate=0 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -o "${source_dir}/bathroom/bathroom/groundtruth.exr"
