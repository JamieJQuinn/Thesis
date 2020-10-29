#!/usr/bin/env bash

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null && pwd  )"

# Be sure to checkout correct branch or fail trying
DIR=$(pwd)

resolution=512
time_end=40.0
dt_snapshots=2.0

twisting_velocity=0.15
ramp_time=2.0

for n in 4; do
  for m in 0; do
    for visc in '-switching'; do
      for switching_param in 150; do
        folder=v-${n}r-$m$visc
        cp -rv $SCRIPT_DIR/lare3d $folder
        cd $folder
        sed -i -e 's/\(visc3 = 1.0e-\)[0-9]*/\1'$n'/' src/control.f90
        sed -i -e 's/\(eta_background =\) .*/\1'0.0_num'/' src/control.f90
        #sed -i -e 's/\(eta_background = 1.0e-\)[0-9]*/\1'$m'/' src/control.f90
        sed -i -e 's/\(switching_param = \).*/\1'${switching_param}'/' src/control.f90
        sed -i -e 's/\(nx_global = \).*/\1'${resolution}'/' src/control.f90
        sed -i -e 's/\(ny_global = \).*/\1'${resolution}'/' src/control.f90
        sed -i -e 's/\(nz_global = \).*/\1'${resolution}'/' src/control.f90
        sed -i -e 's/\(t_end = \).*\(_num\)/\1'${time_end}'\2/' src/control.f90
        sed -i -e 's/\(dt_snapshots = \).*\(_num\)/\1'${dt_snapshots}'\2/' src/control.f90
        sed -i -e 's/\(twisting_velocity = \).*\(_num\)/\1'${twisting_velocity}'\2/' src/control.f90
        sed -i -e 's/\(ramp_time = \).*\(_num\)/\1'${ramp_time}'\2/' src/control.f90
        #z_min_max=3.0
        #sed -i -e 's/\(z_min = \).*\(_num\)/\1'-${z_min_max}'\2/' src/control.f90
        #sed -i -e 's/\(z_max = \).*\(_num\)/\1'${z_min_max}'\2/' src/control.f90
        #sed -i -e 's/\(initial = \).*/\1IC_RESTART/' src/control.f90
        #sed -i -e 's/\(restart_snapshot = \).*/\130/' src/control.f90
        cd ..
      done
    done
  done
done
