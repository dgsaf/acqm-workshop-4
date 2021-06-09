#!/bin/bash

# input, output directory
input="data/input"
output="data/output"

# input file
input_file="${input}/data.in"

# constant parameters
dr="0.001"
rmax="200"
lmin="0"

# cons_input ()
# will overwrite existing ${input} file
# usage: cons_input zproj lmax energy
cons_input () {
  echo "$3" > ${input_file}
  echo "${rmax} , ${dr}" >> ${input_file}
  echo "$1 , ${lmin} , $2" >> ${input_file}
}

# cons_str ()
# usage: cons_str zproj lmax energy
cons_str () {
  printf "z-%+i.l-%i.en-%.1f" $1 $2 $3
}

# store_output ()
# usage: store_output zproj lmax energy
store_output () {
  local str=$(cons_str $1 $2 $3)

  mv "${output}/ics.txt" "${output}/ics.${str}.txt"
  mv "${output}/dcs.txt" "${output}/dcs.${str}.txt"
}

# cons_linear_set ()
# usage: cons_linear_set start final
cons_linear_set () {
  local set=""
  for (( ii = $1 ; ii <= $2 ; ii++ )) ; do
    set="${set} ${ii}"
  done
  echo "${set}"
}

# cons_quadratic_set ()
# usage: cons_quadratic_set start final n
cons_quadratic_set () {
  local a=$(awk '{print ($2 - $1)/(($3*$3) - 1.0)}' <<< "$1 $2 $3")
  local set=""
  for (( ii = 1 ; ii <= $3 ; ii++ )) ; do
    set="${set} $(awk '{printf "%.4f", $1 + $2*(($3*$3) - 1.0)}' <<< "$1 ${a} ${ii}")"
  done
  echo "${set}"
}

# varying parameters
zproj_set="-1 +1"
lmax_set=$(cons_linear_set 0 5)
energy_set=$(cons_quadratic_set 0.1 50.0 20)

# compile
make bin/main

# perform jobs
echo "potential scattering jobs"
echo "constant parameters: "
echo "> dr = ${dr}"
echo "> rmax = ${rmax}"
echo "> lmin = ${lmin}"
echo "varying parameters: "
echo "> zproj_set: ${zproj_set}"
echo "> lmax_set: ${lmax_set}"
echo "> energy_set: ${energy_set}"

# energy jobs (for electron, and positron)
echo ""
echo "energy jobs:"

# set lmax to highest
lmax=$(echo "${lmax_set}" | awk '{print $NF}')

for zproj in ${zproj_set} ; do
  for energy in ${energy_set} ; do
    # announce job parameters
    echo ""
    echo "> lmax = ${lmax} , zproj = ${zproj} , energy = ${energy}"

    # skip if job already performed
    str=$(cons_str ${zproj} ${lmax} ${energy})
    if [ -f "${output}/ics.${str}.txt" ] ; then
      echo "> job already performed"
      continue
    fi

    # construct input file
    cons_input ${zproj} ${lmax} ${energy}

    # execute job
    t_start=$(date +%s)
    bin/main > /dev/null
    t_end=$(date +%s)
    printf "> %is elapsed\n" $((t_end - t_start))

    # relocate ics, dcs output files
    store_output ${zproj} ${lmax} ${energy}
  done

  # collect all ics across energies into one file for given zproj
  collect=$(printf "${output}/collected_ics.z-%+i.l-%i.txt" ${zproj} ${lmax})
  > ${collect}
  for energy in ${energy_set} ; do
    str=$(cons_str ${zproj} ${lmax} ${energy})
    ics_file="${output}/ics.${str}.txt"

    partial_ics=$(grep -v '^[ ]*#' ${ics_file} | awk '{printf "%lf ", $2}' | tr -s '\n' ' ')
    total_ics=$(tail -n 1 ${ics_file} | awk '{printf $3}')

    echo "$energy ${total_ics} ${partial_ics}" >> ${collect}
  done
done

# lmax job (for electron)
echo ""
echo "lmax jobs:"

# set lmax higher
lmax_set=$(cons_linear_set 0 9)

# set zproj to -1
zproj="-1"

# set energy to 25.0 eV
energy="25.0"

for lmax in ${lmax_set} ; do
  # announce job
  echo ""
  echo "zproj = ${zproj} , energy = ${energy} , lmax = ${lmax} "

  # skip if job already performed
  str=$(cons_str ${zproj} ${lmax} ${energy})
  if [ -f "${output}/ics.${str}.txt" ] ; then
    echo "job already performed"
    continue
  fi

  # construct input file
  cons_input ${zproj} ${lmax} ${energy}

  # execute job
  t_start=$(date +%s)
  bin/main > /dev/null
  t_end=$(date +%s)
  printf "> %is elapsed\n" $((t_end - t_start))

  # relocate ics, dcs output files
  store_output ${zproj} ${lmax} ${energy}
done
