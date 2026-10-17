#!/bin/bash
# Regression: compare the current build with the original 2.4b (+readouts) on fixed searches.
# usage: tools/regress.sh <dir with qfind.c/common.h> <path to original-with-readouts binary>
# Build the reference first:  gcc -std=c11 -fopenmp -O3 -o origr original-with-readouts/qfind.c
D=$1; O=$2; F="-std=c11 -fopenmp -march=native -O3"
gcc $F -o $D/r_d $D/qfind.c || exit 1
gcc $F -DPACKLISTS -o $D/r_p $D/qfind.c || exit 1
gcc $F -DLA_VERIFY -o $D/r_v $D/qfind.c || exit 1
gcc $F -DFACTORED_SORTED -o $D/r_fs $D/qfind.c || exit 1
ex(){ grep -E "Queue full|^[0-9bo\$!]+$|^x =|found|Maximum depth" | sed -E 's/^[0-9]+\/[0-9]+\/[0-9]+ [0-9:]+ //; s/ +\[.*\]//'; }
while read -r args; do
  [ -z "$args" ] && continue
  b=$($O $args 2>&1 | ex | md5sum); r=OK
  for x in r_d r_p r_v; do [ "$($D/$x $args 2>&1 | ex | md5sum)" = "$b" ] || r="$r DIFF($x)"; done
  [ "$($D/r_fs $args --tables factored 2>&1 | ex | md5sum)" = "$b" ] || r="$r DIFF(factored-sorted)"
  echo "$r | $args"
done <<'EOT'
-r B3/S23 -v c/5 -w 8 -s odd -f 1 -t 1
-r B3/S23 -v c/4 -w 6 -s even -f 1 -t 1
-r B3/S23 -v 2c/5 -w 6 -s odd -f 1 -t 1
-r B3/S23 -v c/3 -w 5 -s asym -f 1 -t 1
-r B3/S23 -v c/6 -w 7 -s even -f 1 -t 1
-r B36/S23 -v c/4 -w 6 -s odd -f 1 -t 1
-r B3/S23 -v c/7 -w 6 -s odd -f 1 -t 1
-r B3/S23 -v 2c/7 -w 6 -s even -t 1 -f 3
-r B3/S23 -v 2c/6 -w 6 -s even -f 2 -t 1 --disable-subperiod
-r B3/S23 -v c/4 -w 5 -s gutter -f 1 -t 1
-r B3/S23 -v c/3 -w 6 -s even -o odd -f 1 -t 1
-r B3/S23 -v c/4 -w 7 -s odd -t 1
EOT
