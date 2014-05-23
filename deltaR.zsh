for i in $(seq 1 100); do
      ./montegrappa pol/unfold.pol pot/1PGB.pot provatot.par 2> err.err;
      grep DELTAR err.err > DELTAR.dat;
      cat DELTAR.dat | awk '{sum+=$2+$3+$4} END {print sum/NR}';
done
