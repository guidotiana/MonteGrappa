#!/bin/bash

if [ $# -ne 7 ]; then
 echo "Generate input file for montegrappa from 5C/hiC data."
 echo " Input file in format: i j averagecount stdev"
 echo " Beads numbering i and j starts from zero."
 echo
 echo "Filename of the 5C/HiC data?"
 read ffi
 echo "Number of beads?"
 read n
 echo "Normalization constant?"
 read norm
 echo "Rootname for output files?"
 read rootname
 echo "Energy scale for initial interaction matrix (in temperature units)?"
 read escale
 echo "Interaction range (in units of interbeads distance)?"
 read rint
 echo "Hardcore radius  (in units of interbeads distance)?"
 read rhc
else
 ffi=$1
 n=$2
 norm=$3
 rootname=$4
 escale=$5 
 rint=$6
 rhc=$7
fi

dist=1
temp=1
echo 

#
# write polymer.pol
#

nm=$(( $n -1 ))
echo "[ backbone ]" > $rootname.pol 
for i in $(seq 0 $nm)
do
 x=`echo $i "/sqrt(2)*" $dist  | bc -l`;
 h=`echo $i"%2" | bc`;
 y=`echo $h"/sqrt(2)*"$dist  | bc -l`;
 echo $i $i "   P  " $i "   DNA   " $i "   0    " $x "  "$y"    0.00     1" >> $rootname.pol
done;

#
# write parameters.par
#

echo -e "nchains\t\t1" >$rootname.par
echo -e "nrun\t\t100" >>$rootname.par
echo -e "Temp\t\t"$temp >>$rootname.par
echo -e "nstep\t\t500000000" >>$rootname.par
echo -e "nprinttrj\t100000" >>$rootname.par
echo -e "nprinte\t\t100000" >>$rootname.par
echo -e "nprintlog\t100000" >>$rootname.par
echo -e "flip\t\t1" >>$rootname.par
#echo -e "mflip\t\t10" >>$rootname.par
echo -e "pivot\t\t10" >>$rootname.par
echo -e "seed\t\t-1" >>$rootname.par
echo -e "traj\t\t"$rootname  >>$rootname.par
echo "nosidechain"  >>$rootname.par
echo "nodihpot"  >>$rootname.par
echo "noangpot"  >>$rootname.par
echo -e "op_minim\tsample"   >>$rootname.par
echo -e "op_file\t\t"$rootname.op   >>$rootname.par
echo -e "op_deltat\t100000"  >>$rootname.par
echo -e "op_itermax\t1000"  >>$rootname.par
echo -e "op_step\t\t0.1"  >>$rootname.par
echo -e "op_T\t\t1.0"  >>$rootname.par
 
#
# process op file
#

nn=`wc -l $ffi` 
set -- $nn
echo "ndata "$1 > $rootname.op

awk '{ 
	x=$3/"'"$norm"'"; 
	if (x>1) x=1;
	sigmax = $4/"'"$norm"'";
	if (sigmax<0.01) sigmax=0.01 
	print $1,$2,0,x,sigmax; 
     }' $ffi > tmp.op

cat tmp.op >> $rootname.op

#
# write potential.pot
#

echo "[ global ]" > $rootname.pot
echo "hardcore    "$rhc >>  $rootname.pot

echo "[ pairs ]" >> $rootname.pot
awk 'BEGIN{
	escale="'"$escale"'";
	temp="'"$temp"'";
	rint="'"$rint"'";
	rhc="'"$rhc"'";
     }{ 
	ll=0+$2-$1;
	if (ll<0) ll=-ll;
	p=0+$4;
	if (p<0.95 && p>0.) e = -escale*temp*1.5*log(ll)-escale*temp*log(p/(1-p)) ;
	else if (p>0) e =  -escale*temp*1.5*log(ll)-escale*temp*2.94 ;
	else e=0;
	print $1,$2,e,rint,rhc;
     }'  tmp.op >> $rootname.pot 

		

#
# end
#

echo "I have generated for you:"
echo $rootname.pol
echo $rootname.par
echo $rootname.pot
echo $rootname.op
echo

