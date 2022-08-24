#!/bin/bash

a=0.26
b=0.40
tm1=0.26
tm2=0.40
st=$((`echo "$a < $b"| bc`))

echo "$st"

stt=$((`echo "$tm1 < $tm2"| bc`))
echo "$stt"

comp=$((`echo "$tm1 < $tm2"| bc`))
echo "$comp"

if [ $((`echo "$a > $b" | bc`)) ]; then
	echo "bigger"
fi

tmscore1=0.2539
tmscore2=0.3596
#tm1=$(echo $tmscore1 | bc)
#tm2=$(echo $tmscore2 | bc)
tm1=0.26
tm2=0.26
echo "TM1=$tm1\t TM2=$tm2"
comp=$((`echo "$tm1 < $tm2"| bc`))
echo "$comp"
compeq=$((`echo "$tm1 == $tm2" | bc`))
echo "COMP= $comp \t COMPEQ=$compeq"

if [ $compeq ]; then
	echo "ONE $compeq"
fi
#echo $(( a < b ? a : b ))
