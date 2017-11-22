#/bin/bash


if [ ! -f $1 ]
then
  exit -1
fi
 
if [ -f $1.h5 ]
then
  rm $1.h5
fi

sdds2stream $1 -col=t > tmp 

exec<tmp
a=0
while read line
do a=$(($a+1));
done


h5import tmp -dims $a -p t -t TEXTFP -s 64 -o $1.h5

sdds2stream $1 -col=p > tmp
h5import tmp -dims $a -p p -t TEXTFP -s 64 -o $1.h5

sdds2stream $1 -col=x > tmp
h5import tmp -dims $a -p x -t TEXTFP -s 64 -o $1.h5
sdds2stream $1 -col=y > tmp
h5import tmp -dims $a -p y -t TEXTFP -s 64 -o $1.h5

sdds2stream $1 -col=xp > tmp
h5import tmp -dims $a -p xp -t TEXTFP -s 64 -o $1.h5
sdds2stream $1 -col=yp > tmp
h5import tmp -dims $a -p yp -t TEXTFP -s 64 -o $1.h5

rm tmp