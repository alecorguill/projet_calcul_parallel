#!/bin/bash

file=speedup.csv
rm -f $file

# for i in {1..3..1}
# do
#     /usr/bin/time -f "%e" mpirun --mca pml ob1 -np $i run config.cfg >> $file
# done

make
for i in {2..8..1}
do
    echo -n "$i " >> $file  
    { /usr/bin/time -f "%e" mpirun --mca pml ob1 -np $i run config.cfg ; } 2>> $file
done
