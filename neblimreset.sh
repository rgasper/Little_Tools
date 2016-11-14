#!/bin/bash
#this script takes an NEB and remakes it using the chosen image as the new limits
#        usage:
#       nebresetlims.sh initial final num_images cut_iteration_number
#       cut iteration is just for archive labeling purposes in case of multiple resets
usage="Usage: \n nebresetlims.sh initial final num_images cut_iteration_number"
argerr="ERROR: please input 4 arguments"

#error checking inputs
if [ -z $1 ]
        exit $E_MISSING_POS_PARAM
fi

if [ -z $2 ]
then
        echo ${argerr}
        echo ${usage}
        exit $E_MISSING_POS_PARAM
fi

if [ -z $3 ]
then
        echo ${argerr}
        echo ${usage}
        exit $E_MISSING_POS_PARAM
fi

if [ -z $4 ]
then
        echo ${argerr}
        echo ${usage}
        exit $E_MISSING_POS_PARAM
fi

#inputs fine

echo "Resetting your NEB. This will take a minute."

cp -rf $1 initial_cut_$4_$1
cp -rf $2 final_cut_$4_$2
rm -rf 0*
nebmake.pl initial_cut_$4_$1/CONTCAR final_cut_$4_$2/CONTCAR $3
nebavoid.pl 0.9 #this is my personal default, behaves well for hydrocarbon reactions on surfaces

echo "Done. Be sure to copy outcars into final and initial image directories!"
#EOF
