#!/bin/bash
#this script takes an NEB and interpolates, adding extra images in between the two you pick.
#Currently requires one to manually move the directories first
# example: if you have 00 01 02 03 04 05 and want to add 2 images in between 3 and 4
#               must first move 4 and 5 to 06 and 07;
#it works via a bit of chicanery with henkelman's nebmake, but it works

usage="Usage: \n nebinterp.sh leftdir rightdir num_images"
argerr="ERROR: Please input 3 arguments"

#error checking inputs

if [ -z $1 ]
then
        echo ${argerr}
        echo ${usage}
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
#inputs fine
#tell the user
echo "Interpolating between images $1 and $2 with $3 new images. This may take a minute."
cat >> notes <<!
Interpolated between images $1 and $2 with $3 new images.
!
#nebmake.pl in false environment
mkdir tmpdir
cd tmpdir
cp -rf ../$1 init
cp -rf ../$2 fin
nebmake.pl init/CONTCAR fin/CONTCAR $3

echo "Done. Please copy the interpolated configurations into the proper directories."
#EOF
