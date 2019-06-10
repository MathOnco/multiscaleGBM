#!/usr/bin/env bash
mkdir ../classes
rm ../classes/*.class
javac Main.java -d ../classes

cd ../classes

#variables
tN=1;#tumor id, 1=nodular, 2=intermediate, 3=diffuse, 4=fit to all metrics (heterogeneous), 5=fit to all metrics (homogeneous)
tx=0;#0=grow only, 1=anti-proliferative Tx, 2=anti-migratory Tx, 3=AP+AM Tx

#toggle collection of movies, tracks, and single cell phenotypes
mov=1;#movie: 0=off, 1=on
tracks=1;#tracks: 0=off, 1=on
phenos=0;#phenotypes: 0=off, 1=on

#import parameters
parameters=`cat '../source/input/params.txt' | tail -n +$tN | head -1`

java Main ${parameters} $tN $tx $mov $tracks $phenos


