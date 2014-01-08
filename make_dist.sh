#!/bin/bash

prog=SurfaceReconstructionLibrary;


cd /tmp;
cvs -d :ext:oph@serpentor.bu.edu:/home/cvs export -D now PROJECTS/"$prog";
cd PROJECTS/"$prog"/doc;

version=`cat libsr.doxy | sed -n /PROJECT_NUMBER/p | sed "s/^.* = //"`;
name="$prog"-"$version";

doxygen libsr.doxy;
cd ../;

rm -rf notes
rm TODO;
rm make_dist.sh;
cd ../;

mv "$prog" "$name"
tar -czvf "$name".tgz "$name";
cp "$name".tgz $HOME/projects/html/external/software/"$prog"/;

cd $HOME/projects/html/external/software/"$prog"/;
cvs add "$name".tgz;
cvs commit -m "" "$name".tgz;

echo "now shell into eslab and update the website";