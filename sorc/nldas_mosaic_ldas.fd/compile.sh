#!  /bin/csh -fx

cd make/MAKDEP
gmake
cd ..
cp ./MAKDEP/makdep makdep 
gmake clean
gmake
