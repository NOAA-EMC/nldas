#!  /bin/csh -fx

cd make/MAKDEP
gmake
cd ..
gmake clean
gmake
