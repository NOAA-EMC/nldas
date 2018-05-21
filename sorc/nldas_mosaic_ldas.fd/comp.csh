#!  /bin/csh -fx

rm LISM temp* 
cd make/MAKDEP
gmake
cd ..
rm LIS
gmake
cp LIS ../LISM
cd ..

