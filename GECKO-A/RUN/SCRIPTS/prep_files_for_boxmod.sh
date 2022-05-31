#! /bin/bash 

#path to the box model
box=../../../BOX1.20_BA_20jan/

if [ $# -ne 2 ] 
then
  echo argument missing
  exit
fi

[ -e ../OUT/reactions.dum ] || echo 'no mechanism in this directory' 
[ -e ../OUT/reactions.dum ] || exit 

cd ../
echo
echo '--------------------------------'
echo 'writing files for simulation ...'
echo '--------------------------------'
echo 'make the *.mech file'
cp ./SCRIPTS/makemech.sh ./OUT/
cd ./OUT 
./makemech.sh $1
rm ./makemech.sh

echo '... cp  dictionary.out *.dict'
cp dictionary.out $1.dict

echo ... cp pvap.nan.dat $1.nan.sat
cp pvap.nan.dat $1.nan.sat
echo ... cp pvap.sim.dat $1.sim.sat
cp pvap.sim.dat $1.sim.sat
echo ... cp pvap.myr.dat $1.myr.sat
cp pvap.myr.dat $1.myr.sat

echo ... cp henry_gromhe.dat $1.hlc
cp henry_gromhe.dat $1.hlc

echo ... cp kicovi.dat $1.koh
cp kicovi.dat $1.koh

echo ... cp Tg.dat $1.Tg
cp Tg.dat $1.Tg

echo ... cp henry.dep $1.dep
cp henry.dep $1.dep
mv $1.dep $box'INPUT'

echo
echo --------------------------------
echo searching code name of the input formula
echo --------------------------------
name=$(cut -d" " -f2  listprimary.dat)
echo $name

echo
echo --------------------------------
echo move files .....
echo --------------------------------
echo ... mv $1'.*' $box'CHEMDAT'
mv $1.* $box'CHEMDAT'

dirsav=$box'CHEMDAT/COUNTERS/'$2
[ -e $dirsav ] || mkdir $dirsav && echo 'Directory already exists'

echo ... cp 'pero*' $box'CHEMDAT/COUNTERS/'$2'/'
cp ./pero* $box'CHEMDAT/COUNTERS/'$2'/'
[ -e indat1.dim ] && cp indat*.dim $box'CHEMDAT/COUNTERS/'$2'/'

# ---------------------------------------------------------------------
#----------  NOW in boxmod directory
# ---------------------------------------------------------------------

echo
echo ---------------------------------
echo interpreting $1.mech
echo ---------------------------------
cd $box'CHEMDAT'
./runintp.sh $1 > $1_interp.out

echo
echo ---------------------------------
echo creating $2 directory for simulation
echo ---------------------------------

cd ..

[ -e $box$2 ] && echo 'Directory already exists' || cp -r XXXXX $2

cd $2

echo ... modifying start.sh
sed -i s/xxxxx/$1/g start.sh
sed -i s/XXXXX/$2/g start.sh

echo ... modifying start.key
sed -i "s/xxxxx/G$name/g" start.key
echo
echo ... modifying precu.nam
sed -i "s/xxxxx/G$name/g" precu.nam
echo

echo ... modifying plot.key
sed -i "s/xxxxx/G$name/g" plot.key
echo
echo enter the number of carbon of the parent species
read numc
let numc2=$numc*4
sed -i "s/XX/$numc2/g" plot.key

echo ... done
