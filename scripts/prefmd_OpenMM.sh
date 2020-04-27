#!/bin/bash
# carry out structure refinement for a single chain on one GPU
# Michael Feig, Michigan State University, 2017

if [ $# -eq 0 ]; then
    echo "usage: prefmd.sh [options] INPUT_PDB"
    exit -1
fi

CHARMMEXECSINGLE=$CHARMMEXEC      # serial executable
CHARMMEXECMPI=$CHARMMEXEC         # MPI executable


topfile=`pwd`/top.rtf # see README.md
parfile=`pwd`/par.prm # see README.md

cpus=8
gpu=0
mdruns=5
mdsteps=15000000    # 30ns
#tmpdir=/tmp
tmpdir=.
kcons=0.05
bottom=0.0

while :; do
  case $1 in 
    -c|--cpus) 
      if [ -n "$2" ]; then
        cpus=$2
        shift
      fi
      ;;
    -g|--gpu)
      if [ -n "$2" ]; then
        gpu=$2
        shift
      fi
      ;;
    --mdsteps)
      if [ -n "$2" ]; then
        mdsteps=$2
        shift
      fi
      ;;
    --kcons)
      if [ -n "$2" ]; then
        kcons=$2
        shift
      fi
      ;;
    --flat_bottom)
      if [ -n "$2" ]; then
        bottom=$2
        shift
      fi
      ;;
    --mdruns)
      if [ -n "$2" ]; then
        mdruns=$2
        shift
      fi
      ;;
    --tmpdir)
      if [ -n "$2" ]; then
        tmpdir=$2
        shift
      fi
      ;;
    --)
      shift
      break
      ;;
    *)
      break
  esac
  shift
done

if [[ $bottom < 0.0 ]]; then
    echo "Width of a restraint-free region for flat-bottom harmonic restraints should be positive!"
    exit -1
fi
          
inp=$1

hostname=`uname -n`
ttag=$hostname-$$
tag="${ttag,,}"

pwd=`pwd`

/bin/rm -rf $tmpdir/$tag
mkdir -p $tmpdir/$tag
cd $tmpdir/$tag

# initial prep
/bin/cp $pwd/$inp init0.pdb
convpdb.pl -out generic -setchain "A" -setall -cleanaux -nsel protein -orient init0.pdb > init.pdb
export CHARMMEXEC=$CHARMMEXECSINGLE
locprefmd.sh init.pdb > iniref.pdb
convpdb.pl -center iniref.pdb > iniref.center.pdb
convpdb.pl -segnames iniref.center.pdb > iniref.seg.pdb
convpdb.pl -nsel CA -readseg iniref.seg.pdb > caref.pdb

# equilibration
export CHARMMEXEC="mpirun -np $cpus $CHARMMEXECMPI"
equiCHARMM.pl -par param=22x,xpar=$parfile,xtop=$topfile -neutralize -cutoff 9 -cons cab iniref.seg.pdb 0:9999_0.5 iniref.center.pdb

grep CRYSTAL -A 2 md.equi.3.restart | tail -n+2 | sed -e "s/D/E/g" | awk '{printf " %s %s %s", $1, $2, $3}' | awk '{printf "X= %f\nY= %f\nZ= %f\n", $1, $3, $6}' > boxsize

#boxsize=`cat boxsize`
boxx=$(grep X boxsize | awk '{print $2}')
boxy=$(grep Y boxsize | awk '{print $2}')
boxz=$(grep Z boxsize | awk '{print $2}')
boxshape='ortho'
boxsize="boxx=$boxx,boxy=$boxy,boxz=$boxz"

# to generate PDB files for OpenMM runs
trX=$(echo $boxx/2 | bc -l)
trY=$(echo $boxy/2 | bc -l)
trZ=$(echo $boxz/2 | bc -l)
mv iniref.seg.pdb iniref.seg.0.pdb
convpdb.pl -readseg -segnames -translate $trX $trY $trZ iniref.seg.0.pdb > iniref.seg.pdb
convpdb.pl -readseg -segnames -translate $trX $trY $trZ md.equi.3.pdb > md.prod.init.pdb

genPSF.pl -par param=22x,xpar=$parfile,xtop=$topfile -crdout md.prod.crd md.prod.init.pdb > md.prod.psf
awk 'BEGIN {show=1;} /MOLNT/ {show=0;} /NCRTERM/ {show=1;} show == 1 {print}' md.prod.psf > md.prod.openmm.psf

par="param=22x,xpar=$parfile,xtop=$topfile,dyntstep=0.002,dynsteps=$mdsteps,dynoutfrq=5000,dyntemp=298,lang=1,langfbeta=0.01,langupd=0,openmm,dyneqfrq=0,boxshape=$boxshape,$boxsize,periodic,dyntrfrq=0"

if [[ $bottom > 0.0 ]]; then
    cons=`echo $kcons $bottom | awk '{printf("ca iniref.seg.pdb 0:9999_%1.5f_%02.5f", $1, $2);}'`
else
    cons=`echo $kcons | awk '{printf("ca iniref.seg.pdb 0:9999_%1.5f",$1);}'`
fi

convpdb.pl -nsel protein -setchain "A" -setall md.prod.init.pdb > solute.pdb
convpdb.pl -nsel CA solute.pdb > solute.ca.pdb
nsol=`grep ATOM solute.pdb | wc -l`
awk '/ATOM/ {print $2}' solute.ca.pdb > ca.list

for run in $(seq 1 $mdruns); do
  mdOpenMM.py -par $par -cons $cons -restout md.$run.restart -trajout md.$run.dcd -log md.$run.log -final md.$run.pdb -psf md.prod.openmm.psf md.prod.crd md.prod.init.pdb
  mdconv -out solute.$run.dcd -atoms 1:$nsol md.$run.dcd 
  mdconv -atomlist ca.list -out solca.$run.dcd solute.$run.dcd 
done

if [ -r solute.10.dcd ]; then
  mdconv -out solute.dcd solute.?.dcd solute.??.dcd
  mdconv -out solca.dcd solca.?.dcd solca.??.dcd
else
  mdconv -out solute.dcd solute.?.dcd
  mdconv -out solca.dcd solca.?.dcd
fi

processDCD.pl -step 10 -ensdir ens -ens sample solute.pdb solute.dcd
gunzip ens/*/*/sample.pdb.gz

# calculating iRMSD/iGDT
tpwd=`pwd`
ensrun.pl -overwrite -cpus $cpus -update 200 -nocompress -noinp -dir ens -set irmsd0:2,igdtha0:5 sample tmscore.pl sample.pdb $tpwd/iniref.seg.pdb
ensrun.pl -overwrite -cpus $cpus -update 200 -nocompress -noinp -dir ens -set rwplus sample rwplus.sh sample.pdb

# generate average structures
nstruct=`ensfiles.pl -dir ens sample | wc -l`
skip=`echo $nstruct | awk '{print int($1/10+0.5)}'`
#echo nstruct $nstruct skip $skip
averagerefine.pl ens irmsd0 rwplus $skip

convpdb.pl -center average.pdb | convpdb.pl -segnames > avg.center.pdb

export CHARMMEXEC="mpirun -np $cpus $CHARMMEXECMPI"
equiCHARMM.pl -par param=22x,xpar=$parfile,xtop=$topfile -neutralize -cutoff 9 -prefix avg -densitysteps 2000 -equi 50:1000=100:1000=200:1000=250:1000=298:5000 -cons ca avg.center.pdb 0:9999_10 avg.center.pdb

export CHARMMEXEC=$CHARMMEXECSINGLE
convpdb.pl -nsel protein -setchain " " -setall avg.heat.298.pdb > avg.protein.pdb
locprefmd.sh avg.protein.pdb > avg.clean.pdb
locprefmd.sh avg.clean.pdb > avg.clean2.pdb

convpdb.pl -setchain " " -setall -out generic avg.clean2.pdb 
cd $pwd

/bin/rm -rf $tmpdir/$tag
