#!/usr/bin/perl
#
use strict;
use Math::Trig;
use POSIX;

my $ens=shift @ARGV;
my $prop1=shift @ARGV;
my $prop2=shift @ARGV;
my $skipinitial=shift @ARGV;

$ens="ens" if (!defined $ens);
$prop1="irmsd0" if (!defined $prop1);
$prop2="rwplus" if (!defined $prop2);
$skipinitial=10 if (!defined $skipinitial);

my $tavgirmsd=0;
my $tavgirmsd2=0;
my $tavgrwplus=0;
my $tavgrwplus2=0;
my $navg;

#printf "skipinitial %d\n",$skipinitial;

my @data=();
open INP,"getprop.pl -dir $ens -prop $prop1,$prop2 sample |";
while (<INP>) {
  chomp;
  my @f=split(/\s+/);

  if ($f[0]>$skipinitial) {
    my $trec={};
    $trec->{inx}=$f[0];
    $trec->{irmsd}=$f[1];
    $trec->{rwplus}=$f[2];
    $trec->{fname}=sprintf("%s/%d/%d/sample.pdb",$ens,int($trec->{inx}/100),$trec->{inx}-int($trec->{inx}/100)*100);
#    printf "selected %s\n",$trec->{fname};
    push(@data,$trec); 

    $tavgirmsd+=$trec->{irmsd};
    $tavgirmsd2+=$trec->{irmsd}*$trec->{irmsd};
    $tavgrwplus+=$trec->{rwplus};
    $tavgrwplus2+=$trec->{rwplus}*$trec->{rwplus};
    $navg++;
  }
}
close INP;

$tavgirmsd/=$navg;
$tavgrwplus/=$navg;

my $tval=($tavgirmsd2-$navg*($tavgirmsd*$tavgirmsd))/($navg-1);
my $sdevirmsd=($tval<0)?0:sqrt($tval);

$tval=($tavgrwplus2-$navg*($tavgrwplus*$tavgrwplus))/($navg-1);
my $sdevrwplus=($tval<0)?0:sqrt($tval);

my $avgirmsd=0;
my $avgirmsd2=0;
my $avgrwplus=0;
my $avgrwplus2=0;
$navg=0;

foreach my $d ( @data ) {
  if ($d->{rwplus}<$tavgrwplus+3.0*$sdevrwplus && $d->{irmsd}<$tavgirmsd+3.0*$sdevirmsd) {
    $avgirmsd+=$d->{irmsd};
    $avgirmsd2+=$d->{irmsd}*$d->{irmsd};
    $avgrwplus+=$d->{rwplus};
    $avgrwplus2+=$d->{rwplus}*$d->{rwplus};
    $navg++;
  }
}

$avgirmsd/=$navg;
$avgrwplus/=$navg;

#printf "navg: %d\n",$navg;

my $tval=($avgirmsd2-$navg*($avgirmsd*$avgirmsd))/($navg-1);
my $sdevirmsd=($tval<0)?0:sqrt($tval);

$tval=($avgrwplus2-$navg*($avgrwplus*$avgrwplus))/($navg-1);
my $sdevrwplus=($tval<0)?0:sqrt($tval);

my $rad=1;
my $gamma=35;
my $theta=240;

my $qx = cos(deg2rad($theta ));
my $qy = sin(deg2rad($theta ));
my $qz = 0.0;
my $qsize=sqrt($qx*$qx + $qy*$qy + $qz*$qz);

my $naccept=0;
open OUT,"| pdb2traj.pl -out subset.dcd -f -";
foreach my $d ( @data ) {
  my $s1=($d->{irmsd}-$avgirmsd)/$sdevirmsd;
  my $s2=($d->{rwplus}-$avgrwplus)/$sdevrwplus;

  my $accept = 0;
  if($s1*$s1+$s2*$s2 >= $rad*$rad &&
    (acos(($s1*$qx + $s2*$qy) / (sqrt($qsize*($s1*$s1 + $s2*$s2)))) * 180.0 / 3.1415) < $gamma) {
    $d->{select}=1;
    printf OUT "%s\n",$d->{fname};
    $naccept++;
  } else {
    $d->{select}=0;
  }
}  
close OUT;

#printf "naccept: %d\n",$naccept;

system("analyzeCHARMM.pl -avg -fitsel CA -comp iniref.seg.pdb -pdb solute.pdb subset.dcd | convpdb.pl -setchain ' ' -setall > average.pdb");

