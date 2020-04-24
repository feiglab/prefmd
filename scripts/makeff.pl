#!/usr/bin/env perl

# this utility combines CHARMM force field files for proteins and water into a single file
#
# input (in the current directory):
#    top_all36_prot.rtf
#    par_all36_prot.rtf
#    toppar_water_ions.str
#
# output:
#    top.rtf
#    par.prm
#

use strict;

my $topin="top_all36_prot.rtf";
my $parin="par_all36_prot.prm";

my $wstream="toppar_water_ions.str";

my $topout="top.rtf";
my $parout="par.prm";

my $readtop=0;
my @wmass=();
my @rtop=();
my %par;
my $section="none";
my $lastsection="none";
open INP,"$wstream";
while (<INP>) {
  if (/^MASS/ && !$readtop) {
    push(@wmass,$_);
  } elsif (/default first/) {
    $readtop=1;
  } elsif (/END/ && $readtop == 1) {
    $readtop=2;
  } elsif ($readtop == 1) {
    push(@rtop,$_);
  } elsif ($readtop == 2) {
    if (/^ATOM/) {
      $section="atoms";
    } elsif (/^BOND/) {
      $section="bonds";
    } elsif (/^ANGLE/) {
      $section="angles";
    } elsif (/^DIHEDRAL/) {
      $section="dihedrals";
    } elsif (/^IMPROPER/) {
      $section="impropers";
    } elsif (/^CMAP/) {
      $section="cmap";
    } elsif (/^NONBON/) {
      $section="nonbonded";
      if (/ -\s*$/) {
        <INP>;
      }
    } elsif (/^NBFIX/) {
      $section="nbfix";
    } elsif (/^HBOND/) {
      $section="hbond";
    } elsif (/^END/) {
      $section="done";
    }

    if ($section ne $lastsection) {
      $par{$section}=() if (!exists $par{$section});
    } else {
      if (!/^\!/ && $section ne "none" && $section ne "done") {
        push(@{$par{$section}},$_);
      }
    }

    $lastsection=$section;
  }
}
close INP;

open OUT,">$topout";
open INP,"$topin";
my @tmass=();
my $topmassout=0;
while (<INP>) {
  if (/^MASS/) {
    push(@tmass,$_);
  } elsif (/^DECL/ && !$topmassout) {
    foreach my $t (@tmass) {
      print OUT $t;
    }
    printf OUT "\n"; 
    foreach my $t (@wmass) {
      print OUT $t;
    } 
    printf OUT "\n"; 
    print OUT $_;
    $topmassout=1; 
  } elsif (/^END/) {
    foreach my $t (@rtop) {
      print OUT $t;
    }
    print OUT $_;
  } else {
    print OUT $_;
  }
}
close INP;
close OUT;

open OUT,">$parout";
open INP,"$parin";
$section="none";
$lastsection="none";
my %have;

foreach my $a ( @{$par{atoms}} ) {
  if ($a=~/MASS\s+\S+\s+(\S+)\s+/) {
    $have{$1}=1;
  }
}

while (<INP>) {
  my $line=$_;
  if (/^ATOM/) {
    $section="atoms";
  } elsif (/^BOND/) {
    $section="bonds";
  } elsif (/^ANGLE/) {
    $section="angles";
  } elsif (/^DIHEDRAL/) {
    $section="dihedrals";
  } elsif (/^IMPROPER/) {
    $section="impropers";
  } elsif (/^CMAP/) {
    $section="cmap";
  } elsif (/^NONBON/) {
    $section="nonbonded";
  } elsif (/^NBFIX/) {
    $section="nbfix";
  } elsif (/^HBOND/) {
    $section="hbond";
  } elsif (/^END/) {
    $section="done";
  }

  if ($section eq "atoms" && /MASS\s+\S+\s+(\S+)\s+/) {
    $have{$1}=1;
  }

  if ($lastsection ne $section) {
    if (exists $par{$lastsection} && $#{$par{$lastsection}}>=0) {
      foreach my $t ( @{$par{$lastsection}} ) {
        if ($lastsection eq "nbfix") {
          if ($t=~/^\s*(\S+)\s*(\S+)\s*/) {
            if ($have{$1} && $have{$2}) {
              print OUT $t;
            }
          }
        } else {
          print OUT $t;
        }
      }
    } 

    if ($lastsection eq "nonbonded" && ($section eq "hbond" || $section eq "done")) {
      if (exists $par{nbfix} && $#{$par{nbfix}}>=0) {
        print OUT "\n";
        print OUT "NBFIX\n";
        foreach my $t ( @{$par{nbfix}} ) {
          if ($t=~/^\s*(\S+)\s*(\S+)\s*/) {
            if ($have{$1} && $have{$2}) {
              print OUT $t;
            }
          }
        }
        print OUT "\n";
      }
    }
  } 

  $lastsection=$section;
   
  print OUT $line;
}
close INP;
close OUT;

