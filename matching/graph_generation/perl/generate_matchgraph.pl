#!/usr/bin/env perl
##############################################################################
# SCRIPT NAME:	load_wmda.pl
# DESCRIPTION:	load WMDA imputation results to graph csv for import
# PARAMETERS:	none
# OUTPUT:	
# TABLES:       
#
# DATE WRITTEN: 2015-07-22
# WRITTEN BY:   Martin Maiers
#
# REVISION HISTORY: 
# REVISION DATE		REVISED BY	DESCRIPTION 
# ------- ----------	--------------	-------------------------------------
# 1.0     2015-07-22    M. Maiers       Orig
# 1.1     2016-02-04    M. Maiers       SLGs
# 1.2     2016-09-05    M. Maiers       MUG-1s
# 1.3     2017-08-10    M. Maiers       ALLELEs
#
#
#       COPYRIGHT (C) 2015 NATIONAL MARROW DONOR PROGRAM.  
#               ALL RIGHTS RESERVED        
##############################################################################
use strict;    # always
use warnings;  # or else
use REST::Client;
use JSON;
use Data::Dumper;
use MIME::Base64;

# nodes
my %S;  # hash of Subjects id_typ => id
my %G;  # hash of Genotypes attributes: genotypeId, gf
my %H;  # hash of Haplotypes
my %MUUG;   # hash of MUUGs attributes: mass, homloci
my %MUUG_1; # hash of MUUG_1s attributes: 
my %SLG;    # hash of SLGs attributes: locus, mass, hom
my %ALLELE; # hash of SLGs attributes: locus, mass, hom

# edges
my %GL; # Genotype likelihoods: link between subject to Genotype with prob
my %ML; # MUUG likelihoods: link between subject to MUUG with property prob
my %SL; # SLG likelihoods: link between subject to SLG with property prob
my %M1; # M1: link between MUUG and MUUG_1 
my %SA; # SA: link between SLG and ALLELE 
my %CT; # hash of COLLAPSES_TO
my %HH; # hash of HAS_HAPLOTYPE

my $v =1; # verbose

# dump ids with a null allele for further review
#open NULLD, ">did.N" or die "$!: did.N";


# storage
my $graphdir = "./output/graph";

# nodes filesnames
my $subject_file   = "$graphdir/Subject.csv";
my $haplotype_file = "$graphdir/Haplotype.csv";
my $genotype_file  = "$graphdir/Genotype.csv";
my $muug_file      = "$graphdir/MUUG.csv";
my $muug_1_file    = "$graphdir/MUUG_1.csv";
my $slg_file       = "$graphdir/SLG.csv";
my $allele_file    = "$graphdir/ALLELE.csv";

# edges filenames
my $ml_file        = "$graphdir/MUUG_LIKELIHOOD.csv";
my $sl_file        = "$graphdir/SLG_LIKELIHOOD.csv";
my $gl_file        = "$graphdir/GENOTYPE_LIKELIHOOD.csv";
my $ct_file        = "$graphdir/COLLAPSES_TO.csv";
my $hh_file        = "$graphdir/HAS_HAPLOTYPE.csv";
my $m1_file        = "$graphdir/M1.csv";
my $sa_file        = "$graphdir/SA.csv";


# open node files
open SUBJECT_FILE,   ">$subject_file" or die "$!: $subject_file";
open HAPLOTYPE_FILE, ">$haplotype_file" or die "$!: $haplotype_file";
open GENOTYPE_FILE,  ">$genotype_file" or die "$!: $genotype_file";
open MUUG_FILE,      ">$muug_file" or die "$!: $muug_file";
open MUUG_1_FILE,    ">$muug_1_file" or die "$!: $muug_1_file";
open SLG_FILE,       ">$slg_file" or die "$!: $slg_file";
open ALLELE_FILE,    ">$allele_file" or die "$!: $allele_file";

# open edge files
open ML_FILE,        ">$ml_file" or die "$!: $ml_file";
open SL_FILE,        ">$sl_file" or die "$!: $sl_file";
open GL_FILE,        ">$gl_file" or die "$!: $gl_file";
open CT_FILE,        ">$ct_file" or die "$!: $ct_file";
open HH_FILE,        ">$hh_file" or die "$!: $hh_file";
open M1_FILE,        ">$m1_file" or die "$!: $m1_file";
open SA_FILE,        ">$sa_file" or die "$!: $sa_file";

# headers
print SUBJECT_FILE "subjectId:ID(Subject),id,id_typ,:LABEL\n";
my $subjectId=0;
print HAPLOTYPE_FILE "haplotypeId:ID(Haplotype),name,:LABEL\n";
my $haplotypeId=0;
print GENOTYPE_FILE "genotypeId:ID(Genotype),name,mass:double,:LABEL\n";
my $genotypeId=0;
print MUUG_FILE "muugId:ID(MUUG),name,mass:double,homloci,:LABEL\n";
my $muugId=0;
print SLG_FILE "slgId:ID(SLG),locus,name,mass:double,hom,:LABEL\n";
my $slgId=0;
print MUUG_1_FILE "muug_1Id:ID(MUUG_1),name,:LABEL\n";
my $muug_1Id=0;
print ALLELE_FILE  "alleleId:ID(ALLELE),name,:LABEL\n";
my $alleleId=0;


print ML_FILE ":START_ID(Subject),prob:double,:END_ID(MUUG)\n";
print SL_FILE ":START_ID(Subject),prob:double,:END_ID(SLG)\n";
print GL_FILE ":START_ID(Subject),prob,:END_ID(Genotype)\n";
print CT_FILE ":START_ID(Genotype),:END_ID(MUUG)\n";
print HH_FILE ":START_ID(Genotype),:END_ID(Haplotype)\n";
print M1_FILE ":START_ID(MUUG),:END_ID(MUUG_1)\n";
print SA_FILE ":START_ID(SLG),:END_ID(ALLELE)\n";


# imputation results using same HF set as matching 
my $don_file = "../../../multi_race_impute/output/don.gl.txt_out";
my $pat_file = "../../../multi_race_impute/output/pat.gl.txt_out";
my $file = "cat $don_file $pat_file |";
open FILE, $file or die "$!: $file";

# parse the imputation output of the subjects
print STDERR "parse imputation output\n" if $v;
while(<FILE>) {
  chomp;
  my ($id, $hap1, $freq1, $pop1, $hap2, $freq2, $pop2) = split /,/;
  my $g = join ('+', sort $hap1, $hap2);
  my $id_typ = substr($id, 0, 1);
  my $idnum = int (substr ($id, 1));

  #uncomment for ids up to 30
  #next if $idnum >30;

  # don't use $gf from imputation since it doesn't do HW
  my $hwgf=0;
  if($hap1 eq $hap2) { $hwgf = $freq1 * $freq2; } else { $hwgf = 2 * $freq1 * $freq2;}

  # 
  # generate hash for  subject and genotype
  # 
  countSubject($id, $id_typ, $g, $hwgf);
}

# 
# output subject nodes
# 
print STDERR "output Subject nodes\n" if $v;
foreach my $id (keys %S) {
  my $id_typ = $S{$id}{id_typ};
  my $subjectId = $S{$id}{subjectId};
  # Subject.csv
  print SUBJECT_FILE join ',', $subjectId, $id, $id_typ, "Subject"; 
  print SUBJECT_FILE  "\n";
}

# 
# output haplotype nodes
# 
print STDERR "output Haplotype nodes\n" if $v;
foreach my $h (keys %H) {
  # Haplotype.csv
  print HAPLOTYPE_FILE join ',', $H{$h}{haplotypeId}, $h, "Haplotype";
  print HAPLOTYPE_FILE "\n";
}


#
# compute and then output normalized genotype likelihoods
#
print STDERR "compute and output GL edges\n" if $v;
foreach my $id (keys %GL) {
  my $gf_total=0;
  foreach my $g (keys %{$GL{$id}}) {
    $gf_total+=$GL{$id}{$g}{gf};
  }
  foreach my $g (keys %{$GL{$id}}) {
    # GENOTYPE_LIKELIHOOD.csv
    # subject, prob, genotypeId
    my $ngf = $gf_total ? $GL{$id}{$g}{gf}/$gf_total : 0;
    print GL_FILE join ',', $id, $ngf, $G{$g}{genotypeId};
    print GL_FILE "\n";
    $G{$g}{mass}+= $ngf;
  }
}
# 
# output genotype nodes
# 
print STDERR "compute and output Genotype nodes\n" if $v;
print STDERR "compute and output Haplotype nodes\n" if $v;
foreach my $g (keys %G) {
  # Genotype.csv
  print GENOTYPE_FILE join ',', $G{$g}{genotypeId}, $g, $G{$g}{mass}, "Genotype";
  print GENOTYPE_FILE "\n";
  my ($h1, $h2) = split /\+/, $g;
  print HH_FILE join (',', $G{$g}{genotypeId}, $H{$h1}{haplotypeId}), "\n";
  print HH_FILE join (',', $G{$g}{genotypeId}, $H{$h2}{haplotypeId}), "\n";
}


#
# compute and then output normalized muug likelihoods
#
print STDERR "compute and output ML edges\n" if $v;
foreach my $id (keys %ML) {
  my $mf_total=0;
  # $ML{$id}{$muug}{mf}+= $gf; # absolute gf adds to MUUG frequency 
  foreach my $muug (keys %{$ML{$id}}) {
    $mf_total+=$ML{$id}{$muug}{mf};
  }

  # MUUGs

  foreach my $muug (keys %{$ML{$id}}) {
    # MUUG_LIKELIHOOD.csv
    # subject, nmf, muugId
    my $nmf = $mf_total ? $ML{$id}{$muug}{mf}/$mf_total : 0;
    print ML_FILE join ',', $id, $nmf, $MUUG{$muug}{muugId};
    print ML_FILE "\n";
    $MUUG{$muug}{mass}+= $nmf;

    # SLGs
    foreach my $slg (split /\^/, $muug) {
      $SLG{$slg}{mass} += $nmf;
      $SLG{$slg}{locus} = (split /\*/, $slg)[0];
      $SLG{$slg}{hom}   = getHom($slg);
      $SL{$id}{$slg}{prob}+= $nmf; #
    }

    # MUG_1
    my @g_new;
    # convert muug to a matrix
    my @loci = split /\^/, $muug;
    for (my $i =0; $i<=$#loci; $i++) {
      my @genotype = split /\+/, $loci[$i]; 
      for (my $j =0; $j<=$#genotype; $j++) {
        $g_new[$i][$j] = $genotype[$j];
      }
    }
    # exchange each position for a wildcard and register the MUUG_1
    for (my $i =0; $i<=$#loci; $i++) {
      my @genotype = split /\+/, $loci[$i]; 
      for (my $j =0; $j<=$#genotype; $j++) {
        my $a_orig = $g_new[$i][$j];
        $g_new[$i][$j] = "X";
        
        my @m_1;
        for (my $i =0; $i<=$#loci; $i++) {
          push @m_1, join ('+', sort @{$g_new[$i]});
        }
        my $muug_1 = join ('^', @m_1);
        # create MUUG_1 node
        if (!defined $MUUG_1{$muug_1}) {
          $muug_1Id++; #accession a new muug_1Id
          $MUUG_1{$muug_1}{muug_1Id}=$muug_1Id;
        } 
        # create link between MUUG and MUUG_1
        $M1{$muug}{$muug_1}++;
        $g_new[$i][$j] = $a_orig;
      }
    }
  }
}

print STDERR "compute and output SLG edges\n" if $v;
foreach my $id (keys %SL) {
  foreach my $slg (keys %{$SL{$id}}) {
    print SL_FILE join ',', $id, $SL{$id}{$slg}{prob}, $SLG{$slg}{slgId};
    print SL_FILE "\n";
  }
}

print STDERR "output M1 edges\n" if $v;
foreach my $muug (keys %M1) {
  foreach my $muug_1 (keys %{$M1{$muug}}) {
    print M1_FILE join ',', $MUUG{$muug}{muugId}, $MUUG_1{$muug_1}{muug_1Id};
    print M1_FILE "\n";
  }
}

# generate ALLELE nodes and SA SLG-ALLELE edges
foreach my $slg (keys %SLG) {
  my ($a1, $a2) = split /\+/, $slg;
  if (!defined $ALLELE{$a1}) {
    $alleleId++; #accession a new alleleId
    $ALLELE{$a1}{alleleId}=$alleleId;
  } 
  if (!defined $ALLELE{$a2}) {
    $alleleId++; #accession a new alleleId
    $ALLELE{$a2}{alleleId}=$alleleId;
  } 
  # create link between SLG and ALLELE
  $SA{$slg}{$a1}++;
  $SA{$slg}{$a2}++;
}
print STDERR "output SA edges\n" if $v;
foreach my $slg (keys %SA) {
  foreach my $allele (keys %{$SA{$slg}}) {
    print SA_FILE join ',', $SLG{$slg}{slgId}, $ALLELE{$allele}{alleleId};
    print SA_FILE "\n";
  }
}

#
# output collapses_to links and set MUUG attributes
#
print STDERR "output CT edges\n" if $v;
foreach my $muug (keys %CT) {
  # each muug has its intrinsic GF (mass)
  my $muug_gf = 0;

  # determine the number of loci (e.g. 5)
  my $numloci =  scalar (split /\^/, $muug);

  my $homloci = getHomLoci($muug);
  $MUUG{$muug}{homloci} = $homloci;

  foreach my $g (keys %{$CT{$muug}}) {
    print CT_FILE join ',', $G{$g}{genotypeId}, $MUUG{$muug}{muugId};
    print CT_FILE "\n";
  }
}
# 
# output MUUG nodes
# 
print STDERR "output MUUG nodes\n" if $v;
foreach my $m (keys %MUUG) {
  # MUUG.csv
  print MUUG_FILE join ',', $MUUG{$m}{muugId}, $m, $MUUG{$m}{mass}, $MUUG{$m}{homloci}, "MUUG";
  print MUUG_FILE "\n";
}

# 
# output SLG nodes
# 
print STDERR "output SLG nodes\n" if $v;
foreach my $s (keys %SLG) {
  # SLG.csv
  print SLG_FILE join ',', $SLG{$s}{slgId}, $SLG{$s}{locus}, $s, $SLG{$s}{mass}, $SLG{$s}{hom}, "SLG";
  print SLG_FILE "\n";
}

# 
# output MUUG_1 nodes
# 
print STDERR "output MUUG_1 nodes\n" if $v;
foreach my $s (keys %MUUG_1) {
  # MUUG_1.csv
  print MUUG_1_FILE join ',', $MUUG_1{$s}{muug_1Id}, $s, "MUUG_1";
  print MUUG_1_FILE "\n";
}

# 
# output ALLELE nodes
# 
print STDERR "output ALLELE nodes\n" if $v;
foreach my $s (keys %ALLELE) {
  # ALLELE.csv
  print ALLELE_FILE join ',', $ALLELE{$s}{alleleId}, $s, "ALLELE";
  print ALLELE_FILE "\n";
}


exit 0;

sub countSubject {
  my ($id, $id_typ, $g, $gf) = @_;

  if (!defined $S{$id}) { $S{$id}{subjectId}=++$subjectId; }
  $S{$id}{id_typ} = $id_typ;


  #print STDERR "countSubject: $id $id_typ $g\n" if $v;
  my ($h1, $h2) = split /\+/, $g;

  # convert haplotypes to array of loci
  my @l1 = split /\~/, $h1;
  my @l2 = split /\~/, $h2;
  my @muug;
  my @slg;
  my $hasnull=0;

  # for each locus, build up the mug by adding an allele pair
  for (my $i=0; $i<=$#l1; $i++) {
    # 
    # 
    # this is where we move from genotype from phenotype
    # due to superstition and weak-minded fear of the letter N
    # lacking any clinical evidence
    #
    my $mom = $l1[$i];
    my $dad = $l2[$i];
    my $momnull++ if $mom=~/N$/;
    my $dadnull++ if $dad=~/N$/;

    if ($momnull && !$dadnull) {
      $mom = $dad;
    } elsif ($dadnull && !$momnull) {
      $dad = $mom;
    }
    push @muug, join ('+', sort $mom, $dad);
    $hasnull++ if $l1[$i]=~/N$/;
    $hasnull++ if $l2[$i]=~/N$/;
  }
  my $muug = join '^', @muug;
  #print STDERR join ('	', $g, $h1, $h2, $muug), "\n" if $v;
  #print NULLD join ('	', $id, $id_typ, $g), "\n" if $hasnull;

  # ML
  $ML{$subjectId}{$muug}{mf}+= $gf; # absolute gf adds to MUUG frequency 
  $GL{$subjectId}{$g}{gf}= $gf; # absolute gf adds to gf 

  if (!defined $H{$h1}) { $H{$h1}{haplotypeId}=++$haplotypeId; }
  if (!defined $H{$h2}) { $H{$h2}{haplotypeId}=++$haplotypeId; }

  if (!defined $G{$g}) { $G{$g}{genotypeId}=++$genotypeId; }

  if (!defined $MUUG{$muug}) {
    $muugId++; #accession a new muugId
    $MUUG{$muug}{muugId}=$muugId;
  }

  foreach my $slg (@muug) {
    if (!defined $SLG{$slg}) {
      $slgId++; #accession a new slgId
      $SLG{$slg}{slgId}=$slgId;
    }
  }

  # COLLAPSES_TO relationship
  $CT{$muug}{$g}++;
}

sub getHomLoci {
  my $m = shift;
  my $ret = 0;
  foreach my $loc (split /\^/, $m) {
    $ret++ if getHom($loc);
  }
  return $ret;
}

sub getHom {
  my $loc = shift;
  my $ret = 0;
  my ($g1, $g2) = split /\+/, $loc;
  $ret++ if $g1 eq $g2;
  return $ret;
}
__END__
