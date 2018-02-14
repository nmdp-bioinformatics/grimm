#!/usr/bin/perl
##############################################################################
# SCRIPT NAME:	g2gl.pl
# DESCRIPTION:	genotype to genotype-list
#
# DATE WRITTEN: 2017-02-08
# WRITTEN BY:   Martin Maiers
#
##############################################################################
use strict;    # always
use warnings;  # or else
use MAC;
use Getopt::Std;

my %opts;
getopts('i:o:', \%opts);
die "$0: -i inputfile -o outputfile" unless defined $opts{i} && defined $opts{o};



# set to absolute path of your installation
my $top = "../..";

# set to HF file
my $hf_file = "$top/graph_generator/data/wmda/freqs.txt";

# set to input genotypes file
my $gfile = $opts{i};
my $ofile = $opts{o};

##############################################################################
# parse file of haplotype frequencies
##############################################################################
open HF_FILE, $hf_file or die "$!: $hf_file";
my %A;  # hash table of alleles
while(<HF_FILE>) {
  chomp;
  #A*01:01g~C*01:02g~B*15:01g~DRB1*01:01~DQB1*05:01;6e-05
  my ($hap, $freq) = split /\;/;
  foreach my $allele (split /\~/, $hap) {
    $allele=~s/g//;      # remove little-g
    my $who_allele = "HLA-".$allele;
    $A{$who_allele}+=$freq;  # note that this allele exists; sum freq just because
  }
}
close HF_FILE;

##############################################################################
# parse input file
##############################################################################
open OFILE, ">$ofile" or die "$!: $ofile";
open GFILE, $gfile or die "$!: $gfile";
while(<GFILE>) {
  chomp;
  #D000001%A*11:XX+A*02:XX^C*UUUU+C*UUUU^B*13:XX+B*44:XX^DRB1*07:XX+DRB1*01:XX^DQB1*UUUU+DQB1*UUUU
  my ($id, $g) = split /\%/;
  my @gl_array = ();
  foreach my $loctyp (split /\^/, $g) {
    my ($t1, $t2) = split /\+/, $loctyp;
    my ($l1, $a1) = split /\*/, $t1;
    my ($l2, $a2) = split /\*/, $t2;

    if ($l1 ne $l2) {
      warn "mismatching loci: $l1 vs $l2 for id: $id";
      next;
    }

    if ($a1 eq "UUUU" || $a2 eq "UUUU") {
      # if either typing is UUUU then skip this loctyp
      next; 
    }

    my $who_typ1 = "HLA-".$t1;
    my $who_typ2 = "HLA-".$t2;

    # expand
    my @al1 = expandTyp($id, $t1);
    my @al2 = expandTyp($id, $t2);
  
    # reduce
    my $gl = join ('+', 
    sort  join ('/', reduceTyp(@al1)), join ('/', reduceTyp(@al2)));
      
    $gl=~s/HLA\-//g;
    push @gl_array, $gl;
    }
  print OFILE join ('%', $id, join ('^', @gl_array)), "\n";
}
close OFILE;
close GFILE;

exit 0;


sub reduceTyp {
  my (@al) = @_;
  my @ret;
  
  foreach (@_) {
    push @ret, $_ if defined $A{$_};
  }
  if (@ret) {
    return @ret;
  } else {
    return "";
  }
}

sub expandTyp {
  my ($id, $t) = @_;
  my $who_typ = "HLA-".$t;
  my @al;

  # expand
  my $rc = MAC::decode($who_typ, \@al);
  if ($rc ne 200) {
    warn "decode failed with code: $rc for $who_typ for id: $id";
    next;
  }
  return @al;
}
