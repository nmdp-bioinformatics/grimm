#!/usr/bin/env perl
##############################################################################
# SCRIPT NAME:	compare_results.pl
# DESCRIPTION:	
#
# DATE WRITTEN: 2017-08-11
# WRITTEN BY:   Martin Maiers
#
##############################################################################
use strict;    # always
use warnings;  # or else

my $outfile = "output/cmp.txt";
open OUTFILE, ">$outfile" or die "$!: $outfile";


my $file= "../../graph_generator/data/wmda/set3.consensus.txt";
open FILE, $file or die "$!: $file";
#open FILE, "head -100000 $file|" or die "$!: $file";

my %C;
while(<FILE>) {
  chomp;
  my (@data) = split /\;/;
  my ($r, $d, $a_mg, $a_p, $c_mg, $c_p, $b_mg, $b_p,
    $drb1_mg, $drb1_p, $dqb1_mg, $dqb1_p, $p9of10, $p10of10) = @data;
  $C{$r}{$d} = $_;
}
close FILE;

$file= "output/mr.txt";
open FILE, $file or die "$!: $file";
#open FILE, "head -100000 $file|" or die "$!: $file";

my %N;
while(<FILE>) {
  chomp;
  my (@data) = split /\;/;
  my ($r, $d, $a_mg, $a_p, $c_mg, $c_p, $b_mg, $b_p,
    $drb1_mg, $drb1_p, $dqb1_mg, $dqb1_p, $p9of10, $p10of10) = @data;
  $N{$r}{$d} = $_;
}
close FILE;

my @loci = qw/A B DRB1 DQB1 C/;

foreach my $r (sort keys %C) {
  foreach my $d (sort keys %{$C{$r}}) {
    if (not exists $N{$r}{$d}) {
       print OUTFILE join("\t", "MISSING", $r, $d), "\n";
       next;
    }
    my $crd = parse_line($C{$r}{$d});
    my $nrd = parse_line($N{$r}{$d});
    if ($crd->{p10} != $nrd->{p10}) {
      print OUTFILE join ("\t", "DIFF10", $r, $d,
            $crd->{p10}, $nrd->{p10}, abs($crd->{p10} - $nrd->{p10})), "\n";
    }
    if ($crd->{p9} != $nrd->{p9}) {
      print OUTFILE join ("\t", "DIFF9", $r, $d,
            $crd->{p9}, $nrd->{p9}, abs($crd->{p9} - $nrd->{p9})), "\n";
    }
    foreach my $loc (@loci) {
      my $mg = $loc."_mg";
      my $p  = $loc."_p";
      if ($crd->{$mg} ne $nrd->{$mg}) {
        print OUTFILE join ("\t", "DIFF".$loc."_MG", $r, $d,
              $crd->{$mg}, $nrd->{$mg}), "\n";
      }
      if ($crd->{$p} != $nrd->{$p}) {
        print OUTFILE join ("\t", "DIFF".$loc."_P", $r, $d,
              $crd->{$p}, $nrd->{$p}, abs($crd->{$p} - $nrd->{$p})), "\n";
      }
    }
  }
}
close OUTFILE;
exit 0;

sub parse_line {
  my $line = shift;
  my (@data) = split /\;/, $line;
  my ($r, $d, $a_mg, $a_p, $c_mg, $c_p, $b_mg, $b_p,
    $drb1_mg, $drb1_p, $dqb1_mg, $dqb1_p, $p9of10, $p10of10) = @data;

  my $data = {};
  $data->{A_mg} = $a_mg;
  $data->{A_p} = $a_p;
  $data->{B_mg} = $b_mg;
  $data->{B_p} = $b_p;
  $data->{C_mg} = $c_mg;
  $data->{C_p} = $c_p;
  $data->{DQB1_mg} = $dqb1_mg;
  $data->{DQB1_p} = $dqb1_p;
  $data->{DRB1_mg} = $drb1_mg;
  $data->{DRB1_p} = $drb1_p;
  $data->{p10} = $p10of10;
  $data->{p9} = $p9of10;
  return $data;
}
