#!/usr/bin/env perl
##############################################################################
# PACKAGE NAME:	MAC.pm
# DESCRIPTION:	validation against mac service
#
# DATE WRITTEN: 2016-02-19
# WRITTEN BY:   Martin Maiers
#
#       COPYRIGHT (C) 2016 NATIONAL MARROW DONOR PROGRAM.  
#               ALL RIGHTS RESERVED        
##############################################################################
package MAC;
use strict;
use warnings;
use LWP::UserAgent;

my $imgtHlaRelease = "3.27.0";
my $expand = "false";


sub set_expand {
  my $e = shift;
  $expand = $e;
}

sub set_imgtHlaRelease {
  my $r = shift;
  $imgtHlaRelease = $r;
}

sub decode {
    my $typing = shift;
    my $rlist  = shift;

    my $ua = new LWP::UserAgent;
    $ua->agent("MAC_Client/0.1");
    my $base_url = "https://hml.nmdp.org/mac/api";
    my $url = "$base_url/decode?imgtHlaRelease=$imgtHlaRelease&expand=$expand&typing=$typing";
  my $response = $ua->request(new HTTP::Request("GET", $url));
  my $rcode = $response->code;
  my $content = $response->content;
  if ($rcode == 200) {  # OK
    my @allele_list= split ("/", $content);
    push @{$rlist}, @allele_list;
  } elsif ($rcode == 400) { # Bad Request
    push @{$rlist}, $content;
  } else {
    die "System error: code=$rcode $content\n";
  }
  return $rcode;
}

1;
