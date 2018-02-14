#!/usr/bin/env perl
##############################################################################
# SCRIPT NAME:	match_results.pl
# DESCRIPTION:	run 10/10, 9/10 and 5 single-locus queries to generate WMDA
#               match results
#
# the WMDA results look like this
# P000001;D000001;M;0;M;0;M;0;M;0;M;0;0;0
# 1000 patients x 10000 donors = 1M rows
#
#
#
# DATE WRITTEN: 2017-08-10
# WRITTEN BY:   Martin Maiers
#
##############################################################################
use strict;    # always
use warnings;  # or else
use REST::Client;
use JSON;
use JSON::Parse 'parse_json';
use MIME::Base64;
use Math::Round;

my $v=1; #verbose
my %D;
my %P;
my %M;


my $outfile = "./output/mr.txt";
open OUTFILE, ">$outfile" or die "$!: $outfile";
#
# setup neo4j client
#
my $service_url = 'http://localhost:7474'; # location of the service
my $service_up = encode_base64(join (':', "neo4j", "ontological"));
my $client = REST::Client->new({ host => $service_url,});
$client->addHeader('Content-Type', 'application/json;charset=UTF-8');
$client->addHeader('Accept', 'application/json');
$client->addHeader('Authorization', "Basic $service_up");


#
# patients
#
my $cypher = qq/match (p:Subject {id_typ:'P'}) return  p.id/;
foreach my $row (doCypher($client, $cypher)) {
  my ($pid) = @{$row->{row}};
  $P{$pid}++;
}

#
# donors
#
$cypher = qq/match (d:Subject {id_typ:'D'}) return  d.id/;
foreach my $row (doCypher($client, $cypher)) {
  next unless defined $row;
  my ($did) = @{$row->{row}};
  $D{$did}++;
}

#

#
# 10/10 query
#
$cypher = qq/match (d:Subject {id_typ:'D'})-[dml:MUUG_LIKELIHOOD]- 
	(m:MUUG)- [rml:MUUG_LIKELIHOOD]-(r:Subject {id_typ:'P'})
        return  d.id, r.id, round(100*sum(dml.prob*rml.prob))/;

foreach my $row (doCypher($client, $cypher)) {
  my ($did, $pid, $p10of10) = @{$row->{row}};
  $M{$pid}{$did}{p10of10} = $p10of10;
}

#
# 9/10 query
#
$cypher = qq/match (d:Subject {id_typ:'D'})-[dml:MUUG_LIKELIHOOD]-
	(dm:MUUG)-[:M1]-(:MUUG_1)-[:M1]-(rm:MUUG)-
	[rml:MUUG_LIKELIHOOD]-(r:Subject {id_typ:'P'})
        return  d.id, r.id, round(100*sum(dml.prob*rml.prob))/;


foreach my $row (doCypher($client, $cypher)) {
  my ($did, $pid, $p9of10) = @{$row->{row}};
  $M{$pid}{$did}{p9of10} = $p9of10;
}

#
# single locus queries
#
my @loci = qw/A C B DRB1 DQB1/;

foreach my $loc (@loci) {
  $cypher = qq/match (d:Subject {id_typ:'D'})-
        [dsl:SLG_LIKELIHOOD]- (s:SLG {locus:'$loc'})- [rsl:SLG_LIKELIHOOD]-
        (r:Subject {id_typ:'P'})
        return  d.id, r.id, 100*sum(dsl.prob*rsl.prob)/;

  foreach my $row (doCypher($client, $cypher)) {
    my ($did, $pid, $p) = @{$row->{row}};
    my $mg = "M";
    $mg = "P" if $p;
    $mg = "A" if $p >=99.9999999999;  # tune this
    my $rounded = round($p);
    $M{$pid}{$did}{$loc}{p} = $rounded;
    $M{$pid}{$did}{$loc}{mg} = $mg;
  }
}


foreach my $pid (sort keys %P) {
  foreach my $did (sort keys %D) {
    my $a_mg = defined $M{$pid}{$did}{A}{mg} ? $M{$pid}{$did}{A}{mg} : 'M';
    my $c_mg = defined $M{$pid}{$did}{C}{mg} ? $M{$pid}{$did}{C}{mg} : 'M';
    my $b_mg = defined $M{$pid}{$did}{B}{mg} ? $M{$pid}{$did}{B}{mg} : 'M';
    my $drb1_mg = defined $M{$pid}{$did}{DRB1}{mg}  ? $M{$pid}{$did}{DRB1}{mg} : 'M';
    my $dqb1_mg = defined $M{$pid}{$did}{DQB1}{mg} ? $M{$pid}{$did}{DQB1}{mg} : 'M';
    my $a_p = defined $M{$pid}{$did}{A}{p} ? $M{$pid}{$did}{A}{p} : 0;
    my $c_p = defined $M{$pid}{$did}{C}{p} ? $M{$pid}{$did}{C}{p} : 0;
    my $b_p = defined $M{$pid}{$did}{B}{p} ? $M{$pid}{$did}{B}{p} : 0;
    my $drb1_p = defined $M{$pid}{$did}{DRB1}{p}  ? $M{$pid}{$did}{DRB1}{p} : 0;
    my $dqb1_p = defined $M{$pid}{$did}{DQB1}{p} ? $M{$pid}{$did}{DQB1}{p} : 0;
    my $p9of10 = defined $M{$pid}{$did}{p9of10} ? $M{$pid}{$did}{p9of10} : 0;
    my $p10of10 = defined $M{$pid}{$did}{p10of10} ? $M{$pid}{$did}{p10of10} : 0;
    print OUTFILE join (';', 
 	$pid, $did, 
        $a_mg, $a_p, 
        $c_mg, $c_p, 
        $b_mg, $b_p, 
        $drb1_mg, $drb1_p, 
        $dqb1_mg, $dqb1_p, 
        $p9of10, $p10of10,), "\n";
  }
}


exit 0;


sub doCypher {
  my ($client, $query) = @_;
  print STDERR "cypher: $cypher\n" if $v;
  my $request = { statements => [ { statement => $cypher } ] };
  my $json_request = JSON::to_json($request);
  $client->POST('/db/data/transaction/commit', $json_request, {});
  my $json_response = $client->responseContent;
  #my $response = JSON::from_json($json_response);
  my $response = parse_json($json_response);
  # reference to response
  return undef unless defined $$response{results}[0]{data};
  return @{$$response{results}[0]{data}};
}
