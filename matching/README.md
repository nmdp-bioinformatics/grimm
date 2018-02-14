
# Match Results
the WMDA results look like this
P000001;D000001;M;0;M;0;M;0;M;0;M;0;0;0
1000 patients x 10000 donors = 1M rows

The single-locus match grade values are A, M, P
The single-locus 2-allele match likelihoods range from 0-100
  A => always 100
  M => always 0
  P => ranges 0-100

## Perl code


 - generate_matchgraph.pl 
creates the graph in ./graph
 - bulk_load_neo4j.sh
loads the graph from ./graph
 - match_results.pl
computes the results in a format that can be compared to ../graph_generator/data/wmda/set3.consensus.txt


# Walk-through of Running Match

## Pre-requisites
  - Generate Imputation files for both Patient/Donors. See (Multi-race Imputation)[../../multi_race_impute]
  - Perl  
	- Install Perl Brew from http://perlbrew.pl/
	```
	curl -L https://install.perlbrew.pl | bash
	```
	
	Add perlbrew to your bashrc

	```
	source ~/perl5/perlbrew/etc/bashrc
	```

	Install the latest stable build of Perl
	
	```
	perlbrew init
	perlbrew install perl-5.26.0
	```

	Use the installed Perl
	
	```
	perlbrew list
	perlbrew use perl-5.26.0

	perl --version
	```

	Install the Perl libraries used by matching.

	```
	cpan install REST::Client
	cpan install JSON
	cpan install JSON::Parse
	cpan install Math::Round
	```
## Matching

Generate the CSV files for loading into Neo4J

```
cd graph-imputation-match/matching/graph_generation/perl
mkdir -p output/graph
time ./generate_matchgraph.pl
./generate_matchgraph.pl  141.92s user 14.39s system 98% cpu 2:39.04 total
```

Load the CSV files to Neo4J

```
time ./bulk_load_neo4j.sh

	./bulk_load_neo4j.sh  97.07s user 6.56s system 202% cpu 51.289 total
```

Visit http://localhost:7474/browser/ to verify database is working.

Run the matcher program from `graph-imputation-match/matching/search`.

```
cd ../../search
mkdir -p output
time ./match_results.pl

	./match_results.pl  7447.10s user 61.19s system 98% cpu 2:07:09.69 total

```

This produces a `mr.txt` file that holds the match results.

```
ls -lah output/mr.txt
	-rw-r--r--  1 pbashyal  827464065   383M Aug 26 13:29 output/mr.txt

```

Compare the results from the WMDA consensus results.

```
time ./compare_results.pl

	./compare_results.pl  293.36s user 351.14s system 78% cpu 13:37.46 total

```

This will produce a `output/cmp.txt` file. This will list the difference in the Neo4J version versus Consensus rsults.

```
wc -l output/cmp.txt
       0 output/cmp.txt

```
