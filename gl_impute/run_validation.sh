#!/bin/sh

pip uninstall .
python setup.py clean --all
python setup.py install

# Make sure to have a py2neo.json file in tests/ that contains your db credentials 
# ex.
# {
#     "graph":"http://localhost:7474/data/",
#     "password":"new",
#     "user":"neo4j"
# } 
nosetests -v tests

time python bin/impute_file.py -c tests/py2neo.json -f ../g2gl/pat.gl.txt > pat.impute.txt 2> pat.impute.stderr
time python bin/impute_file.py -c tests/py2neo.json -f ../g2gl/don.gl.txt > don.impute.txt 2> don.impute.stderr

cut -f1 -d ' ' pat.impute.txt | sort -u | perl -ne 'END{ my $f = "../g2gl/pat.gl.txt";open(my $fh,"<",$f);while(<$fh>){chomp;my($id,$gl) = split(/\%/,$_);print $_,"\n" if !defined $h{$id};}close $fh;}chomp;$h{$_}++;' > pat.missing.txt
cut -f1 -d ' ' don.impute.txt | sort -u | perl -ne 'END{ my $f = "../g2gl/don.gl.txt";open(my $fh,"<",$f);while(<$fh>){chomp;my($id,$gl) = split(/\%/,$_);print $_,"\n" if !defined $h{$id};}close $fh;}chomp;$h{$_}++;' > don.missing.txt

wc -l *.missing.txt
