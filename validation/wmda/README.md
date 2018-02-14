

g2gl - genotype to genotype-list 

This module starts with file of haplotype frequencies to determine a 
starting set of alleles.

Then for an input genotype 
D002731%A*02:XX+A*02:XX^C*UUUU+C*UUUU^B*40:XX+B*40:XX^DRB1*04:01+DRB1*13:KAKU^DQB1*UUUU+DQB1*UUUU

1) expand to GL string with allele lists
2) remove loci with UUUU
3) reduce down to starting set of alleles


The MAC.pm library will only work if you have set this in the environment:
export PERL_LWP_SSL_VERIFY_HOSTNAME=0
