'''
Created on Feb 8, 2017

@author: mhalagan
'''
from itertools import combinations


class CypherQuery(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        haplo = "A*02:01~C*03:03~B*15:01~DRB1*11:01~DQB1*03:01"
        self.loc_map = {'A': 'A', 'B': 'B', 'C': 'C', 'DRB1': 'R', 'DQB1': 'Q'}
        self.loci_a = sorted([self.getLocus(x) for x in haplo.split("~")])
        self.return_freq = "return abcqr.name,abcqr.frequency"

        self.mapping = {}
        for L in range(1, len(self.loci_a)+1):
            for subset in combinations(self.loci_a, L):
                loci = "~".join(subset)
                if loci not in self.mapping:
                    self.mapping[loci] = list()

                loci_lc = "".join([allele.lower() for allele in subset])
                nodeName = "("+loci_lc+":`"+loci+"`)"
                self.mapping[loci].append(nodeName)

    def getLocus(self, allele):
        """ Get the locus back from any given allele
        """
        loc_allele = allele.split("*")
        return self.loc_map[loc_allele[0]]

    def getMatch(self, typing):
        loci = sorted([self.getLocus(a) for a in typing[0].split("~")])
        loci_lc = "".join([allele.lower() for allele in loci])
        match = "MATCH ("+loci_lc+")"
        return match

    def makeWhereClause(self, haplo):
        # ** This is only used by plan B ** #
        alleles = haplo.split("~")
        where_clause = "WHERE"
        numcp = 5 - len(alleles)
        for i in range(0, numcp):
            c_ind = i + 1
            where_clause = where_clause + " c" + str(c_ind) + ".CP<20000.0 AND"
        where_clause = where_clause + " abcqr.frequency>0.0"
        return where_clause

    def haplosIn(self, haplos):
        # Takes in ambigous list of haplotypes and creates WHERE clause#
        loci = sorted([self.getLocus(a) for a in haplos[0].split("~")])
        haplo_quoted = ",".join(["\"" + haplo + "\"" for haplo in haplos])
        loci_lc = "".join([allele.lower() for allele in loci])
        where_clause = "WHERE " + loci_lc + ".name IN [ " + haplo_quoted + " ]"
        return where_clause

    def getPath(self, typing):

        loci_t = sorted([self.getLocus(x) for x in typing[0].split("~")])
        missing = list()
        for locus in self.loci_a:
            if locus not in loci_t:
                missing.append(locus)

        nodes = list()
        for missing_locus in missing:
            new_loci = loci_t
            new_loci.append(missing_locus)
            loci_lc = "".join(sorted([allele.lower() for allele in new_loci]))
            loci_up = "".join(sorted(new_loci))
            nodeName = "("+loci_lc+":`"+loci_up+"`)"
            nodes.append(nodeName)

        path = list()
        path.append("[c"+str(1)+"]")
        path.append(nodes[0])

        for i in range(1, len(nodes)):
            path.append("[c"+str(i+1)+"]")
            path.append(nodes[i])

        query_path = "-".join(path)
        return(query_path)

    def buildQuery(self, typing):
        alleles = typing[0].split("~")
        # If there are missing loci...
        if len(alleles) < len(self.loci_a):
            match = self.getMatch(typing)
            path = self.getPath(typing)
            hapsin = self.haplosIn(typing)
            return(match + "-" + path + " " + hapsin + " " + self.return_freq)
        else:
            match = self.getMatch(typing)
            hapsin = self.haplosIn(typing)
            return(match + " " + hapsin + " " + self.return_freq)
