'''

@author: yoram
'''

import pandas as pa
import os
import gzip
from hf_cypher.cypherQuery import CypherQuery


class Imputation(object):
    '''
    classdocs
    '''

    def __init__(self, graph, populations):
        """Constructor
        Intialize an instance of `Imputation` with a py2neo graph
        and a `CypherQuery` object.
        """
        # graph is a link to the graph and populations is the list of populations
        self.graph = graph
        self.cypher = CypherQuery()
        self.populations = populations

    def power_find(self, n):
        """produces all powers of 2
        """
        result = []
        binary = bin(n)[:1:-1]
        for x in range(len(binary)):
            if int(binary[x]):
                result.append(x)
        return result

    def gl2haps(self, GL_String):
        # Receives a GL string adn produces a genotype in list structure
        split_hap = GL_String.split('^')
        N_Loci = len(split_hap)
        t1 = []
        t2 = []
        for i in range(N_Loci):
            curr_locus = split_hap[i].split('+')
            if len(curr_locus) == 1:
                return []
            t1.append(curr_locus[0])
            t2.append(curr_locus[1])

        Gen = [t1, t2]
        return {'Genotype': Gen, 'N_Loc': N_Loci}

    def gen_phases(self, gen, n_loci):
        # Generates all phases, but does not handle locus ambiguities
        Phases = []
        N_Phases = 2 ** (n_loci - 1)  # Total Number of phases
        exists = {}
        for i in range(0, N_Phases):
            H1 = []  # Hap lists
            H2 = []
            M1 = self.power_find(i)  # find all the powers of 2 in i
            L = [0] * n_loci  # Iitiated at 0 for all loci
            for m in M1:
                L[m] = 1
                # take a phase and set it to 1,
                # all others go to the other phase.
            for k in range(n_loci):
                H1.append(gen[L[k]][k])
                H2.append(gen[1 - L[k]][k])
            geno = "^".join(["~".join(sorted(H1)), "~".join(sorted(H2))])
            if geno not in exists:
                exists[geno] = 1
                Phases.append([sorted(H1), sorted(H2)])

        return {'Phases': Phases, 'N_Phases': N_Phases}

    def open_ambiguities(self, hap, loc):
        # This opens all allele ambiguities
        hap_new = [] # produces an empty list of haplotypes
        for k in range(len(hap)): #slit a given locus in all haps.
            split_loc = hap[k][loc].split('/')
            hap1 = hap[k]
            if len(split_loc) > 1:
                for i in range(len(split_loc)):
                    hap1[loc] = split_loc[i]
                    hap_new.append(hap1[:])
            else:
                hap_new.append(hap1[:])
        return hap_new

    def comp_hap_prob(self, Hap, N_Loc, epsilon):
        haplo_probs = self.get_haplo_freqs(Hap, epsilon)
        probs = list(haplo_probs.values())
        #probs=haplo_probs.values()
        haplos = list(haplo_probs.keys())
        if not haplo_probs:
            return {'Haps': '', 'Probs': ''}
        return {'Haps': haplos, 'Probs': probs}

    def get_haplo_freqs(self, haplos, epsilon):
        haplo_probs = {}
        all_hap = []
        for hap_cand in haplos:
            haplos_joined = ["~".join(sorted(hap)) for hap in hap_cand]
            all_hap.append(haplos_joined)
        haplo_query1 = self.cypher.buildQuery(haplos_joined)
        fq = pa.DataFrame(self.graph.data(haplo_query1))
        if not fq.empty:
            freq1_dic = fq.set_index('abcqr.name')['abcqr.frequency'].to_dict()
            haplo_probs.update(freq1_dic)

        return haplo_probs

    def get_haplo_freqs_miss(self, haplos, epsilon):
        haplo_probs = {}
        all_hap = []
        for hap_cand in haplos:
            haplos_joined = ["~".join(sorted(hap)) for hap in hap_cand]
            all_hap.append(haplos_joined)
        haplo_query1 = self.cypher.buildQuery(haplos_joined)
        haplo_query1
        fq = pa.DataFrame(self.graph.data(haplo_query1))
        if not fq.empty:
            freq1_dic = fq.set_index('abc.name')['abc.frequency'].to_dict()
            haplo_probs.update(freq1_dic)
        return haplo_probs

    def cal_prob(self, probs1, probs2, epsilon):
        # This is the part where we loop over all race combinations.
        # N*N loop, where N is the number of races/
        places = []
        for i in range(len(probs1)):
            for j in range(len(probs2)):
                if probs1[i]*probs2[j] >= epsilon and probs1[i]*probs2[j] != 0:
                    places.append([i,j])
        # places are the indices of the positions in the population array.
        return places

    def comp_phase_prob(self, phases, N_Loc, epsilon):
        # receives a list of phases and computes haps and
        # probabilties and accumulate cartesian product
        geno_seen = set([])
        hap_total = []
        p_total = []
        pop_res = []
        Prob2 = []
        for i in range(len(phases)):
            P1 = self.comp_hap_prob(phases[i][0], N_Loc, epsilon)
            # This will open locus ambiguities and comp probabilities for Hap1
            Haps1 = P1['Haps']
            Prob1 = P1['Probs']
            if len(Prob1)>0:
                P2 = self.comp_hap_prob(phases[i][1], N_Loc, epsilon)
            # This will do the same for Hap 2;
                Haps2 = P2['Haps']
                Prob2 = P2['Probs']
            for h in range(len(Prob1)):
                for k in range(len(Prob2)):
                        places = self.cal_prob(Prob1[h], Prob2[k], epsilon)
                        for i in range(len(places)):
                            p=(places[i])
                            # avoid reporting the same haplotype pair more than once
                            h1_id = (Haps1[h], self.populations[p[0]])
                            h2_id = (Haps2[k], self.populations[p[1]])
                            geno_id = tuple(sorted([h1_id, h2_id]))
                            if geno_id not in geno_seen:
                                geno_seen.add(geno_id)
                                hap_total.append([geno_id[0][0], geno_id[1][0]])
                                pop_res.append([geno_id[0][1], geno_id[1][1]])
                                # record prob in same order as associated hap & pop
                                # for example, in the WMDA imputation output, this:
                                # D000001,A*02:01~B*44:03~C*02:02~DQB1*02:01~DRB1*07:01,0.00009000,CAU,A*11:01~B*13:02~C*06:02~DQB1*05:01~DRB1*01:01,0.00006000,CAU
                                # should be this:
                                # D000001,A*02:01~B*44:03~C*02:02~DQB1*02:01~DRB1*07:01,0.00006000,CAU,A*11:01~B*13:02~C*06:02~DQB1*05:01~DRB1*01:01,0.00009000,CAU
                                if geno_id[0] == h1_id:
                                    p_total.append([Prob1[h][p[0]], Prob2[k][p[1]]])
                                else:
                                    p_total.append([Prob2[k][p[1]], Prob1[h][p[0]]])
        # p_total returns an array of N*N (N is number of populations), hap_total - pairs of haplotypes.
        # pop_res are the names of the populations
        return {'Haps': hap_total, 'Probs': p_total,'Pops': pop_res}

    def open_phases(self, haps, N_Loc):
        phases=[]
        for j in range(len(haps)):
            H1=[]
            H2=[]
            for k in range(2):
                hap_list = []
                hap_list.append(haps[j][k])
                for i in range(N_Loc):
                    hap_list = self.open_ambiguities(hap_list, i)
                if (k == 0):
                    H1.append(hap_list)
                else:
                    H2.append(hap_list)
            phases.append([H1, H2])
        return phases

    def comp_cand(self, gl_string,epsilon=0.0001):
        # receives a list of phases and computes haps and
        # probabilties and accumulate cartesian productEpsilon=0.0001
        chr = self.gl2haps(gl_string)
        chr1 = self.gen_phases(chr['Genotype'], chr['N_Loc'])
        if chr1 == []:
            return
        phases=self.open_phases(chr1['Phases'],chr['N_Loc'])

        n_res = 0
        min_res = 10
        min_epsilon = 1.e-3
        res = {'Haps': 'NaN', 'Probs': 0}
        while (epsilon > 0) & (n_res < min_res):
            epsilon /= 10
            if (epsilon < min_epsilon):
                epsilon = 0.0
            res = self.comp_phase_prob(phases,  chr['N_Loc'], epsilon)
            n_res = len(res['Haps'])

        return res

    def impute_file(self, fname):
        print("Starting Imputation!\n")
        # TODO: do the right thing if its a gzip
        f = open(fname, 'r')
        filename = os.path.basename(fname)

        fout_name= 'output/' + filename +'_out'
        print("Output Imputation file: " + fout_name)
        fout=open(fout_name,'w')

        fout1_name= 'output/' + filename + '_val'
        fout1=open(fout1_name,'w')

        fout2_name= 'output/' + filename + '_miss'
        fout2=open(fout2_name,'w')


        x = f.readlines()
        for i in range(len(x)):
            x[i] = x[i].strip('\n')
            name_gl = x[i].split('%')
            if (len(name_gl) == 2):
                res = self.comp_cand(name_gl[1], 0.0001)
                haps =  res['Haps']
                probs = res['Probs']
                pops = res['Pops']
                print(i, "Subject:", name_gl[0], len(haps))
                if (len(haps) == 0):
                    fout2.write(x[i] + '\n')
                fout1.write(str(len(haps)) + ',' + str(name_gl[0]) + '\n')
                for j in range(len(haps)):
                    # Write the next format:ID, Haplotype1, Probability 1, Race 1, Haplotype2, Probability 2, Race 2
                    # No header
                    fout.write(str(name_gl[0]) + ',' + str(haps[j][0]) + ',' + str(probs[j][0])+','+str(pops[j][0]) +
                               ',' + str(haps[j][1]) + ',' + str(probs[j][1]) + ',' + str(pops[j][1]) + '\n')
        f.close()
        fout.close()
        fout1.close()
        fout2.close()
