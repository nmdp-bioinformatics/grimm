
# coding: utf-8

# #### TODO
# * ~~Remove little g from allele names~~
# * ~~Add toplinks~~
# * ~~Set of ABCRQ in FULL_LOCI~~
# * ~~Remove ChainMap~~
# * ~~Use namedtuples~~
# * ~~Ignore 0 HF frequency from Input file~~
# * Compare with Matlab output 
# * Document the code
# * Extract a module
#     * Accepts the list of FULL_LOCI and initial Freq. file

from itertools import combinations, count 
from collections import defaultdict, namedtuple
from operator import indexOf, add, truediv
import csv
import os
import pandas as pd

os.makedirs('output/csv', exist_ok=True)

HaplotypeNode = namedtuple('HaplotypeNode', 'node_id, freq_list, allele_list, parents, top_links')
TopLink = namedtuple('TopLink', 'node_id, haplotype')
ParentLoci = namedtuple('ParentLoci', 'loci, haplotype, p')


# #### nCr: n Choose r

def nCr(nchars, r):
    return [''.join(x) for x in combinations(nchars, r)]


FULL_LOCI='ABCQR'


all_combo_list = [FULL_LOCI]
for i in range(len(FULL_LOCI) -1, 0, -1):
    all_combo_list.extend(nCr(FULL_LOCI, i))
len(all_combo_list)


sequence = count(0)


loci_combo_map =  {combo_list: defaultdict(lambda : HaplotypeNode(node_id = next(sequence),
                                                                  freq_list = [],
                                                                  allele_list=[],
                                                                  parents=[],
                                                                  top_links=set())) 
                   for combo_list in all_combo_list }


# #### Verify
# There should be 31 loci combos. The keys should be all the different combinations.

print(len(loci_combo_map))
loci_combo_map.keys()


def make_allele_list(haplotype):
    hl = haplotype.split('~')
    # Remove g from the allele names
    hl = list(map(lambda h: h[:-1] if(h[-1] == 'g') else h, hl))
    hl.sort()
    return hl


# #### Verify
# Test that the make_allele_list works


# #### Load Frequency file
# Load the frequency file into the FULL_LOCI dictionary with the cleaned up haplotype as the key.

haplist_overall = {} # list of haplotypes across all populations
pop_hap_combos = {}
pops = {}
# pops = ['AAFA','AFB','AINDI','AISC','ALANAM','AMIND','CARB','CARHIS','CARIBI','FILII','KORI','JAPI','MENAFC',
#         'NAMER','NCHI','SCAHIS','SCAMB','SCSEAI','VIET','AFA','API','CAU','HIS','NAM']
# pops = ['VIET', 'KORI', 'JAPI', 'FILII']
# pops = ['VIET', 'FILII']
freqfile = 'data/wmda/hpf.csv'
with open(freqfile) as f:
    for hap_line in f:
        haplotype, pop, freq = hap_line.split(',')
        if haplotype == "hap":
            continue
        freq = float(freq)
        # Ignore lines with 0 freq
        if freq == 0.0:
            continue
        #if freq < 0.0001: # REMOVE THESE TWO LINES TO GENERATE COMPLETE EDGE
        #    continue
        hap_list = make_allele_list(haplotype)
        haplotype = '~'.join(hap_list)
        pop_haplotype = pop + '-' + haplotype
        pops[pop] = 1
        haplist_overall[haplotype] = 1
        pop_hap_combos[pop_haplotype] = freq

for haplotype in haplist_overall:
    # make frequency array
    freqs = []
    for pop in pops:
        # check to see if population + haplotype combo exists
        pop_haplotype = pop + '-' + haplotype
        if pop_haplotype in pop_hap_combos:
            freqs.append(pop_hap_combos[pop_haplotype])
        else:
            freqs.append(0)
    hap_list = make_allele_list(haplotype)        
    loci_combo_map[FULL_LOCI][haplotype] = HaplotypeNode(node_id = next(sequence),
                                                             freq_list = freqs,
                                                             allele_list = hap_list,
                                                             parents = None,
                                                             top_links= set())


# #### Verify
# Make sure one of the items in the map looks ok.

h = next(iter(loci_combo_map[FULL_LOCI].items()))
print(h)


# Find the index of a loci name in `FULL_LOCI`. ie. `BC` in `ABCRQ` is `[1,2]`

def find_indices(loci):
    return list(map(lambda x: indexOf(FULL_LOCI, x), loci))


# #### Verify
# Print out the indices to see that they make sense.

for loci in all_combo_list:
    print(find_indices(loci))


find_indices('ABCQ')


# Given a list of indices find the corresponding loci combo

def find_loci_name(loci_indices):
    return ''.join([FULL_LOCI[i] for i in sorted(loci_indices)])


find_loci_name([1,2])


# Find all the parents of the loci with loci_indices against loci_allele_list. 

def find_parents(loci_indices, loci_allele_list, p):
    fi = find_indices(FULL_LOCI)
    new_loci = set(fi).difference(loci_indices)
    parents = []
    for i in list(new_loci):
        new_loci = loci_indices.copy()
        new_loci.append(i)
        #print(new_loci)
        loci_name = find_loci_name(new_loci)
        #print(loci_name)
        haplotype = '~'.join([loci_allele_list[i] for i in sorted(new_loci)])
        parents.append(ParentLoci(loci=loci_name, haplotype=haplotype, p=p))
    return parents


# #### Verify
# Find parents of a random haplotype from C

# haplotype = 'A*01:01~B*15:02~C*01:02~DQB1*03:01~DRB1*12:01'
haplotype = 'A*24:02~B*52:01~C*12:02~DQB1*06:01~DRB1*15:02'
# haplotype = 'A*01:01~B*15:01~C*01:02~DQB1*05:01~DRB1*01:01'


list(loci_combo_map[FULL_LOCI].items())[99]


# For each loci block, get all the allele combinations. Create parents link, and accumulate the freqs to create new Haplotype Nodes.

for loci in all_combo_list:
    if loci != FULL_LOCI: # FULL_LOCI is already in the dictionary
        loci_map = loci_combo_map[loci]
        loci_indices = find_indices(loci)
        for full_haplotype_name, full_haplotype_node in loci_combo_map[FULL_LOCI].items():
            allele_list = full_haplotype_node.allele_list
            haplotype = '~'.join([allele_list[i] for i in loci_indices])
            haplotype_node = loci_map[haplotype]            
            
            parents = find_parents(loci_indices, allele_list, full_haplotype_node.freq_list)
            parents_list = haplotype_node.parents
            for parent in parents:
                parents_list.append(parent)
                
            top_links = haplotype_node.top_links
            top_links.add(TopLink(node_id=full_haplotype_node.node_id, haplotype=full_haplotype_name))

            # have to make freq_list with all zeros otherwise map add function returns empty list
            new_freq_list = []
            if not haplotype_node.freq_list:
                for pop in pops:
                    new_freq_list.append(0)
            else:
                new_freq_list = haplotype_node.freq_list
            freq_sum = list(map(add,new_freq_list,full_haplotype_node.freq_list))

            new_node = HaplotypeNode(node_id=haplotype_node.node_id,
                                    freq_list=freq_sum,
                                    allele_list=None,
                                    parents=parents_list,
                                    top_links=top_links)

            loci_map[haplotype] = new_node


haplotype="A*24:02~B*52:01~C*12:02~DQB1*06:01"
loci_combo_map['ABCQ'][haplotype]
# haplotype="A*24:02~B*52:01~C*12:02~DQB1*06:01~DRB1*15:02"
# loci_combo_map[FULL_LOCI][haplotype]
# haplotype="A*24:02"
# loci_combo_map['A'][haplotype]


# #### Verify
# Print out the first ten haplotypes to make sure they look ok.

for loci in all_combo_list:
    print(loci)
    loci_map = loci_combo_map[loci]
    for k,v in list(loci_map.items())[1:10]:
        print(k, v.node_id, v.freq_list, sep="\t")


# #### Verify
# Make sure all the combos add up to 1 (or close to)

# for loci in all_combo_list:
#    print(loci, sum(map(lambda v: v[1], loci_combo_map[loci].values())), sep="\t")


# #### Build Nodes file

header = ['haplotypeId:ID(HAPLOTYPE)', 'name', 'loci:LABEL', 'frequency:DOUBLE[]']
node_file = 'output/csv/nodes.csv'
with open(node_file, mode='w') as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(header)
    for loci in all_combo_list:
        loci_map = loci_combo_map[loci]
        for haplotype, haplotype_node in loci_map.items():
            freq_array = ';'.join(map(str,haplotype_node.freq_list))
            csv_writer.writerow([haplotype_node.node_id, haplotype, loci, freq_array])


# #### Build Edges File

def dividebyzero(a, b):
    if b == 0:
        return 0
    else: 
        return a/b

# header = [':START_ID(HAPLOTYPE)', ':END_ID(HAPLOTYPE)', 'CP:DOUBLE[]', ':TYPE']
# 
# with open(edge_file, mode='w') as csvfile:
#     csv_writer = csv.writer(csvfile)
#     csv_writer.writerow(header)
edge_file = 'output/csv/edges.csv'
fq_l = list()
for loci in all_combo_list:
    if loci != FULL_LOCI:
        # FULL_LOCI is already in the dictionary
        loci_map = loci_combo_map[loci]
        for haplotype, haplotype_node in list(loci_map.items()):
            for parent in haplotype_node.parents:
                # assert(haplotype_node.freq > 0.0)
                # avoid division by zero
                cp = list(map(dividebyzero, parent.p, haplotype_node.freq_list))
                prob_array = sum(cp)
                loci_combo = parent.loci
                hap = parent.haplotype
                parent_id = loci_combo_map[loci_combo][hap].node_id
                l = list([haplotype_node.node_id, parent_id, prob_array, 'CP'])
                fq_l.extend([l])

freq_df = pd.DataFrame(fq_l, columns=["HapNode", "ParentID", "P", "CP"])
ndf = freq_df.join(freq_df.groupby(['HapNode', 'ParentID'])['P'].sum(), on=['HapNode', 'ParentID'], rsuffix='_r')
df2 = ndf[['HapNode', 'ParentID', 'P_r', 'CP']]
df2 = df2.drop_duplicates()
df2.columns = [':START_ID(HAPLOTYPE)', ':END_ID(HAPLOTYPE)', 'CP:DOUBLE[]', ':TYPE']
df2.to_csv(edge_file, index=False, header=True, line_terminator='\n')

# # #### Generate Top Links file
header = [':START_ID(HAPLOTYPE)', ':END_ID(HAPLOTYPE)', ':TYPE']
top_links_file = 'output/csv/top_links.csv'
with open(top_links_file, mode='w') as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(header)
    for loci in all_combo_list:
        if loci != FULL_LOCI: # FULL_LOCI is already in the dictionary
            # 
            loci_map = loci_combo_map[loci]
            for haplotype, haplotype_node in list(loci_map.items()):
                top_links = haplotype_node.top_links
                for top_link in top_links:
                    csv_writer.writerow([haplotype_node.node_id, top_link.node_id, 'TOP'])

