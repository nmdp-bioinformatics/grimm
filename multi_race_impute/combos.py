# '''
# Created on August 7th, 2017
# ** WILL BE DELETED WITH CODE REORG **
# ** WILL BE ADDED TO impute.py      **
# @author: Mike Halagan
# '''
from itertools import permutations
import numpy as np


def cmp_to_key(mycmp):
    'Convert a cmp= function into a key= function'
    class K(object):
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0  
        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0
    return K


def getArrayValue(a):
    """
    B = 1, C = 2, A = 3, DRB1 = 4, DQB1 = 5
    """
    dic = {1: 5, 2: 4, 3: 3, 4: 2, 5: 1}
    cnt = list()
    for i in a:
        ln = 1
        for j in i:
            ln *= dic[j]
        cnt.append(ln)
    return max(cnt)


def cmp_items(a, b):

    if len(a) > len(b):
        return 1
    elif len(a) == len(b):
        a_v = getArrayValue(a)
        b_v = getArrayValue(b)
        if a_v > b_v:
            return -1
        elif a_v == b_v:
            m_a = max(a, key=lambda x: len(x))
            m_b = max(b, key=lambda x: len(x))
            if len(m_a) > len(m_b):
                return -1
            elif a_v == b_v:
                return 0
            else:
                return 1
        else:
            return 0
        return 0
    else:
        return -1


class LociCombos(object):
    '''
    Creates all of the 5-locus combinations and order the list

    * Get the list *

        from impute import LociCombos
        loci_combos = LociCombos().get_combos()
        for combo in loci_combos:
            ...

    * Print the order *

        from loci_combos import LociCombos
        combo_object = LociCombos()
        combo_object.print_order()

    '''
    def __init__(self):

        loci_a = np.array([1, 2, 3, 4, 5])

        self.full_list = list()

        # Loop through all [1 loc] 4-loc combos..
        for combo in permutations(loci_a, 1):

            # create numpy array from loci in combo
            b = np.array(combo)

            # find difference between combo and
            # starting locus set (loci_a)
            # len(diff1) == 4
            diff1 = np.setdiff1d(loci_a, b)

            # Loop through all [1 loc][2 loc][2 loc] combos
            for combo2 in permutations(diff1, 2):
                b2 = np.array(combo2)
                diff2 = np.setdiff1d(diff1, b2)
                for cmb3 in permutations(diff2, 2):
                    tmp = sorted(list([list(sorted(z)) for z in [b, b2, cmb3]]))
                    if tmp not in self.full_list:
                        self.full_list.append(tmp)

            # Loop through all [1 loc][2 loc][1 loc][1 loc] combos
            for combo2 in permutations(diff1, 2):
                b2 = np.array(combo2)
                diff2 = np.setdiff1d(diff1, b2)
                for combo3 in permutations(diff2, 1):
                    b3 = np.array(combo3)
                    diff3 = np.setdiff1d(diff2, b3)
                    tmp = sorted(list([list(sorted(z)) for z in [b, b2, b3, diff3]]))
                    if tmp not in self.full_list:
                        self.full_list.append(tmp)

        # Loop through all [3 loc] 2-loc combos...
        for combo in permutations(loci_a, 3):
            b = np.array(combo)
            diff1 = np.setdiff1d(loci_a, b)

            # Loop through all [3 loc][2 loc] combos
            for combo2 in permutations(diff1, 2):
                tmp = sorted(list([list(sorted(z)) for z in [b, combo2]]))
                if tmp not in self.full_list:
                    self.full_list.append(tmp)

            # Loop through all [3 loc][1 loc][1 loc] combos
            for combo2 in permutations(diff1, 1):
                b2 = np.array(combo2)
                diff2 = np.setdiff1d(diff1, b2)
                for combo3 in permutations(diff2, 1):
                    tmp = sorted(list([list(sorted(z)) for z in [b, b2, combo3]]))
                    if tmp not in self.full_list:
                        self.full_list.append(tmp)

        # Loop through all [4 loc][1 loc] combos
        for combo in permutations(loci_a, 4):
            b = np.array(combo)
            diff1 = np.setdiff1d(loci_a, b)
            tmp = sorted(list([list(sorted(z)) for z in [b, diff1]]))
            if tmp not in self.full_list:
                self.full_list.append(tmp)

        # Sort the combos based on...
        # 1) How many combos are there
        #   - ex. len([["A"],["C~B"],["DRB1~DQB1"]]) = 3 would be
        #  sorted ahead of len([["A"],["C"],["B"],["DRB1~DQB1"]]
        #
        # 2) What is the value of the combo based on the loci
        #  [['A', 'B', 'C', 'DRB1'], ['DQB1']] - length of 2 - value of 120
        #  [['A', 'B', 'C', 'DQB1'], ['DRB1']] - length of 2 - value of 60
        #  So [['A', 'B', 'C', 'DRB1'], ['DQB1']] will be sorted first
        #
        # 3) What's the max size of a list
        #  [['A', 'B', 'C', 'DQB1'], ['DRB1']] max length of 4
        #  [['A', 'B', 'C'], ['DQB1', 'DRB1']] max length of 3
        # Therefore [['A', 'B', 'C', 'DQB1'], ['DRB1']] will be sorted first
        #
        self.full_list = sorted(self.full_list, key=cmp_to_key(cmp_items))

    def get_combos(self):
        """ Get the combo list back
        """
        return self.full_list

    def print_order(self):
        """ Print the order of the loci combos
        """
        for l in self.full_list:
            m_a = max(l, key=lambda x: len(x))
            print(l, len(l), getArrayValue(l), len(m_a), sep="\t")

    def missing_loci(self, loci):
        """
        Missing loci should never be alone,
        they should always be with other loci
        """
        new_list = list()
        for combo in self.full_list:
            n = 0
            for loc in loci:
                x = 0
                for c in combo:
                    if loc in c and len(c) > 1:
                        x += 1

                if x != 0:
                    n += 1
            if n == len(loci):
                new_list.append(combo)
        return(new_list)

