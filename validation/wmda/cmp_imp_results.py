'''
Created 2017-03-07

@author: mmaiers-nmdp
'''
import gzip
import csv
import argparse  # for command line arguments

class Mrimp():
    def __init__(self, i):
        self.i = i
        # genotype dictionary
        self.gd = {}
    def addrow(self, h1, h2, p):
        self.h1 = h1
        self.h2 = h2
        # cannonical form of phased genotype
        g = '+'.join(sorted([h1, h2]))
        self.gd[g]=float(p)
    def gd(self):
        return self.gd

def main():
    # open oldfile (haplogic imputation) 
    # open newfile (New Imputation Method)
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--oldfile",
                        required=False,
                        default = "data/out.callHFS.ngf.fix.2015-07-22.gz",
                        help="Old file",
                        type=str)
    parser.add_argument("-n", "--newfile",
                        required=False,
                        default = "../../g2gl/pat.impute.txt",
                        help="New file",
                        type=str)
    args = parser.parse_args()
    oldfile = args.oldfile
    newfile = args.newfile
    
    nd = {} # new imputation output dict
    with open(newfile, "rt") as nf:
        reader = csv.reader(nf, delimiter=',')
        for line in reader:
            if len(line) < 3:
                break
            [i, h1, f1, p1, h2, f2, p2] = line
            f1 = float(f1)
            f2 = float(f2)
            if not i in nd:
                nd[i] = Mrimp(i)
            # HW
            if h1 != h2:
                gf = 2 * f1 * f2 
            else:
                gf = f1 * f2 
            nd[i].addrow(h1,h2,gf)


    od = {} # old imputation output dict
    with gzip.open(oldfile, "rt") as of:
        reader = csv.reader(of, delimiter='\t')
        for line in reader:
            [i, id_typ, h1, f1, p1, h2, f2, p2] = line[0:8]
            h1 = '~'.join(sorted(h1.split('~')))
            h2 = '~'.join(sorted(h2.split('~')))
            f1 = float(f1)
            f2 = float(f2)

            if h1 == h2:
                p = f1 * f2
            else:
                p = 2 * f1 * f2 
            hfp1 = [h1, f1, p1]
            hfp2 = [h2, f2, p2]
            if not i in od:
                od[i] = Mrimp(i)
            od[i].addrow(h1, h2, p)

    for i in od:
        if not i in nd:
            print("NONIM", i)
        else: 
            ngd = nd[i].gd
            ogd = od[i].gd
            for g in ogd:
                if not g in ngd:
                    print("ONLYOLD", i, ogd[g], g)
            for g in ngd:
                if not g in ogd:
                    print("ONLYNEW", i, ngd[g], g)
                else:  
                    dif = abs(ngd[g]-ogd[g])
                    ndif = dif/ogd[g]
                    if ndif > 1e-6:
                        print("BIGDIFF", i, dif, ngd[g], ogd[g], g)
                    else:
                        print("DIFF", i, dif, ngd[g], ogd[g], g)
          

if __name__ == '__main__':
    """The following will be run if file is executed directly,
    but not if imported as a module"""
    main()


