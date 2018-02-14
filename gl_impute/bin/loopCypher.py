""" docstring for module """
import glob


import argparse  # for command line arguments

from py2neo import Graph
from impute.impute import Imputation

import re

def main():
    """This is run if file is directly executed, but not if imported as 
    module. Having this in a separate function  allows importing the file 
    into interactive python, and still able to execute the 
    function for testing"""
    
    # Create a new neo4j graph object with py2neo
    graph   = Graph("http://localhost:7474/data/",user="neo4j",password="ontological")
    
    # Create a new imutation object 
    # This object is what made with Yoram's code
    imp     = Imputation(graph)
    
    epsilon = 0.0000000000000000001
    donfile = "g2gl/don.gl.txt"
    patfile = "g2gl/pat.gl.txt"
    
    inputfile = open(donfile)
    lines= inputfile.readlines()
    
    
    # compile the regex
    p = re.compile('%')
    for i in lines:
        i =  i[:-1] if i[-1]=='\n' else i
        [id,gl] = p.split(i)
        print ("id=",id,"<")
        print ("gl=",gl,"<")
    # pass in an glstring and epsilion and get back the imputed genotype frequencies
    res = imp.impute(gl, epsilon)
    N_Res=len(res['Haps'])
    
    for i in range(N_Res):
        print(id,i,''.join(res['Haps'][i]),res['Probs'][i])



if __name__ == '__main__':
    """The following will be run if file is executed directly, 
    but not if imported as a module"""
    main()