""" docstring for module """
import argparse  # for command line arguments
from py2neo import Graph
from impute.impute import Imputation
import json


def main():
    """This is run if file is directly executed, but not if imported as
    module. Having this in a separate function  allows importing the file
    into interactive python, and still able to execute the
    function for testing"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file",
                        required=True,
                        help="Input file",
                        type=str)

    parser.add_argument("-l", "--url",
                        default="http://localhost:7474/data/",
                        help="URL for neo4j database - Default ",
                        type=str)

    parser.add_argument("-u", "--user",
                        default="neo4j",
                        help="Username for neo4j database",
                        type=str)

    parser.add_argument("-p", "--password",
                        default="new",
                        help="Password for neo4j database",
                        type=str)

    parser.add_argument("-e", "--epsilon",
                        help="Epsilion used for trimming",
                        default=0.5,
                        type=float)

    parser.add_argument("-c", "--config",
                        help="Option for providing configuration file",
                        type=str)

    parser.add_argument("-v", "--verbose",
                        help="Option for running in verbose",
                        default=False,
                        type=bool)

    args = parser.parse_args()
    inputfile = args.file
    epsilon = args.epsilon
    config = args.config
    password = args.password
    user = args.user
    verbose = args.verbose

    if(config):
        config_data = open(config).read()
        credentials = json.loads(config_data)
        password = credentials['password']
        user = credentials['user']

    # Create a new neo4j graph object with py2neo
    graph = Graph(args.url, user=user, password=password)

    # Create a new imutation object
    # This object is what made with Yoram's code
    imp = Imputation(graph, verbose)
    inputdata = (l.strip() for l in open(inputfile, "r") if l.strip())

    for line in inputdata:
        line.rstrip()
        line.strip('\n')
        [subid, glstring] = line.split("%")
        res = imp.comp_cand(glstring, epsilon)
        len_res = len(res['Haps'])
        for i in range(len_res):
            print(subid, "^".join(res['Haps'][i]), res['Probs'][i])


if __name__ == '__main__':
    """The following will be run if file is executed directly,
    but not if imported as a module"""
    main()
