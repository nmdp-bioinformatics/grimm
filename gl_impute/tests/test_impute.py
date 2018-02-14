

from py2neo import Graph
from impute.impute import Imputation
from nose.tools import assert_equal, assert_true
import inspect
import os
import json

"""
ex. of py2neo.json - used for loading graph credentials

{
    "graph":"http://localhost:7474/data/",
    "password":"new",
    "user":"neo4j"
}

"""
json_file = os.path.join(os.path.dirname(__file__), 'py2neo.json')
config_data = open(json_file).read()
credentials = json.loads(config_data)


def test_constructer():
    graph = Graph(credentials["graph"], user=credentials["user"], password=credentials["password"])
    impute = Imputation(graph, 1)
    assert_true(inspect.isclass(Imputation))
    assert_equal(impute.__class__.__name__, "Imputation")
    assert_equal(str(type(impute)), "<class \'impute.impute.Imputation\'>")


def test_GL2Haps():
    graph = Graph(credentials["graph"], user=credentials["user"], password=credentials["password"])
    impute = Imputation(graph, 1)
    test_gl = "A*01:01/A*01:02+A*02:01/A*02:02/A*02:05^B*08:01/B*08:04+B*07:02/B*07:01^C*07:01+C*07:02"
    expected = {'Genotype': [['A*01:01/A*01:02', 'B*08:01/B*08:04', 'C*07:01'], ['A*02:01/A*02:02/A*02:05', 'B*07:02/B*07:01', 'C*07:02']], 'N_Loc': 3}
    haplos = impute.gl2haps(test_gl)
    assert_equal(haplos['N_Loc'], expected['N_Loc'])
    assert_equal(haplos['Genotype'], expected['Genotype'])
    assert_equal(haplos['Genotype'][0], expected['Genotype'][0])
    assert_equal(haplos['Genotype'][1], expected['Genotype'][1])


def test_impute():
    graph = Graph(credentials["graph"], user=credentials["user"], password=credentials["password"])
    impute = Imputation(graph, 1)
    test_gl = "A*01:01/A*01:02+A*02:01/A*02:02/A*02:05^B*08:01/B*08:04+B*07:02/B*07:01^C*07:01+C*07:02"
    imputed_genotypes = impute.comp_cand(test_gl, .50)
    assert_true(len(imputed_genotypes['Haps']) > 2)



