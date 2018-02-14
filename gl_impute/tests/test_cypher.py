
from hf_cypher.cypherQuery import CypherQuery
from nose.tools import assert_equal, assert_true
import inspect
import os
import json


json_file = os.path.join(os.path.dirname(__file__), 'cypher-expected.json')
file_data = open(json_file).read()
expected = json.loads(file_data)


def test_constructer():
    cypher = CypherQuery()
    assert_true(inspect.isclass(CypherQuery))
    assert_equal(cypher.__class__.__name__, "CypherQuery")
    assert_equal(str(type(cypher)), "<class \'hf_cypher.cypherQuery.CypherQuery\'>")


def test_get_locus():
    cypher = CypherQuery()
    for allele in expected["getLocus"].keys():
        assert_equal(cypher.getLocus(allele), expected["getLocus"][allele], allele)


def test_get_match():
    cypher = CypherQuery()
    for haplo in expected["getMatch"].keys():
        assert_equal(cypher.getMatch([haplo]), expected["getMatch"][haplo], haplo)


def test_get_path():
    cypher = CypherQuery()
    for haplo in expected["getPath"].keys():
        assert_equal(cypher.getPath([haplo]), expected["getPath"][haplo], haplo)


def test_build_query():
    cypher = CypherQuery()
    for haplo in expected["buildQuery"].keys():
        assert_equal(cypher.buildQuery([haplo]), expected["buildQuery"][haplo], haplo)
