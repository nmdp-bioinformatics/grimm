from py2neo import Graph
from ImputeMultiRace import Imputation
from GetBestHapMultiRace import GetBestOption
import json
import os


# Create output dir if it doesn't exist
os.makedirs('output', exist_ok=True)

json_file = 'neo4j.json'
config_data = open(json_file).read()
credentials = json.loads(config_data)

# Add here the pops you impute
pops = ['CAU']
# pops = ['AAFA','AFB','AINDI','AISC','ALANAM','AMIND','CARB','CARHIS','CARIBI','FILII','KORI','JAPI','MENAFC',
#        'NAMER','NCHI','SCAHIS','SCAMB','SCSEAI','VIET','AFA','API','CAU','HIS','NAM']

# Find options above epsilon
graph = Graph(credentials["graph"], user=credentials["user"], password=credentials["password"])

# Impute Donors
fname = '../validation/wmda/data/don.gl.txt'
test = Imputation(graph, pops)
test.impute_file(fname)

# Impute Patients
fname = '../validation/wmda/data/pat.gl.txt'
test = Imputation(graph, pops)
test.impute_file(fname)
