#!/bin/sh

# set this to where the csv and graph will live
DATA_HOME=/home/polina/PycharmProjects/Dictionary

# path to use from laptop
NEO4J_HOME=/usr

#IMPORTCMD=$NEO4J_HOME/bin/neo4j-admin import
IMPORTCMD=$NEO4J_HOME/bin/neo4j-import
NEO4JCMD=$NEO4J_HOME/bin/neo4j
GRAPHPATH=/var/lib/neo4j/data2/databases/graph.db/

GRAPHOLD=/graph/databases/imputewmda.db
CSVPATH=$DATA_HOME

$NEO4JCMD stop

rm -rf $GRAPHPATH/
rm -rf $GRAPHOLD/



# import graph

set -x
$IMPORTCMD --into $GRAPHPATH \
	 --id-type INTEGER \
	--nodes $CSVPATH/Nodes_Data.csv\
	--relationships $CSVPATH/Edges_Data.csv\
	--relationships $CSVPATH/Edges_Data_P2G10.csv\


#--nodes $CSVPATH/All_Nodes.csv \
#--relationships $CSVPATH/All_Top.csv \


# start neo4j
$NEO4JCMD start
