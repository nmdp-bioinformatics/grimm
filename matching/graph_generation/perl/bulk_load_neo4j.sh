#!/usr/bin/env bash

# Set JAVA_HOME
if [ "$JAVA_HOME" ]; then
	echo "Java Home is set at " ${JAVA_HOME}
else
	echo "Exiting. JAVA_HOME env variable is not set."
	exit
fi

# Set NEO4J_HOME 
if [ "$NEO4J_HOME" ]; then
	echo "Neo4j is at " ${NEO4J_HOME}
else
	echo "Exiting. NEO4J_HOME env variable is not set."
	exit
fi

# DATA_HOME will have the csv output and graph database and config dir
DATA_HOME=${PWD}/output
mkdir -p ${DATA_HOME}

NEO4J_IMPORT_CMD=${NEO4J_HOME}/bin/neo4j-import
NEO4J_CMD=${NEO4J_HOME}/bin/neo4j
NEO4J_ADMIN_CMD=${NEO4J_HOME}/bin/neo4j-admin
NEO4J_DEFAULT_PASSWORD=ontological

NEO4J_ACTIVE_DATABASE=graphmatchwmda.db
GRAPH_PATH=${DATA_HOME}/databases/${NEO4J_ACTIVE_DATABASE}
IMPORT_PATH=${DATA_HOME}/graph

${NEO4J_CMD} stop

# Remove the previous database
rm -rf ${GRAPH_PATH}/*

# import CVS into the graph database
${NEO4J_IMPORT_CMD} --into ${GRAPH_PATH} \
        --id-type INTEGER \
        --nodes ${IMPORT_PATH}/MUUG.csv \
        --nodes ${IMPORT_PATH}/MUUG_1.csv \
        --nodes ${IMPORT_PATH}/SLG.csv \
        --nodes ${IMPORT_PATH}/ALLELE.csv \
        --nodes ${IMPORT_PATH}/Subject.csv \
        --relationships:MUUG_LIKELIHOOD ${IMPORT_PATH}/MUUG_LIKELIHOOD.csv \
        --relationships:SLG_LIKELIHOOD ${IMPORT_PATH}/SLG_LIKELIHOOD.csv \
        --relationships:M1 ${IMPORT_PATH}/M1.csv \
        --relationships:SA ${IMPORT_PATH}/SA.csv

# Create conf directory
NEO4J_CONF_DIR=${DATA_HOME}/conf
mkdir -p ${NEO4J_CONF_DIR}
# Copy the config to the conf directory
cp neo4j/conf/neo4j.conf.template ${NEO4J_CONF_DIR}/neo4j.conf
# Update neo4j.conf file
echo dbms.directories.data=${DATA_HOME} >> ${NEO4J_CONF_DIR}/neo4j.conf
echo dbms.active_database=${NEO4J_ACTIVE_DATABASE} >> ${NEO4J_CONF_DIR}/neo4j.conf

# Set default password to $NEO4J_DEFAULT_PASSWORD
rm -f output/dbms/auth
NEO4J_CONF=${NEO4J_CONF_DIR} ${NEO4J_ADMIN_CMD} set-initial-password ${NEO4J_DEFAULT_PASSWORD}

# start neo4j
NEO4J_CONF=${NEO4J_CONF_DIR} ${NEO4J_CMD} start
