# Graph Data Generation and Load into Neo4J
## Pre-requisites

- JDK 8
	- Install JDK 1.8 from Oracle
	- Set JAVA_HOME
	     On MacOS, you can do:
		```
		export JAVA_HOME="$(/usr/libexec/java_home -v 1.8)"
		```

- Python 3
	- On MacOS install with 
		```
		brew install python3
		```

- Install Neo4J
	- On MacOS install with 
		```
		brew install neo4j
		```

	- Setup NEO4J_HOME
            Point NEO4J_HOME to the root of the NEO4J directory.
		```
		export NEO4J_HOME=/usr/local/Cellar/neo4j/3.2.2/libexec
		```

# Generate Data

- Clone the repository 
	```
	git clone https://github.com/nmdp-bioinformatics/grimm`
	cd grimm
	```

- Setup Python3 virtual environment
  Make sure `virtualenv` is installed.
    ```
    pip3 install virtualenv
    ```

    Create Virtual Environment
 
    ```
    virtualenv -p python3 venv
    source venv/bin/activate
    ```
 
    Install pandas library
    ```
    pip3 install pandas
    ```


- Download and prepare wmda data. Python script downloads reference wmda data and untars it in wmda directory
	```
	cd graph_generator/data
	python wmda_download.py
	```

- Generate nodes/edges/toplinks from the reference wmda data. The freqs file is converted to HPF format first.
	```
	cd ..
	python wmda_to_hpf_csv.py
	python generate_neo4j_wmda_hpf.py
	```

- Generated nodes and edges files are in the output/csv directory.

	```
	output
	└── csv
	    ├── edges.csv
	    ├── nodes.csv
	    └── top_links.csv
	```

- Load the nodes/graph into Neo4J database
	```
	./bulk_load_neo4j.sh
	```

- Login to Neo4j http://localhost:7474/
	- Default *username/password* is *neo4j/ontological*
