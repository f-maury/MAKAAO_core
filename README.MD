[![66.67 % FAIR](https://img.shields.io/badge/FAIR_assessment-66.67_%25-orange)](https://fair-checker.france-bioinformatique.fr/assessment/6841584a442568dc9750d7a7) ![Version](https://img.shields.io/badge/version-1.0.0-blue) 

# ⚠️ DISCLAIMER:  this repository is still a work in progress⚠️
* You can currently find in this repository: the MAKAAO Core dataset, the MAKAAO Knowledge Graph, the MAKAAO Lite Knowledge Graph, and the scripts used to generate them from the dataset.
* For licensing reasons, we also cannot include in this repository some of the data file used as inputs by our scripts. If you try to run the scripts out of the box, they will fail because they will not find their required files at the specified paths. In order to run the scripts, you need to obtain the missing files and manually place them at the specified paths. The missing files are:
     - `MRCONSO.RRF` produced by the <a href="https://www.nlm.nih.gov/research/umls/index.html">UMLS</a>
     - `LOINCFR_JeuDeValeurs.xlsx` produced by the French <a href="https://esante.gouv.fr/">ANS</a>, using the <a href="https://loinc.org/">LOINC</a> international release
     - `LOINCFR_JeuDeValeurs_SKOS-RDF.rdf` produced by the French <a href="https://esante.gouv.fr/">ANS</a>, using the <a href="https://loinc.org/">LOINC</a> international release
     - `en_product4.xml` produced by <a href="https://www.orpha.net/">Orphanet</a>
     - <a href="https://www.snomed.org/">SNOMED CT International</a> release OWL file
       
# MAKAAO Knowledge Graph
MAKAAO KG is a knowledge graph about autoantibodies (AAb) and autoimmune diseases (AID). It is part of the <a href=https://makaao.inria.fr>MAKAAO project</a> (MApping Knowledge about Autoimmune diseases, Autoantibodies, and Other clinical signs). This knowledge graph represents human autoimmune antibodies, their molecular targets, the diseases they are associated with, and some clinical signs related to these diseases. It is located in the `graphs/makaao_kg.rdf`.

# MAKAAO Lite Knowledge Graph
MAKAAO Lite KG is a simplified, lighter version of MAKAAO KG. MAKAAO Lite KG is located in `graphs/makaao_lite_kg.rdf`.

# MAKAAO Core dataset
This repository also includes the raw dataset collected from literature and reviewed by curators from which MAKAAO KG and MAKAAO Lite KG are populated. The dataset is located in the `data/makaao_core.csv`, and the normalized tables obtained from this raw dataset are located in the `data/processed_tables/` folder. This dataset includes some simple references to licensed elements like SNOMED CT, or UMLS CUI identifiers and names. Once the dataset is processed, these elements are also present in the resulting OWL files. Please keep in mind that eventhough this dataset underwent a preliminary review by our curators, it is currently still being reviewed by immunologists and experts of relevant pathologies. Thus, information it contains may be corrected in the future.

# Installation
0. If you are only interested in the pre-built knowledge graphs, you can simply download the relevant OWL files in the `graphs` folder.
1. If you want to run the graph-building script by yourself, you need a Python 3.11 environment. Other Python versions may or may not work, but there is no guarantee.
2. Clone the repository on you local machine:
   ```bash
   git clone https://github.com/f-maury/MAKAAO_core
   ```
3. Navigate to the repository, and install required packages:
   ```bash
   pip install -r requirements.txt
   ```
4. Run
   ```bash
   python scripts/build_graph/build_makaao_kg.py
   python scripts/build_graph/build_makaao_kg_lite.py
   ```
5. The newly-built RDF files will replace the pre-built ones, and you will find them in the `graphs` folder.
