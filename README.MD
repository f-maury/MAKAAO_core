[![66.67 % FAIR](https://img.shields.io/badge/FAIR_assessment-66.67_%25-orange)](https://fair-checker.france-bioinformatique.fr/assessment/6841584a442568dc9750d7a7)
![Version](https://img.shields.io/badge/version-1.0.0-blue) 
![SHACL Validation](https://github.com/f-maury/MAKAAO_core/actions/workflows/shacl_validation.yml/badge.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18228827.svg)](https://doi.org/10.5281/zenodo.18228827)

# MAKAAO knowledge graph

## Installation

   1. Download the repository and enter it: 
   ```bash
   git clone https://github.com/f-maury/MAKAAO_core
   cd MAKAAO_core
   ```

   2. Create a Conda environment with the needed packages:
   ```bash
   conda create -n makaao python=3.11
   conda activate makaao
   python -m pip install -r requirements.txt
   ```
   
   3. Browse to `scripts/` folder, and run the MAKAAO scripts in correct order:
   ```bash
   cd scripts
   python 00_xlsx_to_csv.py
   ```
   ...

   4. If everything runs without issue:
      * MAKAAO KG and Lite KG will be located in `kg/` folder
      * Normalized tables will be in `data/processed_tables/`.
      * Enrichment tables will be in `data/enrichment_tables/`.

## Annex
* `tests/` contains the `shacl_shapes.ttl` file used by the GitHub CI workflow to validate some constraints on the knowledge graph.
* Sample data are located in `data/makaao_sample.csv`.


## WARNING

For licensing reasons, some files that you would need to correctly run the scripts are not included in this repository. These are:
   * `data/en_product4.xml`, obtained from Orphanet (https://www.orphadata.com/data/xml/en_product4.xml)
   * `data/MRCONSO.RRF`, obtained from the UMLS (https://www.nlm.nih.gov/research/umls/licensedcontent/umlsknowledgesources.html)
   * `data/makaao_core.xlsx`
