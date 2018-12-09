# Description of data


### Definitions

**NCI**: National Cancer institute

**NSC Number**: Originally known as Cancer Chemotherapy National Service Center number, the NSC number is an identifying number assigned by DTP to an agent or product (e.g. small molecule or biological agent). Most structures have a unique NSC number but over the years a small percentage of structures or agents may have been assigned more than one NSC number. In the case of salts, different salt forms (e.g. HCl, HOAc) are designated with separate NSC numbers.
It can be used to search in the DTP (developmental therapeutics program) open chemical repository for ordering.


**SID**: identifier possible to look up in pubchem

**Diversity Set**: Compounds selected for dissimilar structures (using a pharmacore based method)

**.sdf**: file type (*s*tructure-*d*ata *f*ile) that can be read by RDKit

### Datafiles

* `Chem2D_Jun2016.zip` can be unzipped to reveal the full Open compound dataset
* `DIV5_*` are for the Diversity Set V
* `NSC_PubChemSID.csv` contains the mapping from NSC ids to Pubchem-searchable SID's

### Links:
https://dtp.cancer.gov/databases_tools/data_search_instructions.htm

https://wiki.nci.nih.gov/display/NCIDTPdata/Chemical+Data

https://wiki.nci.nih.gov/display/NCIDTPdata/Compound+Sets

https://dtp.cancer.gov/organization/dscb/obtaining/default.htm

https://pubchem.ncbi.nlm.nih.gov/


### Some helpful commands

```
python main.py --corpus /data/nci_open_training_data --embedding_size 256 --output_dir /data/graph2vec_embeddings/ --epochs 10
```
