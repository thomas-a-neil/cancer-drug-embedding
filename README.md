# Description of data


### Definitions

**NCI**: National Cancer institute

**NSC Number**: Originally known as Cancer Chemotherapy National Service Center number, the NSC number is an identifying number assigned by DTP to an agent or product (e.g. small molecule or biological agent). Most structures have a unique NSC number but over the years a small percentage of structures or agents may have been assigned more than one NSC number. In the case of salts, different salt forms (e.g. HCl, HOAc) are designated with separate NSC numbers.
It can be used to search in the DTP (developmental therapeutics program) open chemical repository for ordering.


**SID**: identifier possible to look up in pubchem

**Diversity Set**: Compounds selected for dissimilar structures (using a pharmacore based method)

**.sdf**: file type (*s*tructure-*d*ata *f*ile) that can be read by RDKit

**PASS** The program PASS (Prediction of Activity Spectra for Substances) was used to calculate predictions for up to 565 different activities for nearly all the structures in the database. PASS calculates the probability for both activity and inactivity of the compound for a given mechanism. These comprise specific enzymatic inhibitory potencies, therapeutic uses for various diseases, toxicities, and others. Counting the activity and inactivity predictions separately (they can be searched for separately), a total of 64,188,212 predicted values are offered on this site. Because the training set that underlies PASS is large but still limited (on the order of 35,000 compounds), the program cannot reliably predict each activity for every compound in the database.

### Datafiles

* `Chem2D_Jun2016.zip` can be unzipped to reveal the full Open compound dataset
* `DIV5_*` are for the Diversity Set V
* `NSC_PubChemSID.csv` contains the mapping from NSC ids to Pubchem-searchable SID's
* `all_nci_pass_data.sdf` contains all the feature data available from the cactus link below


### Links:
https://dtp.cancer.gov/databases_tools/data_search_instructions.htm

https://wiki.nci.nih.gov/display/NCIDTPdata/Chemical+Data

https://wiki.nci.nih.gov/display/NCIDTPdata/Compound+Sets

https://dtp.cancer.gov/organization/dscb/obtaining/default.htm

https://pubchem.ncbi.nlm.nih.gov/

https://cactus.nci.nih.gov/ncidb2.2/

### Some helpful commands

```
python main.py --corpus /data/nci_open_training_data --embedding_size 256 --output_dir /data/graph2vec_embeddings/ --epochs 10
```

graph2vec can be found at https://github.com/MLDroid/graph2vec_tf
