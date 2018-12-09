import pandas as pd

div6_plate_maps_df = pd.read_csv('data/DIV6_Plate Maps_v2.csv')
div5_plate_maps_df = pd.read_csv('data/DIV5_platemap.csv')

div2_nsc_to_pubchem_sid_df = pd.read_csv('data/divii_mlsmr_nsc_to_pubchem_sid.csv', header=None, names=['NSC', 'SID'])
div2_nsc = div2_nsc_to_pubchem_sid_df['NSC'].values

nsc_pubchem_df = pd.read_csv('data/NSC_PubChemSID.csv', header=None, names=['NSC', 'SID'])

all_open_compounds_filepath = 'data/Chem2D_Jun2016.sdf'


# check diversity set 5 NSCs if they can be converted to
# pubmed SIDs
def check_div5_nscs():
    div5_nsc = set(div5_plate_maps_df['NSC'].values)
    pubchem_nsc = set(nsc_pubchem_df['NSC'].values)
    return div5_nsc < pubchem_nsc


def lookup_nsc(nsc):
    return b'<NSC>\r\n%b\r' % str(nsc).encode('utf8')


# check that all diversity 5 compounds are in the
# full open set
def check_div5_in_open():
    with open(all_open_compounds_filepath, 'rb') as f:
        all_open = f.read()

    def includes_nsc(nsc):
        lookup = lookup_nsc(nsc)
        return lookup in all_open

    div5_nsc = set(div5_plate_maps_df['NSC'].values)
    list_div5_nsc = list(div5_nsc)
    check_div5 = [includes_nsc(nsc) for nsc in list_div5_nsc]
    return all(check_div5)
