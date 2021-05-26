README.txt

############################################
OVERVIEW
############################################

This script is meant to parse through and reformat the logic within a KEGG module's "Defintion" field.
It can return a dictionary keyed by Module IDs, the value of each is a dictionary of reactions.
These are the reactions associated with the module.
The reactions contain sets of frozensets of KO IDs.
Each KO within each frozenset is required to catalyze the keyed reaction.
Each frozenset within each set is an alternative way to catalyze the reaction.

See for example: https://www.kegg.jp/dbget-bin/www_bget?M00009

For example, here's the value of the key "M00009":

    {'R00267': {frozenset({'K00031'}), frozenset({'K00030'})},
    'R00268': {frozenset({'K00031'}), frozenset({'K00030'})},
    'R00342': {frozenset({'K00024'}),
                frozenset({'K00025'}),
                frozenset({'K00026'})},
    'R00351': {frozenset({'K05942'}), frozenset({'K01647'})},
    'R00361': {frozenset({'K00116'})},
    'R00405': {frozenset({'K01903', 'K01902'}), frozenset({'K01899', 'K01900'})},
    'R00432': {frozenset({'K01899', 'K01900'})},
    'R00621': {frozenset({'K01616', 'K00382'}),
                frozenset({'K00164', 'K00382', 'K00658'})},
    'R00709': {frozenset({'K00031'}), frozenset({'K00030'})},
    'R00727': {frozenset({'K01899', 'K01900'})},
    'R01082': {frozenset({'K01678', 'K01677'}),
                frozenset({'K01676'}),
                frozenset({'K01679'})},
    'R01197': {frozenset({'K00174', 'K00175'}),
                frozenset({'K00176', 'K00174', 'K00175'}),
                frozenset({'K00174', 'K00175', 'K00177'}),
                frozenset({'K00176', 'K00174', 'K00175', 'K00177'})},
    'R01324': {frozenset({'K01682'}), frozenset({'K01681'})},
    'R01325': {frozenset({'K01682'}), frozenset({'K01681'})},
    'R01899': {frozenset({'K00031'}), frozenset({'K00030'})},
    'R01900': {frozenset({'K01682'}), frozenset({'K01681'})},
    'R02164': {frozenset({'K00235', 'K00237', 'K00234', 'K00236'}),
                frozenset({'K00241', 'K00242', 'K00240', 'K00239'}),
                frozenset({'K18860', 'K00240', 'K00241', 'K00239'}),
                frozenset({'K00241', 'K00240', 'K00239', 'K18859'}),
                frozenset({'K00245', 'K00246', 'K00244'}),
                frozenset({'K00245', 'K00246', 'K00244', 'K00247'}),
                frozenset({'K00241', 'K00240', 'K00239'})},
    'R02570': {frozenset({'K01616', 'K00382'}),
                frozenset({'K00164', 'K00382', 'K00658'})},
    'R03316': {frozenset({'K01616', 'K00382'}),
                frozenset({'K00164', 'K00382', 'K00658'})},
    'R07618': {frozenset({'K01616', 'K00382'}),
                frozenset({'K00164', 'K00382', 'K00658'})},
    'R10343': {frozenset({'K18118'})}}


Within M00009, there are two frozen sets of KOs which allow reaction 'R07618' to be catalyzed:
frozenset({'K01616', 'K00382'}), and frozenset({'K00164', 'K00382', 'K00658'}).

This reaction requires either 1 or 2:
1. 'K01616' and 'K00382'
2. 'K00164' and 'K00382' and 'K00658'

Rules for formatting determined based on the KEGG "Help" link, and based on representations in 
the provided diagrams. Help text below:

    The definition of the module as a list of K numbers for pathway/signature modules and RC numbers 
    for reaction modules. Comma separated K numbers or RC numbers indicate alternatives. Plus signs 
    are used to represent a complex or a combination and a minus sign denotes a non-essential 
    component in the complex.

The script can also provide a combined set of rules for each reaction across all modules. This is 
basically just a union of all sets associated with each reaction.

Modules whose definitions are made up of other modules are not included.

############################################
USAGE
############################################

There are only a few functions meant to be called by the user.

collect_KEGG_files
    Retrieves KEGG module info from either:
    1. existing jsons specified from the paths
    2. the KEGG REST API, and then writes them to the paths specified.
    Currently rn_ko_dict, and ko_rn_dict are unused.

parse_and_format_modules
    WRITES TO "assets/calculated_module_dict.pkl" BY DEFAULT.
    This is the meat of the whole script, encapsulating the parsing and 
    reformatting of the defintions. This gives you a ruleset of KOs for
    each module.

    Returns a dictionary of:
    keys=modules IDs
    values=sets of frozen sets of valid KO combinations for steps of the module


create_dict_of_local_r_to_k_rules
    WRITES TO "assets/dict_of_local_r_to_k_rules.pkl" BY DEFAULT.
    This gives you a ruleset of KOs for each reaction, within each module.

    Returns a dictionary of:
        keys=modules IDs
        values=dictionary of reactions within the module 
        
        The reaction dictionaries are:
        keys=Reaction IDs
        values = sets of frozen sets of valid KO combinations for catalysis of the reaction
    

create_dict_of_global_r_to_k_rules
    WRITES TO "assets/dict_of_global_r_to_k_rules.pkl" BY DEFAULT.
    This gives you a ruleset of KOs for each reaction, regardless of module. It is a union
    for any given reaction's sets across all modules.

    Returns a dictionary of:
        keys=Reaction IDs
        values = sets of frozen sets of valid KO combinations for catalysis of the reaction
