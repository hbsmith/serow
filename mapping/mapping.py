from Bio.KEGG import REST
import copy
import pickle
import pandas as pd
import collections
import json
import pprint
import itertools
import os

"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CORE MAPPING OUTLINE:

1 Reactions which can be directly associated with modules via the KEGG API
    1.1 Reactions found in any modules not containing addition rules
    1.2 Reactions found in any modules containing addition rules
    1.2.1 Reactions found in any modules containing addition rules, and only map to <=1 KOs
        1.2.2 ^^same as 1.2.1, but from modules added since the `module_entry_dict` was last downloaded
        1.2.3 Reactions found in any module containing addition rules and map to >1 KO 
2 Reactions which CANNOT be directly associated with modules via the KEGG API, but 
  CAN be directly associated with KOs via the KEGG API.
    2.1 Reactions that map to 1 KO (via direct rn->ko mappings)
    2.2 Reactions that map to >1 KOs (via direct rn->ko mappings)
3 Reactions which CANNOT be directly associated with modules via the KEGG API, AND
  CANNOT be directly associated with KOs via the KEGG API, but which CAN be directly
  associated to ECs that can be directly associated to KOs via the KEGG API
    3.1 Reactions that map to ECs that only map to 1 KO
    3.2 Reactions that map to ECs that map to >1 KO
4 Reactions which can happen spontaneously


NOTE: Mapping functions are those which produce an object of the final type of mapping output we're interested in:
      dict(Rn -> {frozenset(...), frozenset(...), ...})

NOTE: Every function which starts with "from_csv" is a mapping function.
      Every other mapping function starts with "map_rn_to_ko"

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FUNCTIONS MEANT TO BE CALLED BY USER:

FUNCTION PRE-MAPPING:
- write_links

WRITE-OUT FUNCTIONS FOR POST-MAPPING:
- dump_maps_by_type
- dump_maps_combined
- dump_maps_to_csv

CORE PIPELINE FUNCTIONS TO GENERATE MAPS:
- generate_all_maps
- write_manual_mapping_rsets_to_csv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
def create_a_to_b_dict(a,b):
    """Create dictionary mapping from KEGG database `a` entries to database `b` entries"""

    abbrv_dict = {"pathway":	"path",
                "brite":	"br",	
                "module":	"md",	
                "orthology":	"ko",	
                "genome":	"gn",	
                "compound":	"cpd",	
                "glycan":	"gl",	
                "reaction":	"rn",	
                "rclass":	"rc",	
                "enzyme":	"ec",	
                "network":	"ne",	
                "variant":	"hsa_var",
                "disease":	"ds",	
                "drug":	"dr",	
                "dgroup":	"dg"}

    ab_link = REST.kegg_link(a,b)
    ab_list = [i.split() for i in ab_link]
    ab_dict = dict()
    for pair in ab_list:
        if pair[0].startswith(abbrv_dict[a]):
            aid = pair[0].split(abbrv_dict[a]+":")[1]
            bid = pair[1].split(abbrv_dict[b]+":")[1]
        else:
            aid = pair[1].split(abbrv_dict[a]+":")[1]
            bid = pair[0].split(abbrv_dict[b]+":")[1]
        if aid in ab_dict:
            ab_dict[aid].append(bid)
        else:
            ab_dict[aid] = [bid]

    for aid,blist in ab_dict.items():
        ab_dict[aid] = set(blist)
        
    return ab_dict

def serialize_sets(obj):
    """Internal function for dumping sets as lists for jsons"""
    if isinstance(obj, set):
        return list(obj)

    return obj

def write_links(links_path = "../mydata/links"):
    """
    Write the following a_to_b dicts used in mapping:
        mo_to_ko_dict
        mo_to_rn_dict
        rn_to_ko_dict
        rn_to_mo_dict
    """

    mo_to_ko_dict = create_a_to_b_dict("module","orthology")
    mo_to_rn_dict = create_a_to_b_dict("module","reaction")
    rn_to_ko_dict = create_a_to_b_dict("reaction","orthology")
    rn_to_mo_dict = create_a_to_b_dict("reaction","module")

    with open(os.path.join(links_path,'mo_to_ko_dict.json'), 'w') as f:
        json.dump(mo_to_ko_dict, f, default=serialize_sets)

    with open(os.path.join(links_path,'mo_to_rn_dict.json'), 'w') as f:
        json.dump(mo_to_rn_dict, f, default=serialize_sets)

    with open(os.path.join(links_path,'rn_to_ko_dict.json'), 'w') as f:
        json.dump(rn_to_ko_dict, f, default=serialize_sets)

    with open(os.path.join(links_path,'rn_to_mo_dict.json'), 'w') as f:
        json.dump(rn_to_mo_dict, f, default=serialize_sets)


class Mapping:
    """
    Main class used to do mapping
    
    Module and reaction entry data is expected to be formatted as retrieved from Bio.TogoWS.entry format="json".
    Module entries required to have `orthologs` and `definiton` fields.
    Reaction entries required to have `orthologs` and `comment` fields.
    Links are expected to be formatted like the output of `write_links`.

    :param module_entries_path: path to downloaded module entries file 
    :param reaction_entries_path: path to downloaded reaction entries file 
    :param links_path: path to the downloaded links directory
    """

    def __init__(self,
                 module_entries_path = "module_ko_to_rn/assets/module_entry_dict-2021_03_22.json",
                 reaction_entries_path = "../mydata/reaction.json",
                 links_path = "../mydata/links"):

        ## KEGG entries
        self.module_entry_dict = self.load_json_into_dict(module_entries_path)
        self.reaction_entry_dict = self.load_json_into_dict(reaction_entries_path)

        ## Links inferred from KEGG
        self.mo_to_ko_dict = {k:set(vlist) for k,vlist in self.load_json_into_dict(os.path.join(links_path,"mo_to_ko_dict.json")).items()}
        self.mo_to_rn_dict = {k:set(vlist) for k,vlist in self.load_json_into_dict(os.path.join(links_path,"mo_to_rn_dict.json")).items()}
        self.rn_to_ko_dict = {k:set(vlist) for k,vlist in self.load_json_into_dict(os.path.join(links_path,"rn_to_ko_dict.json")).items()}
        self.rn_to_mo_dict = {k:set(vlist) for k,vlist in self.load_json_into_dict(os.path.join(links_path,"rn_to_mo_dict.json")).items()}

        ## Modules and reactions from modules with addition in definition (e.g. K00001+K00002)
        self.plusmodules, self.noplusmodules = self.set_plusmodules(self.module_entry_dict)
        self.rsetinplusmodules2plus = None ## Maps to at least one KO
        self.rsetinplusmodules1minus = None ## Maps to at least one KO
        self.get_rsets_in_plusmodules() ## Sets above two vars

        self.rn_to_ko_not_in_modules = None

        ## Finalized Reaction:Frozenset(KOs) dicts
        self.maps = dict()
        self.maps["map_rn2ko_viaMO_noaddition"] = None
        self.maps["map_rn2ko_viaMO_addition1minus"] = None
        self.maps["map_rn2ko_viaMO_addition1minus_extras"] = None
        self.maps["map_rn2ko_viaMO_addition2plus"] = None
        self.maps["map_rn2ko_viaKO_1minus"] = None
        self.maps["map_rn2ko_viaKO_2plus"] = None
        self.maps["map_rn2ko_viaEC_1minus"] = None
        self.maps["map_rn2ko_viaEC_2plus"] = None
        self.maps["map_rn2ko_spontaneous"] = None

        ## Write CSVs that need to be manually mapped
        self.csvs_need_mapping = dict()
        self.csvs_need_mapping["rn2ko_viaMO_addition2plus"] = None
        self.csvs_need_mapping["rn2ko_viaKO_2plus"] = None
        self.csvs_need_mapping["to_csv_rn2ko_viaEC_1minus"] = None ## This doesn't really need mapping, but when we mapped we pulled from a spreadsheet josh made
        self.csvs_need_mapping["to_csv_rn2ko_viaEC_2plus"] = None

        ## Read CSVs that were manually mapped
        self.csvs_mapped = dict()

    @staticmethod
    def load_json_into_dict(path):
        """Convenience function to load json into dict"""
        with open(path) as f:
            mydict = json.load(f)
        return mydict

    def parse_spreadsheet_rules(self,rns_rules,rn_col="Reaction",rule_col="Rule",verbose=False):
        """Read manually mapped spreadsheet into dict which can be used to make df"""
        rules_formatted = dict()
        for d in rns_rules:
            try:
                rules_formatted[d[rn_col]] = self.parse_rule(d[rule_col])
            except: 
                if verbose:
                    print(d)
        return rules_formatted

    @staticmethod
    def parse_rule(rule):
        """Parse a single rule from the manually mapped spreadsheet into sets of frozensets"""
        or_rules = rule.split(",")
        all_ands = set()
        for r in or_rules:
            all_ands.add(frozenset(r.split("+")))
        return set(all_ands)  

    @staticmethod
    def set_plusmodules(module_entry_dict):
        """Identify modules with and without addition in definition (e.g. K00001+K00002)"""
        noplusmodules = dict() ## This is actually only no plus in modules now
        for m, v in module_entry_dict.items():
            if "+" not in v["definition"]:
                noplusmodules[m] = copy.copy(v["definition"])

        plusmodules = set(module_entry_dict) - set(noplusmodules)
        return plusmodules, noplusmodules

    def get_rsets_in_plusmodules(self):
        """Identify reactions that map to 2+ KOs, and 1 KO, in modules that contain addition"""

        ## Set of all reactions that are in modules with `+`; Not required to map to a KO
        rsetinplusmodules = set()
        for m in self.plusmodules:
            if m in self.mo_to_rn_dict:
                rsetinplusmodules |= self.mo_to_rn_dict[m]

        ## Set of all reactions that are in modules with `+` AND which map to 2 or more KOs
        ## NOTE: These 2 or more KOs don't necessarily have to both be associated with the module, 
        ##       as the code is currently written I think this is fine, because we assumed these 
        ##       were isozymes and manually mapped these reactions to KOs anyways.
        rsetinplusmodules1minus = set()
        rsetinplusmodules2plus = set()
        for r in rsetinplusmodules:
            if r in self.rn_to_ko_dict:
                if len(self.rn_to_ko_dict[r])>1:
                    rsetinplusmodules2plus.add(r)
                else:
                    rsetinplusmodules1minus.add(r)

        self.rsetinplusmodules2plus = rsetinplusmodules2plus ## Maps to at least 1 KO
        self.rsetinplusmodules1minus = rsetinplusmodules1minus ## Maps to at least 1 KO

    ##########################################################################################
    ## CORE MAPPING FUNCTIONS
    ##########################################################################################
    ## 1 Reactions which can be directly associated with modules via the KEGG API
    ##########################################################################################
    ## 1.1 Reactions found in any modules not containing addition rules
    ##########################################################################################
    def map_rn2ko_viaMO_noaddition(self):
        """
        1.1 Reactions found in any modules not containing addition rules
        """
        rn_to_ko_noplusmodules = dict()
        rn_to_ko_noplusmodules_koinmodule = dict()

        for m in self.noplusmodules:
            if m in self.mo_to_rn_dict:
                for r in self.mo_to_rn_dict[m]:
                    if r in self.rn_to_ko_dict:
                        if r in rn_to_ko_noplusmodules:
                            rn_to_ko_noplusmodules[r] |= copy.copy(self.rn_to_ko_dict[r])
                        else:
                            rn_to_ko_noplusmodules[r] = copy.copy(self.rn_to_ko_dict[r])
                        
                        if r in rn_to_ko_noplusmodules_koinmodule:
                            rn_to_ko_noplusmodules_koinmodule[r] |= copy.copy({k for k in self.rn_to_ko_dict[r] if k in self.mo_to_ko_dict[m]})
                        else:
                            rn_to_ko_noplusmodules_koinmodule[r] =  {k for k in self.rn_to_ko_dict[r] if k in self.mo_to_ko_dict[m]}

        rn_to_ko_noplusmodules_koinmodule = {r:kset for r,kset in rn_to_ko_noplusmodules_koinmodule.items() if len(kset)>0} ## Make sure at least 1 ko
        self.maps["map_rn2ko_viaMO_noaddition"] = {r:{frozenset([k]) for k in kset} for r,kset in rn_to_ko_noplusmodules_koinmodule.items()} ## Format with frozensets

    ##########################################################################################
    ## 1.2 Reactions found in any modules containing addition rules
    ##########################################################################################
    ## 1.2.1 Reactions found in any modules containing addition rules, and only map to <=1 KOs
    ## 1.2.2 ^^same as 1.2.1, but from modules added since the `module_entry_dict` was last downloaded
    ## 1.2.3 Reactions found in any module containing addition rules and map to >1 KO 
    
    def map_rn2ko_viaMO_addition1minus(self):
        """
        1.2.1 Reactions found in any modules containing addition rules, and only map to <=1 KOs
        """
        addition_modules_rules_1minus = dict()
        for r in self.rsetinplusmodules1minus:
            addition_modules_rules_1minus[r] = set()
            if r in self.rn_to_ko_dict:
                for k in self.rn_to_ko_dict[r]:
                    addition_modules_rules_1minus[r].add(frozenset([k]))

        self.maps["map_rn2ko_viaMO_addition1minus"] = addition_modules_rules_1minus ## rn to frozenset(kos)


    def map_rn2ko_viaMO_addition1minus_extras(self):
        """
        1.2.2 same as 1.2.1, but from modules added since the `module_entry_dict` was last downloaded

        The modules should be double checked to make sure they don't contain addition, but we do 
        check to make sure the reactions don't directly link to multiple KOs.
        """
        totalmodulemappedset = {r for m in self.module_entry_dict if m in self.mo_to_rn_dict for r in self.mo_to_rn_dict[m]}
        rns_with_ko_and_mo_mapping = {r for r in self.rn_to_mo_dict if r in self.rn_to_ko_dict}

        addition_modules_rules_1minus_extras = dict()
        for r in (rns_with_ko_and_mo_mapping - totalmodulemappedset):
            if len(self.rn_to_ko_dict[r])==1:
                addition_modules_rules_1minus_extras[r] = set()
                addition_modules_rules_1minus_extras[r].add(frozenset(list(self.rn_to_ko_dict[r])))

        ## These were manually verified to contain no `+` signs in the definitions
        self.maps["map_rn2ko_viaMO_addition1minus_extras"] = addition_modules_rules_1minus_extras

    def to_csv_rn2ko_viaMO_addition2plus(self,OUTPATH):
        """ 
        1.2.3 (PRE-MAPPING): Reactions found in any module containing addition rules and map to >1 KO 

        These reactions require MANUAL ANNOTATION

        :param OUTPATH: path to write CSV which is used to manually verify mapping.
        """
        listdf = []
        for r in self.rsetinplusmodules2plus:
            for m in self.rn_to_mo_dict[r]:
                _dict = dict()
                _dict["Module"] = m
                _dict["Reaction"] = r
                _dict["url"] = "https://www.genome.jp/dbget-bin/www_bget?rn:"+r
                _dict["KOs"] = ",".join(self.rn_to_ko_dict[r])
                if r in self.reaction_entry_dict:
                    _dict["reaction_orthology"] = self.reaction_entry_dict[r]["orthologs"]
                if m in self.module_entry_dict:
                    _dict["module_orthology"] = self.module_entry_dict[m]["orthologs"]
                    _dict["module_definition"] = self.module_entry_dict[m]["definition"]
                listdf.append(_dict)

        rn2ko_viaMO_addition2plus = pd.DataFrame(listdf)
        rn2ko_viaMO_addition2plus.sort_values(by=["Reaction","Module"],inplace=True)
        rn2ko_viaMO_addition2plus.reset_index(drop=True,inplace=True)

        # OUTPATH = "reaction_mapping/module_ko_to_rn/assets/rsetinplusmodules2plusAllModules.csv"
        self.csvs_need_mapping["rn2ko_viaMO_addition2plus"] = rn2ko_viaMO_addition2plus
        rn2ko_viaMO_addition2plus.to_csv(OUTPATH,index=False)

    def from_csv_rn2ko_viaMO_addition2plus(self,path):
        """ 
        1.2.3 (POST-MAPPING): Reactions found in any module containing addition rules and map to >1 KO 

        These reactions were MANUALLY ANNOTATED

        :param path: path to manually mapped CSV
        """
        spreadsheet = pd.read_csv(path, index_col=False)
        rns_rules = spreadsheet[["Reaction","Rule"]].to_dict(orient="records")
        self.maps["map_rn2ko_viaMO_addition2plus"] = self.parse_spreadsheet_rules(rns_rules)

    ##########################################################################################
    ## 2 Reactions which CANNOT be directly associated with modules via the KEGG API, but
    ##   CAN be directly associated with KOs via the KEGG API.
    ##########################################################################################
    ## 2.1 Reactions that map to 1 KO (via direct rn->ko mappings)
    ## 2.2 Reactions that map to >1 KOs (via direct rn->ko mappings)

    def get_rn_to_ko_not_in_modules(self):
        rns_not_in_modules = set(self.rn_to_ko_dict) - set(self.rn_to_mo_dict)
        self.rn_to_ko_not_in_modules = {k:self.rn_to_ko_dict[k] for k in rns_not_in_modules} ## rns_to_kos_not_in_modules
    
    def map_rn2ko_viaKO_1minus(self):
        """
        2.1 Reactions that only map to 1 KO (via direct rn->ko mappings)
        """
        self.get_rn_to_ko_not_in_modules()
        self.maps["map_rn2ko_viaKO_1minus"] = {k:{frozenset(v)} for k,v in self.rn_to_ko_not_in_modules.items() if len(v)==1}
        
    def to_csv_rn2ko_viaKO_2plus(self, OUTPATH):
        """
        2.2 (PRE-MAPPING) Reactions that map to >1 KOs (via direct rn->ko mappings)

        These reactions require MANUAL ANNOTATION

        :param OUTPATH: path to write CSV which is used to manually verify mapping.        
        """
        ## NOTE: Emulates what I think Josh did to make the spreadsheet 'OLD-keggReactionsToKOs.ambig', 
        ##       which is the sheet I then added keywords to.
        ## TODO: Add keyword as part of this initial function

        self.get_rn_to_ko_not_in_modules()
        listdf = []
        
        for r, kos in {k:v for k,v in self.rn_to_ko_not_in_modules.items() if len(v)>1}.items():  
            ofield = self.reaction_entry_dict[r]["orthologs"]
            _kw_dict = copy.copy({"effector":False,
                                "anchor":False,
                                "auxiliary":False,
                                "carrier":False,
                                "pts":False,
                                "subunit":False,
                                "chain":False,
                                "reductase component":False})
            _dict = dict()
            _dict["Reaction"] = r
            _dict["url"] = "https://www.genome.jp/dbget-bin/www_bget?rn:"+r
            _dict["KOs"] = ",".join(kos)

            for k,txt in ofield.items():
                txtlower = txt.lower()

                for keyword in _kw_dict:
                    if keyword in txtlower:
                        _kw_dict[keyword] = True

            _kw_dict["orthology"] = pprint.pformat(ofield) ## Format orthology
            _dict |= _kw_dict ## Add _kw_dict to _dict
            listdf.append(_dict)

        rn2ko_viaKO_2plus = pd.DataFrame(listdf)
        rn2ko_viaKO_2plus.sort_values(by=["Reaction"],inplace=True)
        rn2ko_viaKO_2plus.reset_index(drop=True,inplace=True)
        self.csvs_need_mapping["rn2ko_viaKO_2plus"] = rn2ko_viaKO_2plus
        rn2ko_viaKO_2plus.to_csv(OUTPATH,index=False)

    def from_csv_rn2ko_viaKO_2plus(self, path):
        """
        2.2 (POST-MAPPING) Reactions that map to >1 KOs (via direct rn->ko mappings)

        These reactions were MANUALLY ANNOTATED

        :param path: path to manually mapped CSV       
        """

        spreadsheet = pd.read_csv(path,index_col=False)
        rns_rules = spreadsheet[["Reaction","Rule"]].to_dict(orient="records")
        self.maps["map_rn2ko_viaKO_2plus"] = self.parse_spreadsheet_rules(rns_rules)    

    ##########################################################################################
    ## 3 Reactions which CANNOT be directly associated with modules via the KEGG API, AND
    ##   CANNOT be directly associated with KOs via the KEGG API, but which CAN be directly
    ##   associated to ECs that can be directly associated to KOs via the KEGG API
    ##########################################################################################
    ## 3.1 Reactions that map to ECs that only map to 1 KO
    ## 3.2 Reactions that map to ECs that map to >1 KO
    ## NOTE: THIS DATA WAS ORIGINALLY GENERATED BY JOSH, SO I DON'T HAVE THE SCRIPT WE USED TO MAKE THESE SHEETS.
    def to_csv_rn2ko_viaEC_1minus(self):
        """
        3.1 Reactions that map to ECs that only map to 1 KO
        """
        print("This was done by Josh")

    def from_csv_rn2ko_viaEC_1minus(self, path): ## AKA resolved
        """
        3.1 Reactions that map to ECs that only map to 1 KO
        """
        spreadsheet = pd.read_csv(path, index_col=False)
        rns_rules = spreadsheet[["rn","ko_list"]].to_dict(orient="records")
        self.maps["map_rn2ko_viaEC_1minus"] = self.parse_spreadsheet_rules(rns_rules,rn_col="rn",rule_col="ko_list")

    def to_csv_rn2ko_viaEC_2plus(self):
        """
        3.2 (PRE-MAPPING) Reactions that map to ECs that map to >1 KO
        
        These reactions require MANUAL ANNOTATION
        """
        print("This was done by Josh")

    def from_csv_rn2ko_viaEC_2plus(self, path): ## AKA ambiguous
        """
        3.2 (POST-MAPPING) Reactions that map to ECs that map to >1 KO
        
        These reactions were MANUALLY ANNOTATED

        :param path: path to manually mapped CSV   
        """
        spreadsheet = pd.read_csv(path, index_col=False)
        spreadsheet = spreadsheet[spreadsheet.ko_list != "NA"]
        rns_rules = spreadsheet[["rn","rule"]].to_dict(orient="records")
        self.maps["map_rn2ko_viaEC_2plus"] = self.parse_spreadsheet_rules(rns_rules,rn_col="rn",rule_col="rule")

    ##########################################################################################
    ## 4 Reactions which can happen spontaneously
    ##########################################################################################
    def map_rn2ko_spontaneous(self):
        """
        4 Reactions which can happen spontaneously

        These may have any of the previous mappings). 
        Identified by looking through comments which contain any of the following keywords (case-insensitive):
            "spontaneous"
            "non enzymatic"
            "non-enzymatic"
            "nonenzymatic"
            
        List was manually double checked to make sure none of the comments erroneously identified reactions 
        which may have said something like, "not spontanteous" (no instances of this were found and our full 
        automatically identified list was used). 
        """
        keyword_rn_dict = {"spontaneous":dict(),
                    "non enzymatic":dict(),
                    "non-enzymatic":dict(),
                    "nonenzymatic":dict(),
                    "multi step":dict(),
                    "multi-step":dict(),
                    "multistep":dict()}
        
        for r,v in self.reaction_entry_dict.items():
            comment = v["comment"].lower()   
            for k in keyword_rn_dict:
                if k in comment:
                    keyword_rn_dict[k][r] = copy.copy(comment)

        spontaneous_rns = {**keyword_rn_dict["spontaneous"],
            **keyword_rn_dict["non enzymatic"],
            **keyword_rn_dict["non-enzymatic"],
            **keyword_rn_dict["nonenzymatic"]}

        for r in spontaneous_rns:
            spontaneous_rns[r] = {frozenset(["spontaneous"])}

        self.maps["map_rn2ko_spontaneous"] = spontaneous_rns

    ##########################################################################################    
    ## BASIC STATISTICS
    ##########################################################################################
    def describe_nfrozensets_in_maps(self):
        rnsetsum = 0
        for d,self.rn_to_ko_dict in self.maps.items():
            print(d)
            print(len(self.rn_to_ko_dict))
            print(collections.Counter([len(v) for k,v in self.rn_to_ko_dict.items()]))
            rnsetsum+=len(self.rn_to_ko_dict)
        print("rnsetsum")
        print(rnsetsum)

    def describe_nfrozensets_and_sizefrozensets_in_maps(self):
        for d, self.rn_to_ko_dict in self.maps.items():
            print(d)
            print(collections.Counter([len(v) for k,v in self.rn_to_ko_dict.items()]).most_common())
            metacounter = dict()
            for k,v in self.rn_to_ko_dict.items():        
                for c in collections.Counter([len(fs) for fs in v]).items():
                    if c[0] in metacounter:
                        metacounter[c[0]]+=c[1]
                    else:
                        metacounter[c[0]]=c[1]
            print(metacounter)
            print("---------")

    def describe_overlap_of_maps(self):
        for dpair in list(itertools.combinations(self.maps,2)):
            print(dpair)
            print(len(set(self.maps[dpair[0]]) & set(self.maps[dpair[1]])))

    ########################################################################################
    ## WRITE-OUT FUNCTIONS FOR POST-MAPPING
    ########################################################################################
    @staticmethod
    def dump_maps_by_type(maps, OUTPATH = "final_map/rn2ko_map_by_type.pkl"):
        pickle.dump(maps, open(OUTPATH, 'wb'))

    @staticmethod
    def dump_maps_combined(maps, OUTPATH = "final_map/rn2ko_map_combined.pkl"):
        ## Combined
        combined_dict_of_rules = dict()
        for d, rn_to_ko_dict in maps.items():
            for r, kset in rn_to_ko_dict.items():
                if r in combined_dict_of_rules:
                    combined_dict_of_rules[r] |= copy.copy(kset)
                else:
                    combined_dict_of_rules[r] = copy.copy(kset)

        pickle.dump(combined_dict_of_rules, open(OUTPATH, 'wb'))

    @staticmethod
    def dump_maps_to_csv(maps, OUTPATH = "final_map/rn2ko_map_by_type.csv"):
        ## Spreadsheet; could clean up using pd.pipe maybe
        rules_df = pd.DataFrame.from_dict(maps).reset_index()
        rules_df.rename(columns={"index": "reaction"},inplace=True)
        rules_df = pd.melt(rules_df,ignore_index=False,id_vars=["reaction"],value_vars=list(maps.keys()),var_name="origin")
        rules_df.dropna(inplace=True)
        rules_df = rules_df.explode("value")
        rules_df.rename(columns={"value": "rule"},inplace=True)
        rules_df.reset_index(drop=True,inplace=True)
        rules_df.sort_values(by=["origin","reaction"],inplace=True)
        rules_df.reset_index(drop=True,inplace=True)
        # OUTPATH = "reaction_mapping/combined_rn_to_ko_rules-MOFILTERFIX.csv"
        rules_df.to_csv(OUTPATH,index=False)


    ########################################################################################
    ## CORE PIPELINE FUNCTIONS TO GENERATE MAPS
    ########################################################################################
    def generate_all_maps(self,
                          rn2ko_viaMO_addition2plus_path = "reaction_ko_keywords/KEGG.ReactionsToKOs.Ambiguous.12July2021+module_rns.csv",
                          rn2ko_viaKO_2plus_path = "reaction_ko_keywords/KEGG.ReactionsToKOs.Ambiguous.12July2021-NEWER2.csv",
                          rn2ko_viaEC_1minus_path = "reaction_ko_keywords/KEGG.ReactionsToKOs.Ambiguous.12July2021-resolved.rn2ko.viaEC.csv",
                          rn2ko_viaEC_2plus = "reaction_ko_keywords/KEGG.ReactionsToKOs.Ambiguous.12July2021-ambig.rn2ko.viaEC.csv"):#,reaction_json_path, rsetinplusmodules2plus_to_csv_OUTPATH):
        """Map from all sources, including manually mapped CSVs"""
        self.map_rn2ko_viaMO_noaddition()
        self.map_rn2ko_viaMO_addition1minus()
        self.map_rn2ko_viaMO_addition1minus_extras()
        self.from_csv_rn2ko_viaMO_addition2plus(rn2ko_viaMO_addition2plus_path)
        self.map_rn2ko_viaKO_1minus()
        self.from_csv_rn2ko_viaKO_2plus(rn2ko_viaKO_2plus_path)
        self.from_csv_rn2ko_viaEC_1minus(rn2ko_viaEC_1minus_path)
        self.from_csv_rn2ko_viaEC_2plus(rn2ko_viaEC_2plus)
        self.map_rn2ko_spontaneous()

        for k,v in self.maps.items():
            if v != None:
                print(k, len(v))

        return self.maps

    def write_manual_mapping_rsets_to_csv(self,
                                          rn2ko_viaMO_addition2plus_path = "csvs_need_mapping/rn2ko_viaMO_addition2plus.csv",
                                          rn2ko_viaKO_2plus_path = "csvs_need_mapping/rn2ko_viaKO_2plus.csv",
                                          rn2ko_viaEC_1minus_path = "csvs_need_mapping/rn2ko_viaEC_1minus.csv",
                                          rn2ko_viaEC_2plus_path = "csvs_need_mapping/rn2ko_viaEC_2plus.csv"):
        """Write all CSVs which will need manual mapping"""
        self.to_csv_rn2ko_viaMO_addition2plus(rn2ko_viaMO_addition2plus_path)
        self.to_csv_rn2ko_viaKO_2plus(rn2ko_viaKO_2plus_path)
        self.to_csv_rn2ko_viaEC_1minus(rn2ko_viaEC_1minus_path)
        self.to_csv_rn2ko_viaEC_2plus(rn2ko_viaEC_2plus_path)

if __name__ == "__main__":

    # write_links()

    map = Mapping()
    maps = map.generate_all_maps()
    # map.write_manual_mapping_rsets_to_csv()
    map.dump_maps_by_type(maps, OUTPATH="final_map/v3/rn2ko_map_by_type.pkl")
    map.dump_maps_combined(maps, OUTPATH="final_map/v3/rn2ko_map_combined.pkl")
    map.dump_maps_to_csv(maps, OUTPATH="final_map/v3/rn2ko_map_by_type.csv")
