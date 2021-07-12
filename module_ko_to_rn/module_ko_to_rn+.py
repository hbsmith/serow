import os
import re
import json
import copy
import itertools
import pickle
import pyparsing as pp
from Bio.KEGG import REST #, Enzyme, Compound, Map
import Bio.TogoWS as TogoWS
from tqdm import tqdm
from pprint import pprint

##########################################
## Pyparsing definitions
##########################################
class Operation(object):
    def __init__(self, tokens):
        self._tokens = tokens[0]
        self.assign()

    def assign(self):
        """
        function to copy tokens to object attributes
        """

    def __repr__(self):
        return self.__class__.__name__ + ":" + repr(self.__dict__)
    __str__ = __repr__
    
class UnOp(Operation):
    def assign(self):
        self.op = self._tokens[0]
        self.terms = [self._tokens[1]]
        del self._tokens

class BinOp(Operation):
    def assign(self):
        self.op = self._tokens[1]
        self.terms = self._tokens[0::2]
        del self._tokens

class optionalOrOpUn(UnOp):
    pass

class optionalOrOpBin(BinOp):
    pass

class andOp(BinOp):
    pass

class mandatoryOrOp(BinOp):
    pass

class thenOp(BinOp):
    pass

##########################################
## Object that all parsed KOs lists become
##########################################
class ValidExprs(object):
    def __init__(self, expressions):
        self.expressions = expressions
    def __repr__(self): 
        return str(self._expressions)

    @property
    def expressions(self):
        """Valid Expressions"""
        return self._expressions

    @expressions.setter
    def expressions(self, value):
        if isinstance(value,str):
            self._expressions = [[value]]

        elif all(isinstance(elem, list) for elem in value):
            self._expressions = value

        elif isinstance(value, ValidExprs):
            self._expressions = value
            
        ### Won't throw if list of list of lists but better than nothing...
        else:
            raise ValueError("Not a valid expressions")

##########################################
## Grabbing/Linking KEGG data
##########################################
def create_id_name_dict(db):
    ## Grab list of ids in db
    id_name_dict = dict()
    raw_list = REST.kegg_list(db)
    id_name_list = [s.split('\t') for s in raw_list.read().splitlines()]
    for i in id_name_list:
        id_name_dict[i[0]] = i[1]
    return id_name_dict

def retrieve_entry_info(id_name_dict,db):
    ## Grab each entry in list of ids
    id_entry_dict = {}
    for entry in tqdm(id_name_dict.keys()):

        entry_id = entry.split(":")[1]
        
        try:
            id_entry_dict[entry_id]=json.load(TogoWS.entry(db, entry_id, format="json"))[0]
        except:
            pass

    return id_entry_dict

def create_link_dicts(target_db,source_db):
    ## Get KEGG ko-reaction mapping--this may be useful for verifying module reaction-kos are associated correctly
    ## The ordering of source and target doesn't matter for the content of the output, just the ordering of the output
    raw_links = REST.kegg_link(target_db,source_db)
    link_list = [s.split('\t') for s in raw_links.read().splitlines()]
    target_source_dict = dict()
    source_target_dict = dict()
    for row in link_list:
        s = row[0].split(":")[1]
        t = row[1].split(":")[1]
        if s in source_target_dict:
            source_target_dict[s].add(t)
        else:
            source_target_dict[s] = {t}
            
        if t in target_source_dict:
            target_source_dict[t].add(s)
        else:
            target_source_dict[t] = {s}
        
    return source_target_dict,target_source_dict

def collect_KEGG_files(db, entry_dict_path, source_target_link_path):
    ## Load module entry dict if it exists
    if os.path.isfile(entry_dict_path):
        with open(entry_dict_path, 'r') as f:
            module_entry_dict = json.load(f)
    else:
        module_name_dict = create_id_name_dict(db)
        module_entry_dict = retrieve_entry_info(module_name_dict,db)
        with open(entry_dict_path, 'w') as f:
            json.dump(module_entry_dict, f, indent=4)

    ## Load link dicts if they exist
    if os.path.isfile(source_target_link_path):
        with open(source_target_link_path, 'r') as f:
            source_target_links = json.load(f)
            rn_ko_dict = source_target_links[0]
            ko_rn_dict = source_target_links[1]

    else:
        rn_ko_dict, ko_rn_dict = create_link_dicts("ko","reaction")
        with open(source_target_link_path, 'w') as f:
            json.dump([rn_ko_dict, ko_rn_dict], f, indent=4, default=convert_sets_to_lists)

    return module_entry_dict, rn_ko_dict, ko_rn_dict


##########################################
## Functions to parse and format KEGG data
##########################################
def declareSearchExpr():
    ## Order of operations for parsing
    searchTerm = pp.Combine('K' + pp.Word(pp.nums, exact=5)).leaveWhitespace()

    searchExpr = pp.infixNotation(searchTerm | pp.quotedString.setParseAction(pp.removeQuotes),
        [
        ('-', 1, pp.opAssoc.RIGHT, optionalOrOpUn),
        (pp.White(exact=1), 2, pp.opAssoc.LEFT, thenOp),
        (pp.oneOf('+'), 2, pp.opAssoc.LEFT, andOp),
        (pp.oneOf('-'), 2, pp.opAssoc.LEFT, optionalOrOpBin),
        (pp.oneOf(','), 2, pp.opAssoc.LEFT, mandatoryOrOp),
        
        ])

    return searchExpr

def getTopLevelOp(module_def):
    ## List of mix of parseResultsObjs and KOs to eventually feed into object
    searchExpr = declareSearchExpr()
    return searchExpr.parseString(module_def).asList()[0]

def replaceStrsWithValidExprs(data):
    ## Replace innermost strings with validExprs (consistent object type to later operate on)
    if issubclass(type(data), Operation):
        data.terms = [replaceStrsWithValidExprs(term) for term in data.terms]
        return data
    elif isinstance(data, str):
        return ValidExprs(data)

def DoOp(validExprsList,op,data):
    if op==",":
        return ValidExprMandatoryOr(validExprsList)
    elif op==" ":
        return ValidExprAnd(validExprsList) ## Then is same as mandatory Or
    elif op=="-":
        if isinstance(data,optionalOrOpBin):
            return ValidExprOptionalOrBin(validExprsList)
        elif isinstance(data,optionalOrOpUn):
            return ValidExprOptionalOrUn(validExprsList)
        else:
            raise ValueError("Should be one of two OptionalOr classes")
    elif op=="+":
        return ValidExprAnd(validExprsList)

def ValidExprMandatoryOr(validExprList):
    return sum([ve.expressions for ve in validExprList], [])

def ValidExprOptionalOrUn(validExprList):
    return sum([ve.expressions for ve in validExprList], [[]]) ## Empty list possible too

def ValidExprOptionalOrBin(validExprList):

    newvelist = []
    for ve in validExprList[1:]:        
        newve = copy.copy(ve.expressions)
        newve.append([])
        newvelist.append(newve)
        
    newvelist.insert(0,validExprList[0].expressions)
    
        
    return [sum(i,[]) for i in itertools.product(*newvelist)]

def ValidExprAnd(validExprList):
    return [sum(i,[]) for i in itertools.product(*[ve.expressions for ve in validExprList])]

def combineValidExprs(data):
    ## Recursively go through parsed KOs (operators and operands) to reformat consistently
    if issubclass(type(data), Operation): 
        if all([isinstance(term,ValidExprs) for term in data.terms]): ## Don't transform things at higher levels
            data = ValidExprs(DoOp(data.terms,data.op,data)) ## Transform to ValidExprs type after operating
            return data
        else: 
            while not all([isinstance(term,ValidExprs) for term in data.terms]):
                data.terms = [combineValidExprs(term) for term in data.terms]
            return data

    elif isinstance(data, ValidExprs):
        return data

def moduleParseObjToExpressions(moduleThenObj):
    mNoStr = replaceStrsWithValidExprs(moduleThenObj)
    m2 = combineValidExprs(mNoStr)
    m1 = combineValidExprs(m2)
    return m1.expressions

def convert_sets_to_lists(obj):
    ## Used to write sets to JSONs
    if isinstance(obj, set):
        return list(obj)
    raise TypeError

def parse_and_format_modules(module_entry_dict):
    ## MAIN FUNCTION TO CALL FOR PARSING AND FORMATTING MODULES
    ## Takes the longest amount of time, on order of several minutes
    # count = 0
    calculated_module_dict = dict()
    for mid, d in module_entry_dict.items():
        if not len(re.findall(r'[M]\d{5}',d["definition"]))>0:
            # if mid=="M00004":
            print(mid)
            data = getTopLevelOp(d["definition"])
            data = replaceStrsWithValidExprs(data) ## Transforms strings to "ValidExprs" objects that can be typechecked
            data = combineValidExprs(data) ## Transforms all but the top level to ValidExprs objs
            data = combineValidExprs(data) ## Transforms top level (Returns a single ValidExprs obj)
            calculated_module_dict[mid] = {frozenset(i) for i in data.expressions} - {frozenset()}

            # pprint(calculated_module_dict[mid])

    # Pickle dictionary using protocol 3.
    pickle.dump(calculated_module_dict, open('assets/calculated_module_dict_AND.pkl', 'wb'))
    return calculated_module_dict

##########################################
## Link Module KOs to Reactions
##########################################
def get_r_to_k_rules(mid, module_entry_dict, calculated_module_dict):
    
    ## Parse Orthologs
    local_r_to_k_dict = dict()
    for k,v in module_entry_dict[mid]["orthologs"].items():
        kids = re.findall(r'[K]\d{5}',k)
        rids = re.findall(r'[R]\d{5}',v)
        for rid in rids:
            if rid not in local_r_to_k_dict:
                local_r_to_k_dict[rid] = copy.copy(kids)
            else:
                local_r_to_k_dict[rid] += kids

    ## Coerce into rule sets
    local_r_to_k_dict_validated = dict()
    ## If all KOs are independent, can take a little shortcut
    if all([len(i)==1 for i in calculated_module_dict[mid]]):
        print("shortcut")
        for rid, kids in local_r_to_k_dict.items():
            local_r_to_k_dict_validated[rid] = {frozenset([i]) for i in kids}

    else:
        print("one by one")
        for rid, kids in local_r_to_k_dict.items():

            kids_set = set(kids)

            ## If R only corresponds to one KO, don't have to check the calculated_module_dict
            if len(kids)==1:
                local_r_to_k_dict_validated[rid] = {frozenset([i]) for i in kids}

            else:

                for fs in calculated_module_dict[mid]:

                    if fs.issubset(kids_set):

                        if rid not in local_r_to_k_dict_validated:
                            local_r_to_k_dict_validated[rid] = {fs}
                        else:
                            local_r_to_k_dict_validated[rid].add(fs)

    return local_r_to_k_dict_validated

def create_dict_of_local_r_to_k_rules(calculated_module_dict, module_entry_dict):
    dict_of_local_r_to_k_rules = {}
    for mid in calculated_module_dict:
        dict_of_local_r_to_k_rules[mid] = get_r_to_k_rules(mid, module_entry_dict, calculated_module_dict)

    pickle.dump(dict_of_local_r_to_k_rules, open("assets/dict_of_local_r_to_k_rules_AND.pkl","wb"))
    return dict_of_local_r_to_k_rules

def create_dict_of_global_r_to_k_rules(dict_of_local_r_to_k_rules):
    dict_of_global_r_to_k_rules = dict()
    for mid in dict_of_local_r_to_k_rules:
        for rid,ruleset in dict_of_local_r_to_k_rules[mid].items():
            if rid not in dict_of_global_r_to_k_rules:
                dict_of_global_r_to_k_rules[rid] = ruleset
            else:
                dict_of_global_r_to_k_rules[rid].update(ruleset)

    pickle.dump(dict_of_global_r_to_k_rules, open("assets/dict_of_global_r_to_k_rules_AND.pkl","wb"))
    return dict_of_global_r_to_k_rules

##########################################
## Main
##########################################
def main():
    db = "module"
    entry_dict_path= "assets/module_entry_dict-2021_03_22.json"
    source_target_link_path = "assets/ko_rn_link_dicts-2021_03_10.json"

    if not os.path.exists("assets"):
        os.makedirs("assets")

    ##########################################
    ## Run everything, assuming no existing files
    ##########################################
    ## Retrieve/write KEGG files if don't exist; otherwise load them
    ## I don't actually use the rn_ko_dict and ko_rn_dict
    # module_entry_dict, rn_ko_dict, ko_rn_dict = collect_KEGG_files(db, entry_dict_path, source_target_link_path)
    with open(entry_dict_path) as f:
        module_entry_dict = json.load(f)

    ## Choose to load or calculate
    # calculated_module_dict = pickle.load(open('assets/calculated_module_dict.pkl', 'rb'))
    calculated_module_dict = parse_and_format_modules(module_entry_dict)

    # ## Choose to load or calculate
    # # dict_of_local_r_to_k_rules = pickle.load(open("assets/dict_of_local_r_to_k_rules.pkl","rb"))
    dict_of_local_r_to_k_rules = create_dict_of_local_r_to_k_rules(calculated_module_dict, module_entry_dict)

    # ## Choose to load or calculate
    # # dict_of_global_r_to_k_rules = pickle.load(open("assets/dict_of_global_r_to_k_rules.pkl","rb"))
    dict_of_global_r_to_k_rules = create_dict_of_global_r_to_k_rules(dict_of_local_r_to_k_rules)


if __name__ == '__main__':
    main()