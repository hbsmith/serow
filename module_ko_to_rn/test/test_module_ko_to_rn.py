import os
import json
import unittest
import sys
sys.path.append("..")
from module_ko_to_rn import *

class TestParseOrthologs(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.db = "module"
        self.entry_dict_path = "assets/module_entry_dict-2021_03_22.json"
        self.source_target_link_path = "assets/ko_rn_link_dicts-2021_03_10.json"
        self.searchExpr = declareSearchExpr()

        ## Load module entry dict if it exists
        if os.path.isfile(self.entry_dict_path):
            with open(self.entry_dict_path, 'r') as f:
                self.module_entry_dict = json.load(f)
        else:
            self.module_name_dict = create_id_name_dict(self.db)
            self.module_entry_dict = retrieve_entry_info(self.module_name_dict,self.db)
            with open(self.entry_dict_path, 'w') as f:
                json.dump(self.module_entry_dict, f, indent=4)

        ## Load link dicts if they exist
        if os.path.isfile(self.source_target_link_path):
            with open(self.source_target_link_path, 'r') as f:
                source_target_links = json.load(f)
                self.rn_ko_dict = source_target_links[0]
                self.ko_rn_dict = source_target_links[1]

        else:
            self.rn_ko_dict, self.ko_rn_dict = create_link_dicts("ko","reaction")
            with open(self.source_target_link_path, 'w') as f:
                json.dump([self.rn_ko_dict, self.ko_rn_dict], f, indent=4, default=convert_sets_to_lists)


    def test_expected_ko_sets(self):

        #= some no sign modules: M00001, M00007
        #= some (+) only modules: M00019, 'M00417','M00934'
        #= some (-) only modules: M00014, 'M00574','M00917'
        #= some both modules: M00009, M00374, M00846

        ## Based on Togo entries
        expected_module_dict = {
            "M00001": {
                frozenset(["K00844"]),
                frozenset(["K12407"]),
                frozenset(["K00845"]),
                frozenset(["K00886"]),
                frozenset(["K08074"]),
                frozenset(["K00918"]),
                frozenset(["K01810"]),
                frozenset(["K06859"]),
                frozenset(["K13810"]),
                frozenset(["K15916"]),
                frozenset(["K00850"]),
                frozenset(["K16370"]),
                frozenset(["K21071"]),
                frozenset(["K00918"]),
                frozenset(["K01623"]),
                frozenset(["K01624"]),
                frozenset(["K11645"]),
                frozenset(["K16305"]),
                frozenset(["K16306"]),
                frozenset(["K01803"]),
                frozenset(["K00134"]),
                frozenset(["K00150"]),
                frozenset(["K00927"]),
                frozenset(["K11389"]),
                frozenset(["K01834"]),
                frozenset(["K15633"]),
                frozenset(["K15634"]),
                frozenset(["K15635"]),
                frozenset(["K01689"]),
                frozenset(["K00873"]),
                frozenset(["K12406"]),
            },
            "M00009": {
                frozenset(["K01647"]),
                frozenset(["K05942"]),
                frozenset(["K01681"]),
                frozenset(["K01682"]),
                frozenset(["K00031"]),
                frozenset(["K00030"]),
                frozenset(["K00164","K00658","K00382"]),
                frozenset(["K01616","K00382"]),
                frozenset(["K00174","K00175"]),
                frozenset(["K00174","K00175","K00177"]),
                frozenset(["K00174","K00175","K00177","K00176"]),
                frozenset(["K00174","K00175","K00176"]),
                frozenset(["K01902","K01903"]),
                frozenset(["K01899","K01900"]),
                frozenset(["K18118"]),
                frozenset(["K00234","K00235","K00236","K00237"]),
                frozenset(["K00239","K00240","K00241"]),
                frozenset(["K00239","K00240","K00241","K00242"]),
                frozenset(["K00239","K00240","K00241","K18859"]),
                frozenset(["K00239","K00240","K00241","K18860"]),
                frozenset(["K00244","K00245","K00246"]),
                frozenset(["K00244","K00245","K00246","K00247"]),
                frozenset(["K01676"]),
                frozenset(["K01679"]),
                frozenset(["K01677","K01678"]),
                frozenset(["K00026"]),
                frozenset(["K00025"]),
                frozenset(["K00024"]),
                frozenset(["K00116"])
            },
            "M00417": {
                frozenset(["K02297","K02298","K02299","K02300"])
            },
            "M00574": {
                frozenset(["K11023"]),
                frozenset(["K11024"]),
                frozenset(["K11025"]),
                frozenset(["K11026"]),
                frozenset(["K11027"]),
            }

        }

        ## Testing
        calculated_module_dict = dict()
        for mid, d in self.module_entry_dict.items():
            if mid in expected_module_dict:
                data = getTopLevelOp(d["definition"])
                data = replaceStrsWithValidExprs(data) ## Transforms strings to "ValidExprs" objects that can be typechecked
                data = combineValidExprs(data) ## Transforms all but the top level to ValidExprs objs
                data = combineValidExprs(data) ## Transforms top level (Returns a single ValidExprs obj)
                calculated_module_dict[mid] = {frozenset(i) for i in data.expressions} - {frozenset()}
                # calculated_module_dict[mid] -= {frozenset()}

        for mid in expected_module_dict:
            print("Testing: ",mid)
            self.assertEqual(calculated_module_dict[mid], expected_module_dict[mid])

        ## Running the below just ensures no errors for all module parsing
        # calculated_module_dict = dict()
        # for mid, d in self.module_entry_dict.items():
        #     print(mid)
        #     data = getTopLevelOp(d["definition"])
        #     data = replaceStrsWithValidExprs(data) ## Transforms strings to "ValidExprs" objects that can be typechecked
        #     data = combineValidExprs(data) ## Transforms all but the top level to ValidExprs objs
        #     data = combineValidExprs(data) ## Transforms top level (Returns a single ValidExprs obj)
        #     calculated_module_dict[mid] = {frozenset(i) for i in data.expressions} - {frozenset()}
