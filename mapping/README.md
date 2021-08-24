There are a few small differences between the mapping produced when running `mapping_v5.py` with the stored links/entries, and the version produced at the end of "mapping week". This is due to the links being retrieved live at the time, without saving links that could be version controlled.

Here are the differences (`True` indicates no differences). The new mapping name is given first, the old mapping name second:
  ```
  map_rn2ko_viaMO_noaddition module_exlusive_to_non_addition
  True
  map_rn2ko_viaMO_addition1minus module_addition_1
  True
  new - old:  set()
  old - new:  {'R02289'}
  map_rn2ko_viaMO_addition1minus_extras module_addition_1_extras
  True
  map_rn2ko_viaMO_addition2plus module_addition_2plus_all
  True
  map_rn2ko_viaKO_1minus non-module_singletons
  False
  new - old:  {'R09797', 'R12760', 'R12757', 'R12787', 'R12790', 'R12769', 'R12770', 'R12768', 'R12775'}
  old - new:  {'R11817', 'R11816'}
  map_rn2ko_viaKO_2plus non-module_manual
  True
  map_rn2ko_viaEC_1minus rn2ko_via_ec_resolved
  True
  new - old:  set()
  old - new:  {nan}
  map_rn2ko_viaEC_2plus rn2ko_via_ec_ambig
  True
  new - old:  set()
  old - new:  {nan}
  map_rn2ko_spontaneous spontaneous_rns
  True
  ```
First, `R02289` is no longer present in the `map_rn2ko_viaMO_addition1minus` dataset. At the end of mapping week, KEGG showed it as being associated with only a single KO, `{frozenset({'K15023'})}`. Today, it is associated with both `K25221` and `K15023`. So it now falls into the category of reactions which need to be manually mapped. 

Second, reactions `R11816` and `R11817` are no longer in `map_rn2ko_viaKO_1minus`. They used to be associated with a single KO each:
  ```
  non-module_singletons {frozenset({'K07215'})}
  non-module_singletons {frozenset({'K07215'})}
  ```
Now KEGG shows them as being associated with both `K07226` and `K07215`. So they now fall into the category of reactions which need to be manually mapped. 

Third, there are now a number of new reactions in the `map_rn2ko_viaKO_1minus` category: `{'R09797', 'R12760', 'R12757', 'R12787', 'R12790', 'R12769', 'R12770', 'R12768', 'R12775'}`. 
- 'R09797' used to only be associated to KO via EC, but now is directly associated to `K17948`.
- The rest of the reactions appear to all be new reactions recently added to KEGG

Lastly, some nans (which were never used) were removed from the new mapping.
