#!/usr/bin/env python
'''
    OMAmer to best pairs

    Need to deal with -top_n_fams == 10

    -- Alex Warwick Vesztrocy, November 2023
'''
from collections import defaultdict
from tqdm.auto import tqdm
import itertools
import pandas as pd
import sys


omamer_res_fns = sys.argv[1:]

# Load the hog mapping from the omamer results files
hogs = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

for (sp_i, fn) in tqdm(enumerate(omamer_res_fns), desc='Processing OMAmer result files', total=len(omamer_res_fns)):
    df = pd.read_csv(fn, sep='\t', comment='!')
    # filter out unplaced
    df = df[~df.hogid.isna()]
    
    # keep best scoring hit per-query
    df = df.drop_duplicates('qseqid', keep='first')

    # sort by scores globally
    df = df.sort_values(['subfamily_score', 'qseq_overlap', 'family_normcount', 'family_p'])

    # keep best scoring entry per hog
    df = df.drop_duplicates('hogid', keep='first')
    
    for (_, r) in df.iterrows():
        fam = int(r.hogid.split('.')[0][5:])
        hogs[fam][r.hogid][sp_i].add(r.qseqid.split('|')[1])

for fam in tqdm(hogs, desc='Generating pairs', unit='family'):
    for (hog, hog_members) in hogs[fam].items():
        # first save precise placed subhog.
        species = list(hog_members.items())
        for i in range(len(species)):
            for j in range(i+1, len(species)):
                for gi in species[i][1]:
                    for gj in species[j][1]:
                        print(gi, gj, sep='\t')

        hog_id = hog.split('.')
        if len(hog_id) > 1:
            for i in range(1, len(hog_id)):
                # add pairs with members higher up in the hog hierarchy
                parent_hog_id = '.'.join(hog_id[:i])
                parent_species = hogs[fam].get(parent_hog_id, None)
                if parent_species is not None:
                    parent_species = list(parent_species.items())
                    for i in range(len(species)):
                        for j in range(len(parent_species)):
                            if species[i][0] != parent_species[j][0]:
                                # filtering out the same species in parental. (ideally these shouldn't happen.
                                for gi in species[i][1]:
                                    for gj in parent_species[j][1]:
                                        print(gi, gj, sep='\t')
