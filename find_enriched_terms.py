# -*- coding: utf-8 -*-
"""
Created on Sun Sep 18 13:56:25 2016

@author: ellisrj2
"""

       
def remove_nonsig_david(wdir):
    import os, os.path
    import pandas 
    import itertools
    from collections import Counter
    import csv
    
    path = os.listdir(wdir)
    go_bp = [] #biological pathways
    most_common_terms_ids = []
    enriched_terms = []
    
    for f,i in zip(path,range(len(path))):
        if f[-3:] == 'txt':
            df = pandas.read_table(f)
            df = df[df['Fold Enrichment'] > 2]
            terms = list(set(df['Term']))
            go_bp.append(terms) # append terms to see which ones appear in multiple studies
            df.to_excel(f + '_processed.xlsx', index=False)
        
    allbpterms = itertools.chain(*go_bp)
    allbpterms = list(allbpterms)
    bp_count = Counter(allbpterms)
    most_common_terms = bp_count.most_common(1000)
    most_common_terms_terms = [x[0] for x in most_common_terms] #top 200 records 
    
    for i in range(len(most_common_terms_terms)):
        tilde_index = most_common_terms_terms[i].find('~')
        most_common_terms_ids.append(most_common_terms_terms[i][:tilde_index])
        enriched_terms.append(most_common_terms_terms[i][tilde_index+1:])
        
    most_common_terms_counts = [x[1] for x in most_common_terms] #top 200 counts

    rows = zip(most_common_terms_ids, enriched_terms, most_common_terms_counts)
    
    rows = [list(t) for t in rows]    
    
    with open('most_common_terms.csv', 'wb') as thefile:
        writer = csv.writer(thefile)
        for row in rows:
            writer.writerow(row)
        
        
def remove_background_david(enriched, background):
    import pandas 
    
    df_enriched = pandas.read_excel(enriched)
    df_background = pandas.read_excel(background)
    
    enriched_terms = df_enriched['Term']
    background_terms = df_background['Term']
    
    not_background = [term for term in enriched_terms if term not in background_terms]
    
    thefile = open('enriched_terms.txt', 'w')    
    for item in not_background:
        thefile.write("%s\n" % item)
    