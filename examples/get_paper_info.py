#!/usr/bin/env python3
"""Get paper info from CrossRef API"""
import json, sys
import urllib.request

# Zhao et al. PRC 94 (2016) - ImQMD for U+U
resp = urllib.request.urlopen("https://api.crossref.org/works/10.1103/PhysRevC.94.024601")
data = json.load(resp)
m = data['message']
print("="*60)
print("ZHAO et al. PRC 94, 024601 (2016)")
print("="*60)
print("Title:", m.get('title',[''])[0][:200])
print("Authors:", ', '.join([a.get('family','') for a in m.get('author',[])]))
print("Pages:", m.get('page','N/A'))
print("Cited by:", m.get('is-referenced-by-count',0))

# Shen et al. PRC 66 (2002) - HIVAP
resp2 = urllib.request.urlopen("https://api.crossref.org/works/10.1103/PhysRevC.66.061602")
data2 = json.load(resp2)
m2 = data2['message']
print()
print("="*60)
print("SHEN et al. PRC 66, 061602(R) (2002)")
print("="*60)
print("Title:", m2.get('title',[''])[0][:200])
print("Authors:", ', '.join([a.get('family','') for a in m2.get('author',[])]))
print("Cited by:", m2.get('is-referenced-by-count',0))

# Get references from Wang et al. paper
resp3 = urllib.request.urlopen("https://api.crossref.org/works/10.1016/j.nimb.2019.02.013")
data3 = json.load(resp3)
m3 = data3['message']
refs = m3.get('reference',[])
print()
print("="*60)
print("References from Wang et al. related to ImQMD/HIVAP:")
print("="*60)
for r in refs:
    t = r.get('article-title','') or r.get('unstructured','')
    if any(k in t.lower() for k in ['imqmd','gemini','hivap','quantum molecular','de-excit','evaporation','transfer reaction','zhao','shen']):
        print(f"  - {t[:120]}")
        print(f"    DOI: {r.get('DOI','N/A')}")
