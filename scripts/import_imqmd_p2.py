#!/usr/bin/env python3
"""Import remaining ImQMD papers into Zotero."""
import json, os, sys, time, re, urllib.request, xml.etree.ElementTree as ET

ZOTERO_API_KEY = "QRBdrYsc9mlUwAEJknuzehPn"
ZOTERO_USER_ID = 9566388
BASE = f"https://api.zotero.org/users/{ZOTERO_USER_ID}"
HEADERS = {"Zotero-API-Key": ZOTERO_API_KEY}
PDF_DIR = "/home/zhuhaofan/work/agent work/mnt-sim/papers"
NS = {"a": "http://www.w3.org/2005/Atom"}
os.makedirs(PDF_DIR, exist_ok=True)

def _opener():
    proxy = "http://127.0.0.1:7890"
    handler = urllib.request.ProxyHandler({"https": proxy, "http": proxy})
    return urllib.request.build_opener(handler)
opener = _opener()

# Papers to add (from INSPIRE-HEP)
PAPERS_WITH_ARXIV = [
    ("nucl-th/0402066", "Wang et al. 2004 ImQMD-II"),
    ("1205.1605", "Zhang-Li-Zhou-Tsang 2012 cluster recognition"),
    ("1402.3790", "Zhang-Tsang-Li-Liu 2014 effective mass splitting"),
    ("1701.02082", "Yao-Wang 2017 Kr+Ni MNT"),
]

# Paper without arXiv — use DOI metadata
PAPERS_DOI = [
    ("10.1016/j.physletb.2017.11.060", "Li et al. 2018 Xe+Pt MNT"),
]

COLLECTION_KEY = "VZQRQTT2"  # ImQMD Model

def fetch_arxiv_meta(arxiv_id):
    url = f"https://export.arxiv.org/api/query?id_list={arxiv_id}"
    req = urllib.request.Request(url, headers={"User-Agent": "ZoteroImport/1.0"})
    with urllib.request.urlopen(req, timeout=15) as r:
        data = r.read()
    root = ET.fromstring(data)
    entry = root.find("a:entry", NS)
    if entry is None: return None
    title = entry.find("a:title", NS).text.strip().replace("\n", " ")
    abstract = entry.find("a:summary", NS).text.strip()[:5000]
    published = entry.find("a:published", NS).text[:10]
    authors = []
    for a in entry.findall("a:author", NS):
        name = a.find("a:name", NS).text.strip().split()
        if len(name) > 1:
            authors.append({"creatorType":"author","firstName":" ".join(name[:-1]),"lastName":name[-1]})
        else:
            authors.append({"creatorType":"author","lastName":name[0]})
    doi = ""
    for link in entry.findall("a:link", NS):
        if link.get("title") == "doi": doi = link.get("href", ""); break
    return {"title": title, "abstractNote": abstract, "date": published,
            "creators": authors, "DOI": doi,
            "url": f"https://arxiv.org/abs/{arxiv_id}",
            "extra": f"arXiv: {arxiv_id}\n"}

def fetch_doi_meta(doi):
    """Fetch metadata from CrossRef for papers without arXiv."""
    url = f"https://api.crossref.org/works/{doi}"
    req = urllib.request.Request(url, headers={"User-Agent": "ZoteroImport/1.0"})
    with urllib.request.urlopen(req, timeout=15) as r:
        data = json.loads(r.read())
    msg = data["message"]
    title = msg.get("title", ["Unknown"])[0]
    abstract = msg.get("abstract", "")[:5000] if msg.get("abstract") else ""
    date_parts = msg.get("published-print", {}).get("date-parts", [[2000]])[0]
    year = str(date_parts[0])
    authors = []
    for a in msg.get("author", []):
        given = a.get("given", "")
        family = a.get("family", "")
        if given and family:
            authors.append({"creatorType": "author", "firstName": given, "lastName": family})
        elif family:
            authors.append({"creatorType": "author", "lastName": family})
    return {"title": title, "abstractNote": abstract, "date": year,
            "creators": authors, "DOI": doi,
            "url": f"https://doi.org/{doi}",
            "extra": f"DOI: {doi}\n"}

def download_pdf(arxiv_id):
    safe = arxiv_id.replace("/", "_")
    path = os.path.join(PDF_DIR, f"{safe}.pdf")
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        return path
    url = f"https://arxiv.org/pdf/{arxiv_id}"
    req = urllib.request.Request(url, headers={"User-Agent": "ZoteroImport/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            data = r.read()
        if data[:4] == b"%PDF":
            with open(path, "wb") as f: f.write(data)
            return path
    except Exception as e:
        print(f"  ⚠ PDF download failed: {e}")
    return None

def create_item(meta, collection_key):
    item = {"itemType": "journalArticle", "title": meta["title"],
            "creators": meta["creators"], "abstractNote": meta["abstractNote"],
            "date": meta["date"], "DOI": meta["DOI"], "url": meta["url"],
            "extra": meta["extra"], "collections": [collection_key],
            "tags": [{"tag": "ImQMD"}, {"tag": "文献调研导入"}],
            "relations": {}}
    body = json.dumps([item], ensure_ascii=False).encode("utf-8")
    req = urllib.request.Request(f"{BASE}/items", data=body,
        headers={**HEADERS, "Content-Type": "application/json"}, method="POST")
    with opener.open(req, timeout=15) as r:
        result = json.loads(r.read())
        successful = result.get("successful", {})
        if successful:
            return successful[list(successful.keys())[0]]["data"]["key"]
    return None

def attach_pdf(item_key, pdf_path):
    att = [{"itemType": "attachment", "linkMode": "linked_file",
            "path": os.path.abspath(pdf_path), "contentType": "application/pdf",
            "parentItem": item_key, "title": os.path.basename(pdf_path),
            "relations": {}}]
    body = json.dumps(att, ensure_ascii=False).encode("utf-8")
    req = urllib.request.Request(f"{BASE}/items", data=body,
        headers={**HEADERS, "Content-Type": "application/json"}, method="POST")
    try:
        with opener.open(req, timeout=30) as r:
            return bool(json.loads(r.read()).get("successful"))
    except: return False

def main():
    print("=" * 60)
    print("Import remaining ImQMD papers")
    print("=" * 60)

    # arXiv papers
    for arxiv_id, label in PAPERS_WITH_ARXIV:
        print(f"\n📄 {label} ({arxiv_id})")
        time.sleep(2)
        meta = fetch_arxiv_meta(arxiv_id)
        if not meta:
            print("  ❌ Failed to fetch metadata")
            continue
        print(f"  Title: {meta['title'][:80]}...")

        pdf = download_pdf(arxiv_id)
        print(f"  📥 PDF: {'OK' if pdf else 'N/A'}")

        item_key = create_item(meta, COLLECTION_KEY)
        if not item_key:
            print("  ❌ Failed to create item")
            continue
        print(f"  ✅ Created: {item_key}")

        if pdf:
            time.sleep(1)
            attach_pdf(item_key, pdf)
            print("  📎 PDF attached")
        time.sleep(2)

    # DOI-only papers
    for doi, label in PAPERS_DOI:
        print(f"\n📄 {label} (DOI: {doi})")
        time.sleep(2)
        meta = fetch_doi_meta(doi)
        if not meta:
            print("  ❌ Failed")
            continue
        print(f"  Title: {meta['title'][:80]}...")
        print("  ⚠ No arXiv PDF available (DOI only)")

        item_key = create_item(meta, COLLECTION_KEY)
        if item_key:
            print(f"  ✅ Created: {item_key}")
        else:
            print("  ❌ Failed")
        time.sleep(2)

    print(f"\n{'='*60}")
    print("Done.")

if __name__ == "__main__":
    main()
