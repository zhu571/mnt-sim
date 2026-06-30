#!/usr/bin/env python3
"""Import ImQMD literature papers into Zotero with PDFs."""
import json, os, sys, time, re, urllib.request, xml.etree.ElementTree as ET

ZOTERO_API_KEY = "QRBdrYsc9mlUwAEJknuzehPn"
ZOTERO_USER_ID = 9566388
BASE = f"https://api.zotero.org/users/{ZOTERO_USER_ID}"
HEADERS = {"Zotero-API-Key": ZOTERO_API_KEY}
PDF_DIR = "/home/zhuhaofan/work/agent work/mnt-sim/papers"
COLLECTION_NAME = "ImQMD Model"
NS = {"a": "http://www.w3.org/2005/Atom"}
os.makedirs(PDF_DIR, exist_ok=True)

# Build proxy opener for Zotero API
def _opener():
    proxy = "http://127.0.0.1:7890"
    handler = urllib.request.ProxyHandler({"https": proxy, "http": proxy})
    return urllib.request.build_opener(handler)
opener = _opener()

# ── Paper list: (arxiv_id, display_label) ──
PAPERS = [
    ("nucl-th/0201079",  "Wang-Li-Wu 2002 ImQMD original"),
    ("2005.12877",       "Zhang et al. 2020 QMD review"),
    ("1405.5271",        "Wang-Ou-Zhang-Li 2014 initialization"),
    ("1009.1928",        "Zhang et al. 2012 NN cross-sections"),
    ("2103.13218",       "Chen-Zhang-Li 2021 Pauli blocking"),
    ("2403.00343",       "Chen et al. 2024 novel Pauli blocking"),
    ("1605.07393",       "Zhao et al. 2016 U+U MNT"),
    ("0811.3107",        "Tsang et al. 2009 symmetry energy"),
    ("1309.7534",        "Wei-Wang-Ou 2014 light clusters"),
]

# ── Helpers ──
def build_pdf_path(arxiv_id):
    safe = arxiv_id.replace("/", "_")
    return os.path.join(PDF_DIR, f"{safe}.pdf")

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

def download_pdf(arxiv_id):
    """Download PDF from arXiv; return path or None."""
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
        print(f"  ⚠ PDF download failed for {arxiv_id}: {e}")
    return None

def get_or_create_collection(name):
    req = urllib.request.Request(f"{BASE}/collections?limit=100", headers=HEADERS)
    with opener.open(req, timeout=15) as r:
        for col in json.loads(r.read()):
            if col["data"]["name"] == name:
                return col["data"]["key"]
    body_data = {"name": name}
    body = json.dumps([body_data], ensure_ascii=False).encode("utf-8")
    req = urllib.request.Request(f"{BASE}/collections", data=body, headers=HEADERS, method="POST")
    with opener.open(req, timeout=15) as r:
        result = json.loads(r.read())
        successful = result.get("successful", {})
        if successful:
            return successful[list(successful.keys())[0]]["data"]["key"]
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
    """Linked file attachment."""
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

def check_existing(arxiv_id):
    """Check if paper already in Zotero by arXiv ID in extra field."""
    url = f"{BASE}/items/top?limit=100&itemType=journalArticle"
    while url:
        req = urllib.request.Request(url, headers=HEADERS)
        with opener.open(req, timeout=15) as r:
            data = json.loads(r.read())
        for item in data:
            extra = item["data"].get("extra", "") or ""
            if arxiv_id in extra:
                return True
        url = None
        # Check Link header for pagination
    return False

# ── Main ──
def main():
    print("=" * 60)
    print(f"Importing {len(PAPERS)} ImQMD papers into Zotero")
    print(f"Collection: {COLLECTION_NAME} | PDF dir: {PDF_DIR}")
    print("=" * 60)

    collection_key = get_or_create_collection(COLLECTION_NAME)
    print(f"\n📂 Collection key: {collection_key}")
    time.sleep(1)

    for arxiv_id, label in PAPERS:
        print(f"\n📄 {label} ({arxiv_id})")
        if check_existing(arxiv_id):
            print("  ⏭  Already in Zotero, skipping")
            continue

        time.sleep(2)
        meta = fetch_arxiv_meta(arxiv_id)
        if not meta:
            print("  ❌ Failed to fetch arXiv metadata")
            continue
        print(f"  Title: {meta['title'][:80]}...")

        pdf = download_pdf(arxiv_id)
        if pdf:
            print(f"  📥 PDF: {os.path.basename(pdf)}")
        else:
            print("  ⚠ No PDF")

        item_key = create_item(meta, collection_key)
        if not item_key:
            print("  ❌ Failed to create Zotero item")
            continue
        print(f"  ✅ Created: {item_key}")

        if pdf:
            time.sleep(1)
            if attach_pdf(item_key, pdf):
                print(f"  📎 PDF attached")
            else:
                print(f"  ⚠ PDF attach failed")

        time.sleep(2)

    print(f"\n{'='*60}")
    print("Done.")

if __name__ == "__main__":
    main()
