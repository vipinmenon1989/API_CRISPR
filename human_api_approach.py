import os
import re
import csv
import json
import requests
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import Entrez
from urllib.error import HTTPError
from google import genai
from google.genai import types
from dotenv import load_dotenv
load_dotenv()
# =================================================================
# 1. CONFIGURATION
# =================================================================
for env in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"]:
    os.environ[env] = "1"

GEMINI_API_KEY = os.getenv("GEMINI_API_KEY")
Entrez.api_key = os.getenv("ENTREZ_API_KEY")
Entrez.email = os.getenv("ENTREZ_EMAIL")
UMD_COOKIES = os.getenv("UMD_COOKIES")

HEADERS = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) Chrome/120.0.0.0 Safari/537.36",
    "Cookie": UMD_COOKIES
}

# --- 1. THE "MUST HAVE" PAPERS (Hard-Coded) ---
# We force these into the list first.
FORCE_PMIDS = [
    "33514696", # Doench/Li (Nature Biotech)
    "36737458", # Wei Li @ UMD (Nature Comms, DeepCas12a)
    "30962622"  # Chen (Nature Methods)
]

# --- 2. THE WEI LI & PI SEARCH QUERIES ---
SEARCH_QUERIES = [
    # TARGET 1: Wei Li at UMD (Last Author)
    '("Li W"[Last Author] AND "Maryland"[Affiliation]) AND ("Cas12" OR "Cas13" OR "Deep learning")',
    
    # TARGET 2: Wei Li with Doench (First Author - Historical)
    '("Li W"[First Author] AND "Doench JG"[Last Author])',
    
    # TARGET 3: The Other Giants (Strict Last Author)
    '("Doench JG"[Last Author]) AND ("Screen" OR "Library")',
    '("Sanjana NE"[Last Author]) AND ("Cas12" OR "Cas13")',
    '("Zhang F"[Last Author]) AND ("Cas12" OR "Cas13")',
    '("Kleinstiver BP"[Last Author]) AND ("Cas12" OR "Cas13")',
    '("Kim HH"[Last Author]) AND ("Cas12" OR "Cas13")'
]

# --- 3. THE ROYAL JOURNALS (Strict Filter) ---
ROYAL_JOURNALS = [
    "Nature", "Nature biotechnology", "Nature methods", "Nature communications", 
    "Nature medicine", "Nature genetics", "Nature structural & molecular biology",
    "Science", "Science advances", "Science translational medicine",
    "Cell", "Molecular cell", "Cell reports", "Cell systems", "Cell genomics",
    "Molecular therapy", "Nucleic acids research"
]

def analyze_paper(pmid):
    try:
        time.sleep(0.5) 
        
        # 1. Fetch Metadata
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        if not records['PubmedArticle']: return None
        
        article = records['PubmedArticle'][0]['MedlineCitation']['Article']
        journal = article['Journal']['Title']
        title = article['ArticleTitle']
        
        # --- JOURNAL POLICE ---
        is_royal = any(j.lower() in journal.lower() for j in ROYAL_JOURNALS)
        if not is_royal and pmid not in FORCE_PMIDS:
            return None # Skip garbage (unless it's forced)

        # 2. Fetch Text
        text_content = "ABSTRACT_ONLY"
        
        link_handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
        link_record = Entrez.read(link_handle)
        if link_record[0]['LinkSetDb']:
            pmc_id = link_record[0]['LinkSetDb'][0]['Link'][0]['Id']
            fetch_handle = Entrez.efetch(db="pmc", id=pmc_id, retmode="xml")
            text_content = f"SOURCE: PMC_FULL\n{fetch_handle.read().decode('utf-8')[:40000]}"
        else:
            handle = Entrez.esummary(db="pubmed", id=pmid)
            doi = Entrez.read(handle)[0].get('ELocationID', '').replace('doi: ', '')
            if doi:
                resp = requests.get(f"https://doi.org/{doi}", headers=HEADERS, timeout=8)
                if resp.status_code == 200:
                    clean_text = re.sub('<[^<]+?>', ' ', resp.text)
                    text_content = f"SOURCE: COOKIE_SCRAPE\n{clean_text[:35000]}"

        # 3. AI Analysis
        client = genai.Client(api_key=GEMINI_API_KEY)
        prompt = f"""
        Analyze this paper for a CRISPR SCREEN/DATASET.
        
        CRITICAL: Identify the Senior Author/Lab Head.
        If the author is Wei Li, check if affiliation is University of Maryland or Broad Institute.
        
        INPUT: {text_content[:30000]}
        
        Return JSON: {{"keep": bool, "lab": "str", "method": "str", "n": "str"}}
        """
        
        response = client.models.generate_content(
            model="gemini-2.5-flash-lite",
            contents=prompt,
            config=types.GenerateContentConfig(response_mime_type="application/json")
        )
        data = json.loads(response.text)
        
        # Always keep FORCE_PMIDS, otherwise check AI filter
        if pmid in FORCE_PMIDS or data.get("keep"):
            return {
                "PMID": pmid,
                "Journal": journal,
                "Lab": data.get("lab", "Wei Li (UMD)" if pmid in FORCE_PMIDS else "Unknown"),
                "Method": data.get("method", "CRISPR Screen"),
                "N": data.get("n", "Unknown"),
                "URL": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            }
    except: return None

if __name__ == "__main__":
    all_pmids = []
    print("--- ⚔️ STARTING WEI LI (UMD) & ELITE PI HUNT ---")
    
    # 1. ADD FORCE LIST
    all_pmids.extend(FORCE_PMIDS)

    # 2. RUN SEARCH QUERIES
    for query in SEARCH_QUERIES:
        print(f"   > Searching: {query[:50]}...")
        handle = Entrez.esearch(db="pubmed", term=query, retmax=20)
        ids = Entrez.read(handle)["IdList"]
        all_pmids.extend(ids)

    unique_pmids = list(set(all_pmids))
    print(f"--- 🚀 Analyzing {len(unique_pmids)} Papers ---")

    results = []
    with ProcessPoolExecutor(max_workers=10) as executor:
        futures = {executor.submit(analyze_paper, p): p for p in unique_pmids}
        for future in as_completed(futures):
            res = future.result()
            if res:
                print(f"     ✅ KEPT: {res['Lab']} | {res['Journal']}")
                results.append(res)

    with open("big_boy_wei_li_final.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["PMID", "Journal", "Lab", "Method", "N", "URL"])
        writer.writeheader()
        writer.writerows(results)
    
    print(f"--- DONE. Found {len(results)} Papers. Check 'big_boy_wei_li_final.csv'. ---")
