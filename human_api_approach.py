import os
import re
import csv
import json
import time
import requests
from concurrent.futures import ProcessPoolExecutor, as_completed
from urllib.error import HTTPError
from typing import Dict, List, Optional, Any

# Third-party imports
from Bio import Entrez
from dotenv import load_dotenv
from google import genai
from google.genai import types

# Load environment variables
load_dotenv()

# =================================================================
# 1. CONFIGURATION & CONSTANTS
# =================================================================

# Optimize thread usage for cluster environments (slurm/HPC)
for env in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"]:
    os.environ[env] = "1"

# API Credentials
GEMINI_API_KEY = os.getenv("GEMINI_API_KEY")
Entrez.api_key = os.getenv("ENTREZ_API_KEY")
Entrez.email = os.getenv("ENTREZ_EMAIL")
UMD_COOKIES = os.getenv("UMD_COOKIES")

# Http Headers for scraping (mimics a browser to avoid 403 Forbidden)
HEADERS = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) Chrome/120.0.0.0 Safari/537.36",
    "Cookie": UMD_COOKIES
}

# --- SEARCH PARAMETERS ---
FORCE_PMIDS = [
    "33514696", # Doench/Li (Nature Biotech)
    "36737458", # Wei Li @ UMD (Nature Comms, DeepCas12a)
    "30962622"  # Chen (Nature Methods)
]

SEARCH_QUERIES = [
    '("Li W"[Last Author] AND "Maryland"[Affiliation]) AND ("Cas12" OR "Cas13" OR "Deep learning")',
    '("Li W"[First Author] AND "Doench JG"[Last Author])',
    '("Doench JG"[Last Author]) AND ("Screen" OR "Library")',
    '("Sanjana NE"[Last Author]) AND ("Cas12" OR "Cas13")',
    '("Zhang F"[Last Author]) AND ("Cas12" OR "Cas13")',
    '("Kleinstiver BP"[Last Author]) AND ("Cas12" OR "Cas13")',
    '("Kim HH"[Last Author]) AND ("Cas12" OR "Cas13")'
]

ROYAL_JOURNALS = [
    "Nature", "Nature biotechnology", "Nature methods", "Nature communications", 
    "Nature medicine", "Nature genetics", "Nature structural & molecular biology",
    "Science", "Science advances", "Science translational medicine",
    "Cell", "Molecular cell", "Cell reports", "Cell systems", "Cell genomics",
    "Molecular therapy", "Nucleic acids research"
]

# =================================================================
# 2. HELPER FUNCTIONS (The "Modules")
# =================================================================

def fetch_metadata(pmid: str) -> Optional[Dict[str, Any]]:
    """
    Retrieves metadata (Journal, Title) for a given PMID from PubMed.

    Args:
        pmid: The PubMed ID to query.

    Returns:
        A dictionary with 'Journal' and 'Title' if successful, else None.
    """
    try:
        with Entrez.efetch(db="pubmed", id=pmid, retmode="xml") as handle:
            records = Entrez.read(handle)
        
        if not records.get('PubmedArticle'):
            return None
            
        article = records['PubmedArticle'][0]['MedlineCitation']['Article']
        return {
            'Journal': article['Journal']['Title'],
            'Title': article['ArticleTitle']
        }
    except Exception as e:
        print(f"[WARN] Metadata fetch failed for {pmid}: {e}")
        return None

def fetch_full_text(pmid: str) -> str:
    """
    Attempts to retrieve full text via PMC or DOI scraping.

    Args:
        pmid: The PubMed ID.

    Returns:
        The extracted text string (up to 40k chars) or "ABSTRACT_ONLY".
    """
    try:
        # Strategy A: Check for free PMC full text
        link_handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
        link_record = Entrez.read(link_handle)
        
        if link_record[0]['LinkSetDb']:
            pmc_id = link_record[0]['LinkSetDb'][0]['Link'][0]['Id']
            with Entrez.efetch(db="pmc", id=pmc_id, retmode="xml") as handle:
                return f"SOURCE: PMC_FULL\n{handle.read().decode('utf-8')[:40000]}"
        
        # Strategy B: Fallback to DOI scraping (using cookies)
        handle = Entrez.esummary(db="pubmed", id=pmid)
        doi = Entrez.read(handle)[0].get('ELocationID', '').replace('doi: ', '')
        
        if doi:
            resp = requests.get(f"https://doi.org/{doi}", headers=HEADERS, timeout=10)
            if resp.status_code == 200:
                clean_text = re.sub('<[^<]+?>', ' ', resp.text)
                return f"SOURCE: COOKIE_SCRAPE\n{clean_text[:35000]}"
                
    except Exception:
        pass # Fail silently and return abstract only
        
    return "ABSTRACT_ONLY"

def analyze_with_ai(text: str, pmid: str) -> Dict[str, Any]:
    """
    Uses Gemini 2.0 Flash to extract CRISPR screen metadata from text.

    Args:
        text: The scientific text to analyze.
        pmid: The PMID (for reference).

    Returns:
        A dictionary containing the AI's analysis (Lab, Method, N, Keep boolean).
    """
    client = genai.Client(api_key=GEMINI_API_KEY)
    
    prompt = f"""
    Analyze this paper for a CRISPR SCREEN/DATASET.
    
    CRITICAL OBJECTIVES:
    1. Identify the Senior Author/Lab Head.
    2. If Author is Wei Li, check affiliation (UMD vs Broad).
    3. Extract method (e.g. "Genome-wide Cas12 Screen").
    4. Estimate 'N' (number of samples/guides).

    INPUT TEXT: {text[:30000]}
    
    RETURN JSON format: {{"keep": bool, "lab": "str", "method": "str", "n": "str"}}
    """
    
    try:
        response = client.models.generate_content(
            model="gemini-2.0-flash", # Upgraded from 1.5-flash-lite for better reasoning
            contents=prompt,
            config=types.GenerateContentConfig(response_mime_type="application/json")
        )
        return json.loads(response.text)
    except Exception as e:
        print(f"[ERROR] AI Analysis failed for {pmid}: {e}")
        return {"keep": False}

# =================================================================
# 3. MAIN PIPELINE FUNCTION
# =================================================================

def process_paper(pmid: str) -> Optional[Dict[str, str]]:
    """
    Orchestrates the entire pipeline for a single paper.
    
    Args:
        pmid: The PubMed ID to process.
        
    Returns:
        A clean dictionary row for the CSV if successful, else None.
    """
    # 1. Rate Limit Safety
    time.sleep(0.5)
    
    # 2. Get Metadata
    metadata = fetch_metadata(pmid)
    if not metadata: 
        return None
    
    journal = metadata['Journal']
    
    # 3. Apply "Journal Police" Filter
    # (Unless it's in our FORCE list, we skip low-tier journals)
    is_royal = any(j.lower() in journal.lower() for j in ROYAL_JOURNALS)
    if not is_royal and pmid not in FORCE_PMIDS:
        return None

    # 4. Get Content
    text_content = fetch_full_text(pmid)
    
    # 5. Run AI Analysis
    ai_result = analyze_with_ai(text_content, pmid)
    
    # 6. Final Decision
    if pmid in FORCE_PMIDS or ai_result.get("keep"):
        return {
            "PMID": pmid,
            "Journal": journal,
            "Lab": ai_result.get("lab", "Unknown"),
            "Method": ai_result.get("method", "CRISPR Screen"),
            "N": ai_result.get("n", "Unknown"),
            "URL": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        }
    
    return None

# =================================================================
# 4. EXECUTION ENTRY POINT
# =================================================================
if __name__ == "__main__":
    print("--- ⚔️ STARTING WEI LI (UMD) & ELITE PI HUNT ---")
    
    all_pmids = []
    
    # A. Add Forced PMIDs
    all_pmids.extend(FORCE_PMIDS)

    # B. Run Search Queries
    print(f"--- 🔍 Querying PubMed for {len(SEARCH_QUERIES)} patterns ---")
    for query in SEARCH_QUERIES:
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=20)
            ids = Entrez.read(handle)["IdList"]
            all_pmids.extend(ids)
        except Exception as e:
            print(f"   > Search failed for query '{query[:20]}...': {e}")

    # C. Deduplicate
    unique_pmids = list(set(all_pmids))
    print(f"--- 🚀 Analyzing {len(unique_pmids)} Unique Papers ---")

    # D. Parallel Processing
    results = []
    with ProcessPoolExecutor(max_workers=5) as executor: # Reduced to 5 to avoid API limits
        futures = {executor.submit(process_paper, p): p for p in unique_pmids}
        
        for future in as_completed(futures):
            try:
                res = future.result()
                if res:
                    print(f"     ✅ KEPT: {res['Lab']} | {res['Journal']}")
                    results.append(res)
            except Exception as e:
                print(f"     ❌ CRITICAL ERROR in worker: {e}")

    # E. Save Results
    output_file = "big_boy_wei_li_final.csv"
    if results:
        with open(output_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["PMID", "Journal", "Lab", "Method", "N", "URL"])
            writer.writeheader()
            writer.writerows(results)
        print(f"--- DONE. Found {len(results)} Papers. Saved to '{output_file}'. ---")
    else:
        print("--- DONE. No relevant papers found. ---")