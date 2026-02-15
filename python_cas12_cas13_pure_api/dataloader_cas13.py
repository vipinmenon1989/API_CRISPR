import os

# ==========================================
# 0. CRITICAL: HPC ENVIRONMENT SETUP
# ==========================================
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1" 

import json
import time
import random
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import Entrez
from google import genai
from google.genai import types
from playwright.sync_api import sync_playwright

# ==========================================
# CONFIGURATION
# ==========================================
GEMINI_API_KEY = "AIzaSyCtfb3_K9uwCrAxGQdKl-DOi_sCmkSlwms"
NCBI_API_KEY = "350c406985265581bb3a3b9c8ba345a88608"
INSTITUTIONAL_EMAIL = "vmenon@som.umaryland.edu"

Entrez.email = INSTITUTIONAL_EMAIL
Entrez.api_key = NCBI_API_KEY

# SCALE: 50 Workers is the sweet spot for your 100 CPUs 
# to leave room for the browser overhead.
MAX_WORKERS = 50 
PAPER_DIR = "paper_Cas13"

# ==========================================
# WORKER LOGIC
# ==========================================

def get_metadata_safe(pmid):
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        article = records['PubmedArticle'][0]['MedlineCitation']['Article']
        journal = article['Journal']['Title'].split(':')[0].replace(" ", "_")
        author = article['AuthorList'][0]['LastName']
        year = article['Journal']['JournalIssue']['PubDate'].get('Year', '2026')
        
        folder_name = f"{journal}_{author}_{year}"
        full_path = os.path.join(PAPER_DIR, folder_name)
        
        doi_list = [str(i) for i in article.get('ELocationID', []) if i.attributes.get('EIdType') == 'doi']
        doi = doi_list[0] if doi_list else f"PMID_{pmid}"
        return full_path, article['ArticleTitle'], doi, year
    except: return None, None, None, None

def analyze_ai_high_throughput(title, doi):
    """Aggressive filtering for High-Throughput RNA data."""
    client = genai.Client(api_key=GEMINI_API_KEY) 
    prompt = f"""
    Analyze if this study has ORIGINAL human Cas13 (RNA-targeting) efficiency data (N > 1000).
    Title: {title} | DOI: {doi}
    
    SEARCH FOR:
    - Massively Parallel Reporter Assays (MPRA)
    - Tiling or Pooled RNA-targeting screens
    - DeepCas13/CasRx optimization training data
    - Sanjana Lab / Wei Lei Lab datasets
    
    Return JSON ONLY: {{"is_relevant": bool, "n_size": int, "download_url": "str"}}
    """
    try:
        # Use a faster model to crush the 10,000 paper queue
        response = client.models.generate_content(
            model="gemini-2.5-flash-lite",
            contents=prompt,
            config=types.GenerateContentConfig(response_mime_type="application/json")
        )
        data = json.loads(response.text)
        # Ensure n_size is always an int
        data["n_size"] = int(data.get("n_size", 0) or 0)
        return data
    except: return {"is_relevant": False, "n_size": 0}

def download_stealth_hpc(url, folder):
    """Stealth browser for HPC environments."""
    if not url or "http" not in url: return False
    os.makedirs(folder, exist_ok=True)
    
    # Jitter to bypass rate limits
    time.sleep(random.uniform(0.5, 3.0))

    try:
        with sync_playwright() as p:
            browser = p.chromium.launch(headless=True, args=["--no-sandbox", "--disable-gpu"])
            # Mimic a high-resolution researcher workstation
            context = browser.new_context(
                user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/121.0.0.0 Safari/537.36",
                viewport={'width': 1920, 'height': 1080}
            )
            page = context.new_page()
            
            # Navigate with generous timeout for complex journal pages
            page.goto(url, wait_until="domcontentloaded", timeout=60000)

            # Hunter for Supplemental Buttons
            selectors = [
                "a:has-text('Supplementary Table')",
                "a:has-text('Supplementary Data')",
                "a:has-text('Download PDF')",
                "a[href*='.xlsx']",
                "a[href*='.xls']"
            ]
            
            for sel in selectors:
                if page.locator(sel).first.is_visible():
                    with page.expect_download(timeout=30000) as download_info:
                        page.locator(sel).first.click()
                    download = download_info.value
                    download.save_as(os.path.join(folder, download.suggested_filename))
                    browser.close()
                    return True
            
            # If no download button, take a snapshot of the abstract/data table
            page.screenshot(path=os.path.join(folder, "data_snapshot.png"))
            browser.close()
            return True
    except: return False

def process_one_paper(pmid):
    folder, title, doi, year = get_metadata_safe(pmid)
    if not title: return None
    
    # 1. Skip check
    if os.path.exists(folder) and len(os.listdir(folder)) > 0:
        return f"⏩ Skipped: {year} | {pmid}"

    # 2. AI Filter
    analysis = analyze_ai_high_throughput(title, doi)
    
    # 3. Download if it passes the N-size threshold
    if analysis.get("is_relevant") and analysis["n_size"] >= 1000:
        url = analysis.get('download_url') or f"https://doi.org/{doi}"
        success = download_stealth_hpc(url, folder)
        if success:
            return f"✅ SUCCESS [{year}] (N={analysis['n_size']}): {folder}"
        else:
            if os.path.exists(folder) and not os.listdir(folder): os.rmdir(folder)
            return f"❌ FAILED [{year}]: {title[:30]}"
    return None

# ==========================================
# MAIN EXECUTION
# ==========================================

if __name__ == "__main__":
    os.makedirs(PAPER_DIR, exist_ok=True)
    all_ids = []
    
    print(f"--- 🔎 Deep Scanning PubMed for 10,000 papers (2014-2026) ---")
    for year in range(2014, 2027):
        query = f'("Cas13" OR "Cas13a" OR "CasRx" OR "Cas13d" OR "C2c2") AND {year}[dp]'
        # Increased retmax for each year to get the maximum hits
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1000)
        record = Entrez.read(handle)
        all_ids.extend(record["IdList"])
        handle.close()
    
    pmid_list = list(set(all_ids))
    print(f"--- 🚀 TOTAL PMIDs FOUND: {len(pmid_list)} ---")
    
    # Running at 50% CPU capacity to manage I/O and Network bandwidth
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(process_one_paper, pmid): pmid for pmid in pmid_list}
        for future in as_completed(futures):
            try:
                res = future.result()
                if res: print(res)
            except Exception as e:
                pass