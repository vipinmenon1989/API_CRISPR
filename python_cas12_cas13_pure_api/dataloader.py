import os
import json
import time
import random
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import Entrez
from google import genai
from google.genai import types
from google.api_core import exceptions
from playwright.sync_api import sync_playwright
from tenacity import retry, stop_after_attempt, wait_exponential_jitter, retry_if_exception_type

# ==========================================
# 0. CRITICAL: KILL NESTED PARALLELISM
# ==========================================
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

# ==========================================
# CONFIGURATION
# ==========================================
GEMINI_API_KEY = "AIzaSyCtfb3_K9uwCrAxGQdKl-DOi_sCmkSlwms"
NCBI_API_KEY = "350c406985265581bb3a3b9c8ba345a88608"
INSTITUTIONAL_EMAIL = "vmenon@som.umaryland.edu"

Entrez.email = INSTITUTIONAL_EMAIL
Entrez.api_key = NCBI_API_KEY

MAX_WORKERS = 50 
PAPER_DIR = "paper"
DB_FILE = "processed_papers.json"

# ==========================================
# PHASE 0: SURGICAL STRIKE (Browser Edition)
# ==========================================
VIP_TARGETS = [
    {
        "folder": "Nature_Biotech_Kim_2018_DeepCpf1", 
        "url": "https://www.nature.com/articles/nbt.4061", # Main Page (Not File Link)
        "file_hint": ".xlsx" 
    },
    {
        "folder": "Nature_Methods_Kim_2017_InVivo",
        "url": "https://www.nature.com/articles/nmeth.4104",
        "file_hint": ".xlsx"
    },
    {
        "folder": "Nature_Biotech_Kleinstiver_2016_Specificity",
        "url": "https://www.nature.com/articles/nbt.3620",
        "file_hint": ".xlsx"
    },
    {
        "folder": "Nature_Biotech_Bin_2018_Cpf1Screen",
        "url": "https://www.nature.com/articles/nbt.4172",
        "file_hint": ".xlsx"
    }
]

def run_surgical_strike():
    print(f"--- ⚔️ PHASE 0: Surgical Strike (Browser Mode) on {len(VIP_TARGETS)} VIP Targets ---")
    
    with sync_playwright() as p:
        # Launch ONE browser for all VIPs (Faster)
        browser = p.chromium.launch(
            headless=True, 
            args=["--disable-blink-features=AutomationControlled", "--no-sandbox", "--disable-gpu"]
        )
        
        for target in VIP_TARGETS:
            folder_path = os.path.join(PAPER_DIR, target["folder"])
            os.makedirs(folder_path, exist_ok=True)
            
            # Check if exists
            if len(os.listdir(folder_path)) > 0:
                print(f"   ⏩ VIP Skipped (Exists): {target['folder']}")
                continue

            print(f"   🔎 Hunting VIP: {target['folder']}...")
            context = browser.new_context(user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36")
            page = context.new_page()
            
            try:
                # 1. Go to Article Page (Establishes Cookies)
                page.goto(target["url"], wait_until="domcontentloaded", timeout=60000)
                
                # 2. Find the Supplementary Link
                # Nature stores these in a "Supplementary information" section
                # We look for ANY link ending in .xlsx or containing "Supplementary"
                try:
                    # Look for xlsx links directly
                    download_promise = page.expect_download(timeout=15000)
                    
                    # Try clicking the first Excel file found
                    found = False
                    
                    # Strategy A: Specific Supp Data Link
                    locators = [
                        f"a[href*='{target['file_hint']}']", # Finds .xlsx
                        "a:has-text('Supplementary Data')",
                        "a:has-text('Supplementary Table')"
                    ]
                    
                    for sel in locators:
                        if page.locator(sel).count() > 0:
                            print(f"      🖱️ Clicking: {sel}")
                            page.locator(sel).first.click()
                            found = True
                            break
                    
                    if found:
                        download = download_promise.value
                        final_path = os.path.join(folder_path, "vip_data.xlsx")
                        download.save_as(final_path)
                        print(f"      ✅ VIP SECURED: {final_path}")
                    else:
                        print(f"      ❌ VIP FAILED: No Excel link found on page.")

                except Exception as e:
                    print(f"      ❌ VIP ERROR (Click): {e}")

            except Exception as e:
                print(f"      ❌ VIP ERROR (Nav): {e}")
            
            context.close()
            time.sleep(2) # Breath between VIPs

        browser.close()

# ==========================================
# DATABASE LOGIC
# ==========================================
def save_to_db(doi, status, n_size=0, year=0):
    try:
        if os.path.exists(DB_FILE):
            with open(DB_FILE, 'r') as f: db = json.load(f)
        else: db = {}
        db[doi] = {"status": status, "n_size": n_size, "year": year, "timestamp": time.ctime()}
        with open(DB_FILE, 'w') as f: json.dump(db, f, indent=4)
    except: pass 

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
        year = article['Journal']['JournalIssue']['PubDate'].get('Year', '2025')
        
        doi_list = [str(i) for i in article.get('ELocationID', []) if i.attributes.get('EIdType') == 'doi']
        doi = doi_list[0] if doi_list else f"PMID_{pmid}"
        
        safe_title = "".join([c for c in article['ArticleTitle'] if c.isalpha() or c.isdigit() or c==' ']).rstrip()
        
        folder_name = f"{journal}_{author}_{year}"
        full_path = os.path.join(PAPER_DIR, folder_name)
        
        return full_path, safe_title, doi, year
    except: return None, None, None, None

def analyze_ai_safe(title, doi):
    client = genai.Client(api_key=GEMINI_API_KEY) 
    prompt = f"""
    Analyze if this study has ORIGINAL human Cas12/Cpf1 efficiency data (N > 2000).
    Title: {title} | DOI: {doi}
    CRITICAL: 
    1. IGNORE citation count.
    2. INCLUDE Deep Learning papers if they have training data.
    Return JSON ONLY: {{"is_relevant": bool, "n_size": int, "download_url": "str", "logic": "str"}}
    """
    time.sleep(random.uniform(0.1, 1.5))
    try:
        response = client.models.generate_content(
            model="gemini-2.5-flash-lite",
            contents=prompt,
            config=types.GenerateContentConfig(response_mime_type="application/json")
        )
        return json.loads(response.text)
    except Exception:
        return {"is_relevant": False}

def download_hunter(url, folder):
    if not url or "http" not in url: return False
    os.makedirs(folder, exist_ok=True)
    try:
        with sync_playwright() as p:
            browser = p.chromium.launch(
                headless=True,
                args=["--disable-blink-features=AutomationControlled", "--no-sandbox", "--disable-gpu"]
            )
            context = browser.new_context(user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36")
            page = context.new_page()

            try:
                response = page.goto(url, wait_until="domcontentloaded", timeout=45000)
            except: response = None

            if response:
                ct = response.headers.get('content-type', '').lower()
                if 'application/pdf' in ct or '.pdf' in page.url or 'excel' in ct or 'sheet' in ct:
                    ext = ".xlsx" if "excel" in ct or "sheet" in ct else ".pdf"
                    with open(os.path.join(folder, f"cas12_data{ext}"), "wb") as f:
                        f.write(response.body())
                    browser.close()
                    return True

            selectors = ["a[data-test='pdf-link']", "a[class*='pdf']", "a:has-text('Download PDF')", "a:has-text('Supplementary Material')"]
            for sel in selectors:
                try:
                    if page.locator(sel).first.is_visible():
                        with page.expect_download(timeout=15000) as download_info:
                            page.locator(sel).first.click()
                        download = download_info.value
                        download.save_as(os.path.join(folder, download.suggested_filename))
                        browser.close()
                        return True
                except: continue
            browser.close()
            return False
    except Exception: return False

# ==========================================
# MAIN EXECUTION
# ==========================================

def process_one_paper(pmid):
    folder, title, doi, year = get_metadata_safe(pmid)
    if not title: return None
    if os.path.exists(folder) and len(os.listdir(folder)) > 0:
        return f"⏩ Skipped (Exists): {year} | {title[:30]}..."
    try:
        analysis = analyze_ai_safe(title, doi)
        if analysis.get("is_relevant") and analysis.get("n_size", 0) >= 2000:
            success = download_hunter(analysis['download_url'], folder)
            if success:
                save_to_db(doi, "Downloaded", analysis['n_size'], year)
                return f"✅ DOWNLOADED [{year}] (N={analysis['n_size']}): {folder}"
            else:
                if os.path.exists(folder) and not os.listdir(folder): os.rmdir(folder)
                save_to_db(doi, "Download_Failed", 0, year)
                return f"❌ Download Failed [{year}]: {title[:30]}..."
    except: pass
    return None 

def search_all_years():
    all_ids = []
    years = range(2014, 2027)
    for year in years:
        print(f"--- 📅 Scanning Year: {year} ---")
        query = f'("Cas12a" OR "Cas12 variants" OR "Cpf1" OR "Cas12") AND ("human" OR "homo sapiens") AND {year}[dp]'
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort='relevance')
            record = Entrez.read(handle)
            all_ids.extend(record["IdList"])
            handle.close()
        except: pass
    return list(set(all_ids))

if __name__ == "__main__":
    os.makedirs(PAPER_DIR, exist_ok=True)
    
    # 1. Surgical Strike (Browser Mode)
    run_surgical_strike()
    
    # 2. Main Search
    pmid_list = search_all_years()
    print(f"--- 🚀 PHASE 1: Spinning up {MAX_WORKERS} Workers for {len(pmid_list)} Papers ---")
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(process_one_paper, pmid): pmid for pmid in pmid_list}
        for future in as_completed(futures):
            result = future.result()
            if result: print(result)