import os
import time
import random
from playwright.sync_api import sync_playwright

VIP_RECOVERY = [
    {
        "name": "Sanjana_Cell_2020_Optimization",
        "url": "https://www.cell.com/cell/fulltext/S0092-8674(20)30607-4",
    },
    {
        "name": "Sanjana_Nature_2018_CasRx",
        "url": "https://www.nature.com/articles/s41586-018-0410-1",
    },
    {
        "name": "Lei_Nature_2021_DeepCas13",
        "url": "https://www.nature.com/articles/s41587-021-00922-2",
    },
    {
        "name": "Sanjana_Massively_Parallel_Cas13",
        "url": "https://www.nature.com/articles/s41467-020-19022-y",
    }
]

PAPER_DIR = "paper_Cas13"

def recover_vips():
    print(f"--- 🚑 Starting VIP Recovery v3 (Ultra-Stealth) ---")
    
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        
        for target in VIP_RECOVERY:
            folder = os.path.join(PAPER_DIR, target["name"])
            os.makedirs(folder, exist_ok=True)
            
            print(f"\n🚀 Targeting: {target['name']}")
            # High-fidelity desktop fingerprint
            context = browser.new_context(
                user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/121.0.0.0 Safari/537.36",
                viewport={'width': 1920, 'height': 1080}
            )
            page = context.new_page()
            
            try:
                # Use domcontentloaded to avoid the 'networkidle' timeout trap
                page.goto(target["url"], wait_until="domcontentloaded", timeout=60000)
                print("   ⏱️ Page loaded. Mimicking human reading time...")
                time.sleep(random.uniform(5, 8))

                # Scroll in small, erratic increments to look human
                for i in range(3):
                    page.mouse.wheel(0, random.randint(400, 800))
                    time.sleep(random.uniform(1, 2))

                # Target specific Nature/Cell supplementary containers
                # Nature: #supplementary-information
                # Cell: .supplementary-material
                found = False
                
                # Broad hunt for Excel/CSV first
                selectors = [
                    "a[href*='.xlsx']", "a[href*='.xls']", "a[href*='.csv']",
                    "a:has-text('Supplementary Table')", "a:has-text('Source Data')"
                ]

                for sel in selectors:
                    elements = page.locator(sel)
                    if elements.count() > 0:
                        print(f"   📊 Found data link: {sel}")
                        with page.expect_download(timeout=30000) as download_info:
                            elements.first.click()
                        download = download_info.value
                        dest = os.path.join(folder, download.suggested_filename)
                        download.save_as(dest)
                        print(f"   ✅ SECURED DATA: {dest}")
                        found = True
                        break
                
                if not found:
                    print("   📄 Data link hidden. Attempting to Print Full Page to PDF...")
                    page.pdf(path=os.path.join(folder, "full_content_snapshot.pdf"))

            except Exception as e:
                print(f"   ⚠️ Failure: {e}")
            
            context.close()
            time.sleep(random.uniform(15, 30)) # Major cooldown to reset IP reputation

        browser.close()

if __name__ == "__main__":
    recover_vips()