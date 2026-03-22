import requests

print("📞 Dialing the MyGene Database (NIH-funded backup server)...")

# 1. The web address of our new, highly stable API Waiter
url = "https://mygene.info/v3/query?q=symbol:PIK3CA&fields=symbol,name,map_location"

# 2. We place the order
response = requests.get(url)

# 3. Check for the 200 Success Code
if response.status_code == 200:
    print("✅ Connection successful! Parsing live data...\n")
    
    live_data = response.json() 
    # MyGene returns a list of "hits". We want the first one [0].
    gene_info = live_data['hits'][0] 
    
    print(f"🧬 Gene Symbol: {gene_info.get('symbol')}")
    print(f"📍 Chromosome Band: {gene_info.get('map_location')}")
    print(f"📖 Official Description: {gene_info.get('name')}")

else:
    print(f"❌ Failed to connect. Error Code: {response.status_code}")
