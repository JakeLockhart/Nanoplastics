import requests

# Example list of UniProt IDs
protein_ids = [
"Q3TGR2",
"E9PV24",
"Q61703",
"A0A338P7H5",
"P06684",
"Q8BH35",
"Q7M084",
"A0A0R4J039",
"G3X8T9",
"Q5FW62",
"A8DUP0",
"B7ZNS9",
"Q00623"
]



# Function to get protein info
def get_protein_info(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        
        # Extract protein name and description (if available)
        protein_name = data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "N/A")
        description = data.get("comments", [])
        function_description = "N/A"
        for comment in description:
            if comment.get("commentType") == "FUNCTION":
                function_description = comment.get("texts", [{}])[0].get("value", "N/A")
                break

        return {
            "id": uniprot_id,
            "name": protein_name,
            "function": function_description
        }
    else:
        return {"id": uniprot_id, "name": "N/A", "function": "Not found or error"}

# Process all protein IDs
for pid in protein_ids:
    info = get_protein_info(pid)
    print(f"{info['id']}: {info['name']}\nFunction: {info['function']}\n")
