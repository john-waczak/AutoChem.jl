CASAPIScrape.jl

using CSV
using DataFrames
using HTTP
using JSON

msl = "/Users/jsadler/Desktop/Julia/AutoChem.jl/src/species/master_species_list.csv"
csv_msl = CSV.File(msl)

msl_df = CSV.read("/Users/jsadler/Desktop/Julia/AutoChem.jl/src/species/master_species_list.csv", DataFrame)

url = "https://commonchemistry.cas.org/api"

r = HTTP.get(url)

molecule_name = "Formaldehyde"
query_url = "$url/search?q=$molecule_name"
response = HTTP.get(query_url)

using Gumbo

if response.status == 200
else
    println("Error: HTTP request failed with status code {response.status}")
end

parsed_html = parsehtml(String(response.body))

using Cascadia

links = eachmatch(sel"your_css_selector_here", parsed_html.root)

for link in links
    link_text = text(link)
    link_href = get(link.attributes, "href", "")

    println("Link Text: $link_text")
    println("Link Href: $link_href")
end




