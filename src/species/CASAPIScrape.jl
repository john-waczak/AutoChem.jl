CASAPIScrape.jl

using HTTP

url = "https://commonchemistry.cas.org/api"

r = HTTP.get(url)

using Gumbo

h = parsehtml(String(r.body))

h.root[1]

using Cascadia

