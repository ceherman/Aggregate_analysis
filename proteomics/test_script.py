#!/usr/bin/python

from urllib import request
from urllib.parse import urlencode, quote_plus

url = 'http://www.uniprot.org/mapping/'

params = {
    'from':'P_REFSEQ_AC',
    'to':'ACC',
    'format':'tab',
    'query':'NP_511114.2 NP_004295.2 NP_031465.2 XP_004934106.1',
}

data = urlencode(params, quote_via=quote_plus).encode("utf-8")
req = request.Request(url, data)

contact = "cherman@udel.edu"
req.add_header('User-Agent', 'Python %s' % contact)

response = request.urlopen(req)
page = response.read(200000)
print(page)
