import bs4
import requests
import os
import re
import numpy as np
import sys

base_url = "https://data.sdss.org/sas/dr17/apogee/spectro/speclib/atmos/marcs/MARCS_v3_2016/"

folder = "mod_z-2.50/"

save_root = "/Volumes/escala_external_drive/apogee/marcs/"

save_folder = f'{save_root}{folder}'

if not os.path.exists(save_folder):
   os.mkdir(save_folder)

r = requests.get(f'{base_url}{folder}')
data = bs4.BeautifulSoup(r.text, "html.parser")

for l in data.find_all("a"):

  try:

      tint = re.search(r'[sp]\d\d\d\d', l["title"]).group()[1:]
      tint = np.float64(tint)

      gint = re.search(r'g\+\d\.\d', l["title"]).group()[2:]
      gint = np.float64(gint)

      cint = re.search(r'c[\+-]\d\.\d\d', l["title"]).group()[1:]
      cint = np.float64(cint)

      print(tint, gint)
      sys.stdout.flush()

  except:
      continue

  if (tint < 3000) or (tint > 6000) or (cint != 0.) or (gint < 0.) or (gint > 5.):
        continue

  if os.path.exists(f'{save_folder}{l["title"]}.gz') or\
     os.path.exists(f'{save_folder}{l["title"]}'):
      continue

  r = requests.get(f'{base_url}{folder}{l["title"]}')

  if r.status_code == 200:

    with open(f'{save_folder}{l["title"]}', 'wb') as f:
       f.write(r.content)
