'''
The WMDA data is part of the supplemental materials from the following
article.

Bochtler W., Gragert L., Patel Z.I., Robinson J., Steiner D., Hofmann J.A.,
Pingel J., Baouz A., Melis A., Schneider J., Eberhard H.P., Oudshoorn M.,
Marsh S.G., Maiers M. & Müller C.R. (2016).
A comparative reference study for the validation of HLA‐matching algorithms
in the search for allogeneic hematopoietic stem cell donors and cord blood
units.
HLA. 2016 Jun;87(6):439-48. doi: 10.1111/tan.12817.

It is distributed under a Creative Commons Attribution‐NonCommercial‐NoDerivs
(CC BY-NC-ND 4.0) license.

https://creativecommons.org/licenses/by-nc-nd/4.0/
'''

import os.path
import urllib.request
from tarfile import TarFile

url='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5089599/bin/TAN-87-439-s006.tgz'
file_name='TAN-87-439-s006.tgz'

# download file
if os.path.isfile(file_name):
    print('Skipping download step.  File found: ', file_name)
else:
    print('Downloading file', file_name, '...')
    urllib.request.urlretrieve(url, file_name)

# unpack it
print('Unpacking file ...')
TarFile.open(file_name).extractall(path="./wmda")

print('Done.')
