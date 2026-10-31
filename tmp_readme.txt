# just a tmp readme for some number field labels scripts to run on Legendre:


# to generate field labels for degree 3
PYTHONPATH=~/sage-extra sage -python scripts/name_lmfdb_fields.py --degrees 3 --source db --only-unnamed --batch 100

# to generate field labels for degreee 4, only for Galois group 4T2
PYTHONPATH=~/sage-extra sage -python scripts/name_lmfdb_fields.py --degrees 4 --galois 4T2 --source db --only-unnamed --batch 100

# to run all the tests:
sage -python -m unittest discover tests -v

