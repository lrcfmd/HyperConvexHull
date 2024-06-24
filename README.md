# HyperConvexHull
# Computation of energy above the convex hull
# for test compositions
# w.r.t. relevant entries in MaterialsProject (MP)
# or prepared references

1. Unzip MP_composition_etotal.pickle.zip 

2. Prepare test compositions file in .csv format with two columns
'composition' and 'energy'

3. (Optional) Prepare references file in .csv format with two columns
'composition' and 'e_total' or 'energy'

4. Run 
python testCH.py -test your_test_file.csv
