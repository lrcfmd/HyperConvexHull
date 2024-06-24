import pymatgen,sys
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from pymatgen.entries.computed_entries import ComputedEntry
#import plotter1
import pandas as pd

def computed(name, energy):
  return ComputedEntry(name, energy)

#df = pd.read_csv('Reference.csv', comment='#') #0
file = sys.argv[1]
df = pd.read_csv(f'{file}', comment='#')

test=[]
ref=[]

for i,j in zip(df.composition, df.e_total):
    test.append(computed(i,j))

PDfull=PhaseDiagram(test)
print(PDfull)

print( "Test new Structures")

for entry in test:
#   if 1000*PDfull.get_e_above_hull(entry) == 0:
      print (entry.composition,',',entry.energy,',', PDfull.get_decomp_and_phase_separation_energy(entry)[1])
#       print (entry.composition,',',entry.energy,',', PDfull.get_decomp_and_phase_separation_energy(entry)[0])
#        print ("%13s,%15.5f,%15.5f," %(entry.composition, entry.energy, 1000*PDfull.get_e_above_hull(entry)))
#        print (entry.composition,',', round(1000*PDfull.get_e_above_hull(entry),1))
#        print (entry.composition, round(1000*PDfull.get_e_above_hull(entry),1))
#       print (entry.composition, entry.energy, 1000*PDfull.get_e_above_hull(entry), PDfull.get_decomp_and_phase_separation_energy(entry))
#        print ("%13s,%15.5f,%15.3f," %(entry.composition, entry.energy, 1000*PDfull.get_e_above_hull(entry)), file=open(f"LiSiSCl_computed_{len(xtalopt+Vlad)}.csv", 'a'))

#plotter = plotter1.get_contour_pd_plot(PDfull)
#plotter.show()
