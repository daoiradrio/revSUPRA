import os
from SUPRAConformer.structure import Structure
from utils.analyzer import Analyzer

mol1 = Structure(os.sys.argv[1])
mol2 = Structure(os.sys.argv[2])
analyzer = Analyzer()
print(analyzer.doubles_advanced(mol1.coords, mol1.bond_partners, mol2.coords, mol2.bond_partners))
