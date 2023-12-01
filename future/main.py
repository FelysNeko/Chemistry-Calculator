from chemistry import Table, Equation, Molecule, balance, lookup


print(Table.periodic)
print(Table.anion)
print(Table.cation)
print(Table.diatomic)
print(Table.activity)


atom = lookup('OH')
print(atom)
print(repr(atom))
print(atom.count)
print(atom.isnull())
print(atom.ispoly())


molecule = Molecule('Mg(OH)2')
print(molecule)
print(repr(molecule))
print(molecule.bond)
print(molecule.mass)
print(molecule.solubility)
print(molecule.count)
print(molecule.isnull())
print(Molecule('MnCl4', 0).rectify(True))
print(Molecule('MnCl4', 0).rectify(False))
print(Molecule.null())


equation = Equation('NaCl')
print(equation)
print(repr(equation))
print(equation.count)
print(equation.isnull())
print(Molecule.null())


re = Equation('C6H12O6', 'O2')
pr = Equation('H2O', 'CO2')
print(*balance(re, pr), sep=' -> ')