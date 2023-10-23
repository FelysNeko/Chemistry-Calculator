from api.chemistry import lookup, Molecule, Equation
''' Uncomment the line(s) you wish to try '''



# ''' Use lookup function to access Atom class '''
# atom = 'Na'
# print(lookup(atom).info)
# print(lookup(atom) == lookup(atom))
# print(lookup('unknown').info)



# ''' Molecule class basic uses '''
# molecule = 'NaCl'
# print(Molecule(molecule).rectify(True).info)
# print(Molecule(molecule).rectify(True).mass)
# print(Molecule(molecule).rectify(True).solubility)
# print(Molecule(molecule).rectify(True).bond)
# print(Molecule(molecule).rectify(True).count)

# ''' More details for rectify method '''
# print(Molecule('MnCl3').rectify(True).info)
# print(Molecule('MgCl').rectify(False).short)

# ''' Advanced feature '''
# print(Molecule.null('0').info)
# bad = Molecule('unknown').rectify(True)
# print(Molecule.isna(bad))
# print(Molecule(molecule).rectify(True) == Molecule(molecule).rectify(True))



# ''' Equation class basic uses '''
# equation = ['Na(OH)', 'HCl']
# print(Equation(equation).rectify(True).info)
# print(Equation(equation).rectify(True).short)
# print(Equation(equation).rectify(True).count)
# print(Equation(equation).rectify(True).reaction)
# print(Equation(equation).rectify(True).prediction.short)

# ''' More details for balance method '''
# re ,pr = Equation(equation).rectify(True).balance(manual=False)
# print(re.short, pr.short)
# re ,pr = Equation(equation).rectify(True).balance(manual=['H2O', 'NaCl'])
# print(re.short, pr.short)

# ''' Advanced feature, please notice the last line does not hold '''
# print(Equation.null('0').info)
# bad = Equation(['unknown', 'again?']).rectify(True)
# print(Equation.isna(bad))
# print(Equation(equation).rectify(True) == Equation(equation).rectify(True))



''' PLAYGROUND START HERE '''
