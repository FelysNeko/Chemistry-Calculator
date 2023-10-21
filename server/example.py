import api.chemistry as chem


# print(chem.lookup('Na').info)

# print(chem.Molecule('NaCl').info)

# print(chem.Molecule('NaCl').mass)

# print(chem.Molecule('NaCl').solubility)

# print(chem.Molecule('NaCl').bond)

# print(chem.Molecule('MnCl3').rectify(True).info)
# print(chem.Molecule('MgCl').rectify(False).short)

# print(chem.Molecule('C6H12O6').count)


# print(chem.Equation(['Na(OH)', 'HCl']).rectify(True).info)

# print(chem.Equation(['Na(OH)', 'HCl']).rectify(True).short)

# print(chem.Equation(['Na(OH)', 'HCl']).rectify(True).reaction)

# print(chem.Equation(['Na(OH)', 'HCl']).rectify(True).prediction.short)

# re ,pr = chem.Equation(['O2', 'C6H12O7']).rectify(True).balance(manual=False)
# print(re.short, pr.short)

# re ,pr = chem.Equation(['Na(OH)', 'HCl']).rectify(True).balance(manual=['H2O', 'NaCl'])
# print(re.short, pr.short)

# error
# print(chem.Equation(['ShaBi', 'DianXiaoJiMaLe']).molecule[1] == chem.Molecule.null('1'))
# print(chem.Molecule.null('1').info)

'''
WORKSPACE START HERE
'''
