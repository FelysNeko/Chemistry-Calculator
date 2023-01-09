# if there exist polyatomic ions in a compound, always use brackets to hold them
from copy import deepcopy

path = 'data.csv'
table = {}
diatomic = ['O', 'N', 'H', 'F', 'Cl', 'I', 'Br']
activity = ['Li', 'K', 'Ba', 'Sr', 'Ca', 'Na', 'Mg', 'Al', 'Zn', 'Cr', 'Fe', 'Cd', 'Co', 
            'Ni', 'Sn', 'Pb', 'H', 'Cu', 'Hg', 'Ag', 'Pt', 'Au', 'F', 'Cl', 'Br', 'I']

# load periodic table
with open(path) as file:
    for each in file.readlines():
        line = each.strip('\n').split(',')
        for i in range(1, len(line)-1):
            line[i] = float(line[i])
        table[line[0]] = line[1:]
print('Database Loaded')


# look up the given atomic and return the parameter
def lookup(name=str) -> list:
    for element, data in table.items():
        if name == element:
            return [element, data]
    print('Atom or ployatomic ion not found')


# read the coefficients in molecole and attach each atom's parameter
def expand(molecule=str) -> list:
    if molecule[0].isdigit() is False:
        molecule = '1' + molecule

    result = []
    # you kick out the polyatomic ions first and then the remainings
    while '(' in molecule or molecule.isdigit() is False:
        if '(' in molecule:
            begin = molecule.find('(')
            end = molecule.find(')')
            information = lookup(molecule[begin + 1:end])
        else:
            for each in molecule:
                if each.isupper() is True:
                    break
            begin = molecule.find(each)
            end = begin + 1 if molecule[-1] != each and molecule[begin + 1].islower() is True else begin
            information = lookup(molecule[begin:end + 1])

        # append the data into result list
        if information is not None:
            quantity = 1
            for i in range(end + 1, len(molecule) + 1):
                if i == len(molecule) or molecule[i].isdigit() is False:
                    break
            if molecule[end + 1:i] != '':
                quantity = int(molecule[end + 1:i])
        else:
            print('Failed to expand the given list')
            return None

        molecule = molecule[:begin] + molecule[i:]
        result.insert(0, [quantity, information])

    result.insert(0, int(molecule))
    return result


# shrink the expanded format into readable format
def shrink(standard=list) -> str:
    molecule = ''
    if standard[0] != 1:
        molecule = str(standard[0])
    for i in range(1, len(standard)):
        if standard[i][1][1][-1] == '0':
            molecule = f'({standard[i][1][0]})' + molecule
        else:
            molecule = f'{standard[i][1][0]}' + molecule
        if standard[i][0] != 1:
            molecule = str(standard[i][0]) + molecule
    return molecule


# determine the interatomic bonds based on whether they are metal/non-mental
def interatomic(molecule=str) -> str:
    molecule = expand(molecule)
    flag = 0

    if molecule is None:
        return None
    for i in range(1, len(molecule)):
        if molecule[i][1][0] == 'C' or molecule[i][1][0] == 'H':
            flag += 1
    if flag >= 2:
        return 'organic'
    elif len(molecule) < 3:
        return None
    elif molecule[1][1][1][2] * molecule[2][1][1][2] > 0:
        return 'covalent'
    elif molecule[1][1][1][2] * molecule[2][1][1][2] < 0:
        return 'ionic'


# determine the reaction type, it does not cover all of them
def reaction(reactant=list) -> str:
    expansion = list(map(expand, reactant))
    if None in expansion:
        return None
    elif len(expansion) == 1 and len(expansion[0]) == 3:
        return 'decomposition'
    elif len(expansion) == 2:
        for i in range(2):
            if 'O2' in reactant[0-i] and len(expansion[0-i]) == 2 and interatomic(reactant[1-i]) == 'organic':
                return 'combustion'
            elif len(expansion[0]) + len(expansion[1]) == 4 and expansion[0][1][1][1][2] * expansion[1][1][1][1][2] < 0:
                return 'combination'
            elif len(expansion[0]) + len(expansion[1]) == 5:
                if len(expansion[0-i]) == 2 and interatomic(reactant[1-i]) == 'ionic':
                    return 'single'
            elif len(expansion[0]) == 3 and len(expansion[1]) == 3:
                for a in range(1, 3):
                    for b in range(1, 3):
                        if expansion[0-i][a][1][0] == 'OH' and expansion[1-i][b][1][0] == 'H':
                            return 'neutralization'
                if interatomic(reactant[0]) == 'ionic' and interatomic(reactant[1]) == 'ionic' and i == 1:
                    return 'double'


# predict the product
def predict(reactant=list):
    def form(raw):
        if len(raw) == 2:
            if raw[1][1][0] in diatomic:
                return [1, [2, raw[1][1][0]]]
            else:
                return [1, [1, raw[1][1][0]]]
        elif interatomic(shrink(raw)) == 'ionic':
            i = abs(raw[1][1][1][2]) if abs(raw[1][1][1][2]) > abs(raw[2][1][1][2]) else abs(raw[2][1][1][2])
            while True:
                if i % raw[1][1][1][2] == 0 and i % raw[2][1][1][2] == 0:
                    break
                else:
                    i += 1
            return [1, [int(i/abs(raw[1][1][1][2])), raw[1][1][0]], [int(i/abs(raw[2][1][1][2])), raw[2][1][0]]]

    def decompose(expansion):
        short = expansion[0][1:]
        for i in range(2):
            if short[0-i][1][0] == 'H' and 'O' in short[1-i][1][0]:
                return [[1, [1, 'O'], [2, 'H']], [1, [2, 'O'], [1, 'C']]]
            elif short[0-i][1][1][2] > 0 and 'ClO3' in short[1-i][1][0]:
                return [form(expand(f'{short[0-i][1][0]}Cl')), [1, [2, 'O']]]
            elif short[0-i][1][1][2] > 0 and 'OH' in short[1-i][1][0]:
                return [form(expand(f'{short[0-i][1][0]}O')), [1, [1, 'O'], [2, 'H']]]
            elif short[0-i][1][1][2] > 0 and 'CO3' in short[1-i][1][0]:
                return [form(expand(f'{short[0-i][1][0]}O')), [1, [2, 'O'], [1, 'C']]]
            elif short[0-i][1][1][0] * short[1-i][1][1][0] != 0:
                return [form([1, [1, short[0-i][1]]]), form([1, [1, short[1-i][1]]])]

    def combine(expansion):
        return [form([1, [1, expansion[0][1][1]], [1, expansion[1][1][1]]])]
    
    def combust(expansion):
            product = [[1, [1, 'O'], [2, 'H']], [1, [2, 'O'], [1, 'C']]]
            for i in range(2):
                for a in range(1, len(expansion[i])):
                    if expansion[i][a][1][0] == 'S':
                        product.append([1, [2, 'O'], [1, 'S']])
                        print('Note: does not support this equation balancing')
                        break
            return product

    def neutralize(expansion):
        molecule = ''
        for i in range(2):
            molecule += f"({expansion[i][2][1][0]})" if expansion[i][1][1][0] == 'OH' or expansion[i][1][1][0] == 'H' else f'({expansion[i][1][1][0]})'
        return [form(expand(molecule)), [1, [1, 'O'], [2, 'H']]]

    def single(expansion):
        if len(expansion[0]) == 3:
            expansion.reverse()
        product = deepcopy(expansion)
        for i in range(1, 3):
            if expansion[0][1][1][1][2] * expansion[1][i][1][1][2] > 0:
                if expansion[0][1][1][0] in activity and expansion[1][i][1][0] in activity:
                    if activity.index(expansion[0][1][1][0]) < activity.index(expansion[1][i][1][0]):
                        product[0] = [1, expansion[1][i]]
                        product[1][i] = expansion[0][1]
                        break
        return [form(product[0]), form(product[1])]

    def double(expansion):
        product = deepcopy(expansion)
        for a in range(1, 3):
            for b in range(1, 3):
                if expansion[0][a][1][1][2] * expansion[1][b][1][1][2] > 0:
                    product[0][a] = [1, expansion[1][b][1]]
                    product[1][b] = [1, expansion[0][a][1]]
                    return [form(product[0]), form(product[1])]

    kind = reaction(reactant)
    expansion = list(map(expand, reactant))

    # ensure the metals are in correct ions
    for i in range(len(expansion)):
        if interatomic(reactant[i]) == 'ionic':
            if expansion[i][1][1][1][2] > 0:
                expansion[i][1], expansion[i][2] = expansion[i][2], expansion[i][1]
            for a in range(2, 5):
                if expansion[i][1][1][1][2] * expansion[i][1][0] + expansion[i][2][1][1][a] * expansion[i][2][0] == 0:
                    expansion[i][2][1][1][2], expansion[i][2][1][1][a] = expansion[i][2][1][1][a], expansion[i][2][1][1][2]
                    break
            else:
                print('Cannot find the correct ion charges')
                return None

    if kind == 'decomposition':
        return decompose(expansion)
    elif kind == 'combination':
        return combine(expansion)
    elif kind == 'combustion':
        return combust(expansion)
    elif kind == 'neutralization':
        return neutralize(expansion)
    elif kind == 'single':
        return single(expansion)
    elif kind == 'double':
        return double(expansion)
    else:
        print('Unknown Reaction')
        return None


def balance(expression, max=21):
    class Data():
        def __init__(self, formula, count) -> None:
            self.formula = formula
            self.count = count

        # break the polyatomic down and add all elements into the count dictionary
        def refresh(self):
            for molecule in self.formula:
                for i in range(1, len(molecule)):
                    temp = expand(molecule[i][1])
                    if len(temp) > 2:
                        for a in range(1, len(temp)):
                            if temp[a][1][0] not in self.count:
                                self.count[temp[a][1][0]] = temp[a][0] * molecule[i][0] * molecule[0]
                            else:
                                self.count[temp[a][1][0]] += temp[a][0] * molecule[i][0] * molecule[0]
                    else:
                        if molecule[i][1] not in self.count:
                            self.count[molecule[i][1]] = molecule[i][0] * molecule[0]
                        else:
                            self.count[molecule[i][1]] += molecule[i][0] * molecule[0]

        def simplify(self):
            result = []
            for each in self.formula:
                molecule = ''
                for element in each[:0:-1]:
                    if lookup(element[1])[1][-1] == '0':
                        molecule += f'({element[1]})'
                    else:
                        molecule += f'{element[1]}'
                    if element[0] != 1:
                        molecule += str(element[0])
                if each[0] != 1:
                    molecule = str(each[0]) + molecule
                result.append(molecule)
            return result


    product = Data(predict(expression), {})
    reactant = Data(list(map(expand, expression)), {})

    # shrink the reactant formula into short version
    if product.formula is not None and reactant.formula is not None:
        for a in range(len(reactant.formula)):
            for b in range(1, len(reactant.formula[a])):
                reactant.formula[a][b][1] = reactant.formula[a][b][1][0]
    else:
        print('Reactant or predicated product seemd to be incorrect')
        return None

    # current does not support three or more molecules balancing
    flag = False
    for a in range(1, max):
        for b in range(1,max):
            for c in range(1,max):
                for d in range(1,max):
                    
                    #test different sets of number
                    reactant.formula[0][0] = a
                    product.formula[0][0] = b
                    if len(reactant.formula) == 2:
                        reactant.formula[1][0] = c
                    if len(product.formula) == 2:
                        product.formula[1][0] = d

                    # refresh count dictionary and empty it after calculation
                    reactant.refresh()
                    product.refresh()
                    if reactant.count == product.count:
                        flag = True
                        break
                    else:
                        reactant.count = {}
                        product.count = {}

                if flag == True:
                    break
            if flag == True:
                break
        if flag == True:
            break
    else:
        print('Calculation failed, you might need to increase the maxium calculation atempts')
        return None

    return [reactant.simplify(), product.simplify()]


# print('>>>', lookup('U'))
# print('>>>', expand('CH4'))
# print('>>>', interatomic('NaCl'))
print('>>>', reaction(['MgO', 'HCl']))
print('>>>', predict(['C2H2' , 'O2']))
print('>>>', balance(['MgO', 'HCl'], 31))
