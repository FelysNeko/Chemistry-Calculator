from copy import deepcopy
from time import time


diatomic = ['O', 'N', 'H', 'F', 'Cl', 'I', 'Br']
activity = ['Li', 'K', 'Ba', 'Sr', 'Ca', 'Na', 'Mg', 'Al', 'Zn', 'Cr', 'Fe', 'Cd', 'Co', 'Ni', 'Sn', 'Pb', 'H', 'Cu', 'Hg', 'Ag', 'Pt', 'Au', 'F', 'Cl', 'Br', 'I']
path = 'data.csv'
table = {}


with open(path) as file:
    for each in file.readlines():
        line = each.strip('\n').split(',')
        for i in range(1, len(line)-1):
            line[i] = float(line[i])
        table[line[0]] = line[1:]


def lookup(name=str) -> list:
    for element, data in table.items():
        if name == element:
            return [element, data]


def expand(molecule=str):
    if molecule[0].isdigit() is False:
        molecule = '1' + molecule

    result = []
    while '(' in molecule or molecule.isdigit() is False:
        if '(' in molecule:
            begin = molecule.find('(')
            end = molecule.find(')')
            information = lookup(molecule[begin+1:end])
        else:
            for each in molecule:
                if each.isupper() is True:
                    break
            begin = molecule.find(each)
            end = begin+1 if len(molecule) > 2 and molecule[begin+1].islower() is True else begin
            information = lookup(molecule[begin:end+1])

        if information is not None:
            quantity = 1
            for i in range(end+1, len(molecule)+1):
                if i == len(molecule) or molecule[i].isdigit() is False:
                    break
            if molecule[end+1:i] != '':
                quantity = int(molecule[end+1:i])
        else:
            return False

        molecule = molecule[:begin] + molecule[i:]
        result.insert(0, [quantity, information])

    result.insert(0, int(molecule))
    return result


class Element():
    def __init__(self, quantity, name, parameter) -> None:
        self.name = name
        self.quantity = quantity
        self.number = parameter[0]
        self.mass = parameter[1]
        self.charge = parameter[2:5]
        self.structure = parameter[5]

        while 0 in self.charge and len(self.charge) > 1:
            self.charge.remove(0)
        if self.structure == '0':
            self.structure = None
        

class Molecule():
    def __init__(self, name) -> None:
        self.atoms = []
        self.coef = 1
        self.bond = None
        self.mass = 0

        if name is not None:
            expansion = expand(name)
            self.coef = expansion[0]
            for each in expansion[1:]:
                temp = Element(each[0], each[1][0], each[1][1])
                self.atoms.append(temp)

        flag = 0
        for atom in self.atoms:
            if atom.name == 'C' or atom.name == 'H':
                flag += 1
        if len(self.atoms) == 1:
            self.bond = 'single'
        elif len(self.atoms) == 2:
            charge = self.atoms[0].charge[0] * self.atoms[1].charge[0]
            self.bond = 'ionic' if charge < 0 else 'covalent'
        if flag >= 2:
            self.bond = 'organic'


    def weigh(self):
        self.mass = 0
        for atom in self.atoms:
            self.mass += (atom.mass * atom.quantity)
        self.mass *= self.coef


    def show(self):
        print(f'Coefficient: {self.coef} | Bond: {self.bond} | Mass: {self.mass}')
        for atom in self.atoms:
            print(f'{atom.__dict__}')


    def rectify(self):
        if 'ionic' in self.bond:
            if self.atoms[0].charge[0] > 0:
                self.atoms[0], self.atoms[1] = self.atoms[1], self.atoms[0]
            charge = abs(self.atoms[0].quantity * self.atoms[0].charge[0])
            for i in range(len(self.atoms[1].charge)):
                if self.atoms[1].charge[i] * self.atoms[1].quantity == charge:
                    self.atoms[1].charge = [self.atoms[1].charge[i]]
                    self.bond = '*ionic'
                    break
        elif self.bond == 'single' and self.atoms[0].name in diatomic:
            self.atoms[0].quantity = 2
        self.weigh()

    def form(self):
        if 'ionic' in self.bond:
            lcm = 1
            while lcm % self.atoms[0].charge[0] or lcm % self.atoms[1].charge[0]:
                lcm += 1
            self.atoms[0].quantity = int(abs(lcm / self.atoms[0].charge[0]))
            self.atoms[1].quantity = int(abs(lcm / self.atoms[1].charge[0]))
            self.bond = '*ionic'
        self.weigh()


class Equation():
    def __init__(self, equation) -> None:
        self.short = tuple(equation)
        self.comp = []
        self.reaction = None
        self.counter = {}

        for each in equation:
            molecule = Molecule(each)
            molecule.rectify()
            self.comp.append(molecule)

    def show(self):
        print('-' * 100)
        print(f'{self.short} | {self.reaction}')
        for molecule in self.comp:
            molecule.show()

    def react(self):
        if len(self.comp) == 1 and len(self.comp[0].atoms) == 2:
            self.reaction = 'decompose'
        elif len(self.comp) == 2:
            if self.comp[0].bond == 'single' and self.comp[1].bond == 'single':
                charge = self.comp[0].atoms[0].charge[0] * self.comp[1].atoms[0].charge[0]
                self.reaction = '*combine' if charge < 0 else 'combine'
            for i in range(2):
                if self.comp[0-i].bond == 'organic' and self.comp[1-i].atoms[0].name == 'O':
                    self.reaction = 'combust'
                elif self.comp[0-i].bond == 'single' and self.comp[1-i].bond == '*ionic':
                    self.reaction = 'single'
                elif self.comp[0-i].bond == '*ionic' and self.comp[1-i].bond == '*ionic':
                    if 'H' in self.short[0-i] and '(OH)' in self.short[1-i]:
                        self.reaction = 'neutralize'
                        break
                    else:
                        self.reaction = 'double'
        else:
            return False
                    
    
    def predict(self):
        def decompose():
            a = f'{self.comp[0].atoms[0].name}{self.comp[0].atoms[0].quantity}'
            b = f'{self.comp[0].atoms[1].name}{self.comp[0].atoms[1].quantity}'
            product = Equation([a, b])
            return product


        def combine():
            product = Equation([f'({self.comp[0].atoms[0].name})({self.comp[1].atoms[0].name})'])
            return product
        
        
        def combust():
            for molecule in self.short:
                product = Equation(['H2O', 'CO2', 'SO2']) if 'S' in molecule else Equation(['H2O', 'CO2'])
            return product
                

        def single():
            if len(self.comp[0].atoms) > len(self.comp[1].atoms):
                self.comp[0], self.comp[1] = self.comp[1], self.comp[0]
            product = deepcopy(self)
            if activity.index(self.comp[0].atoms[0].name) < activity.index(self.comp[1].atoms[1].name):
                product.comp[0].atoms[0], product.comp[1].atoms[1] = product.comp[1].atoms[1], product.comp[0].atoms[0]
            return product
        

        def neutralize():
            if self.comp[0].atoms[0].name == 'OH':
                salt = f'({self.comp[0].atoms[1].name})({self.comp[1].atoms[0].name})'
            else:
                salt = f'({self.comp[0].atoms[0].name})({self.comp[1].atoms[1].name})'
            product = Equation([salt, 'H2O'])
            return product
            

        def double():
            product = deepcopy(self)
            product.comp[0].atoms[0], product.comp[1].atoms[0] = product.comp[1].atoms[0], product.comp[0].atoms[0]
            return product


        if self.react() is False:
            product = None
        elif self.reaction == 'decompose':
            product = decompose()
        elif self.reaction == '*combine':
            product = combine()
        elif self.reaction == 'combust':
            product = combust()
        elif self.reaction == 'single':
            product = single()
        elif self.reaction == 'neutralize':
            product = neutralize()
        elif self.reaction == 'double':
            product = double()
    
        if product is not None:
            for molecule in product.comp:
                molecule.form()

            update = []
            for molecule in product.comp:
                name = str(molecule.coef)
                for atom in molecule.atoms:
                    name += f'({atom.name}){atom.quantity}'
                update.append(name)
            product.short = tuple(update)

        return product


def balance(*args, max=20):
    reactant = Equation([each for each in args])
    product = reactant.predict()


if __name__ == '__main__':
    balance('Na(OH)', 'HCl')
