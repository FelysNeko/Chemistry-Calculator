# PROGRAM IS UNDER RECONSTRUCTION
diatomic = ['O', 'N', 'H', 'F', 'Cl', 'I', 'Br']
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
            return None

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

    def rectify(self):
        if self.bond == 'ionic':
            if self.atoms[0].charge[0] > 0:
                self.atoms[0], self.atoms[1] = self.atoms[1], self.atoms[0]
            charge = abs(self.atoms[0].quantity * self.atoms[0].charge[0])
            for i in range(len(self.atoms[1].charge)):
                if self.atoms[1].charge[i] * self.atoms[1].quantity == charge:
                    self.atoms[1].charge = [self.atoms[1].charge[i]]
                    self.bond = '*ionic'
                    break

    def form(self):
        if self.bond == 'ionic':
            lcm = abs(self.atoms[0].quantity * self.atoms[0].charge[0])
            while lcm % self.atoms[0].charge[0] or lcm % self.atoms[1].charge[0]:
                lcm += 1
            self.atoms[0].quantity = int(abs(lcm / self.atoms[0].charge[0]))
            self.atoms[1].quantity = int(abs(lcm / self.atoms[1].charge[0]))
            self.bond = '*ionic'


class Equation():
    def __init__(self, *args) -> None:
        self.comp = []
        self.type = None

        for each in args:
            molecule = Molecule(each)
            self.comp.append(molecule)


reactant = Equation('Na(OH)', 'HCl')
