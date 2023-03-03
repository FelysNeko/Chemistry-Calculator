from copy import deepcopy
import sys, os
import itertools as it


DIATOMIC = ['O', 'N', 'H', 'F', 'Cl', 'I', 'Br']
ACTIVITY = ['Li', 'K', 'Ba', 'Sr', 'Ca', 'Na', 'Mg', 'Al', 'Zn', 'Cr', 'Fe', 'Cd', 'Co', 'Ni', 'Sn', 'Pb', 'H', 'Cu', 'Hg', 'Ag', 'Pt', 'Au', 'F', 'Cl', 'Br', 'I']
PATH = 'data.csv'
TABLE = {}
TEXT = '''NOTE: THIS PROGRAM ONLY SUPPORT CHEM11 STOICHIOMETRY >_<
-> MODE CODE[1]: Stoichiometry calulator (everything in gram)
-> MODE CODE[2]: Balance equation automatically or manually
-> MODE CODE[3]: Show detailed information of molecle(s)
-> MODE CODE[4]: Calculate the molar mass of a molecule
-> EXIT CODE[0]: Quit program
'''


with open(PATH) as file:
    for each in file.readlines():
        line = each.strip('\n').split(',')
        for i in range(1, len(line)-1):
            line[i] = float(line[i])
        TABLE[line[0]] = line[1:]


def lookup(name):
    for element, data in TABLE.items():
        if name == element:
            return [element, data]


def expand(molecule):
    try:
        if not molecule[0].isdigit():
            molecule = '1' + molecule
    except TypeError:
        print('\033[0;31mError: Invalid Input\033[0m')
        return None
    except IndexError:
        print('\033[0;31mError: Empty Input\033[0m')
        return None

    result = []
    while '(' in molecule or not molecule.isdigit():
        if '(' in molecule:
            begin = molecule.find('(')
            end = molecule.find(')')
            information = lookup(molecule[begin+1:end])
        else:
            for each in molecule:
                if each.isupper() is True:
                    break
            begin = molecule.find(each)
            end = begin+1 if len(molecule)-begin > 1 and molecule[begin+1].islower() is True else begin
            information = lookup(molecule[begin:end+1])

        if information is not None:
            quantity = 1
            for i in range(end+1, len(molecule)+1):
                if i == len(molecule) or not molecule[i].isdigit():
                    break
            if molecule[end+1:i] != '':
                quantity = int(molecule[end+1:i])
        else:
            print('\033[0;31mError: Unable to find the atom or polyatomic molecule\033[0m')
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
        self.mass = 0

        try:
            expansion = expand(name)
            self.coef = expansion[0]
            for each in expansion[1:]:
                temp = Element(each[0], each[1][0], each[1][1])
                self.atoms.append(temp)
        except TypeError:
            print('\033[0;31mError: Unable to the read the input\033[0m')
            return None

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


    def show(self):
        print(f'Coefficient: {self.coef} | Bond: {self.bond} | Mass: {self.mass}')
        for atom in self.atoms:
            print(f'\033[1m{atom.__dict__}\033[0m')


    def weigh(self):
        self.mass = 0
        for atom in self.atoms:
            self.mass += (atom.mass * atom.quantity)
        self.mass *= self.coef


    def digas(self):
        if self.bond == 'single' and self.atoms[0].name in DIATOMIC:
            self.atoms[0].quantity = 2


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
        self.digas()
        self.weigh()


    def form(self):
        if 'ionic' in self.bond:
            lcm = 1
            while lcm % self.atoms[0].charge[0] or lcm % self.atoms[1].charge[0]:
                lcm += 1
            self.atoms[0].quantity = int(abs(lcm / self.atoms[0].charge[0]))
            self.atoms[1].quantity = int(abs(lcm / self.atoms[1].charge[0]))
            self.bond = '*ionic'
        self.digas()
        self.weigh()


class Equation():
    def __init__(self, equation) -> None:
        self.short = tuple(equation)
        self.comp = []
        self.reaction = None
        self.counter = {}

        try:
            for each in equation:
                molecule = Molecule(each)
                molecule.rectify()
                self.comp.append(molecule)
        except TypeError:
            print('\033[0;31mError: Unable to generate equation\033[0m')


    def update(self):
        output = []
        flag = False
        temp = deepcopy(self)
        for molecule in temp.comp:
            name = str(molecule.coef) if molecule.coef != 1 else ''
            if molecule.atoms[0].charge[0] < 0 or molecule.bond == 'organic' and not flag:
                molecule.atoms.reverse()
                flag = True
            for atom in molecule.atoms:
                name += f'({atom.name})' if atom.structure is None else atom.name
                if atom.quantity != 1:
                    name += str(atom.quantity)
            flag = False
            output.append(name)
        self.short = tuple(output)


    def refresh(self):
        self.counter.clear()
        for molecule in self.comp:
            for atom in molecule.atoms:
                if atom.structure is None:
                    polyatomic = expand(atom.name)
                    for each in polyatomic[1:]:
                        count = each[0] * atom.quantity * molecule.coef
                        if each[1][0] not in self.counter:
                            self.counter[each[1][0]] = count
                        else:
                            self.counter[each[1][0]] += count
                else:
                    count = atom.quantity * molecule.coef
                    if atom.name not in self.counter:
                        self.counter[atom.name] = count
                    else:
                        self.counter[atom.name] += count


    def show(self):
        print(f'{self.short} | {self.reaction} | {self.counter}')
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
            if ACTIVITY.index(self.comp[0].atoms[0].name) < ACTIVITY.index(self.comp[1].atoms[1].name):
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
        else:
            product = None
    
        if product is not None:
            for i in range(len(product.comp)):
                product.comp[i].form()

            update = []
            for molecule in product.comp:
                name = str(molecule.coef)
                for atom in molecule.atoms:
                    name += f'({atom.name}){atom.quantity}'
                update.append(name)
            product.short = tuple(update)
        else:
            print('\033[0;31mError: Unable to determine the reaction type\033[0m')

        return product


def balance(equation, product=None, attemps=40):
    try:
        reactant = Equation(equation)
        product = reactant.predict() if product is None else Equation(product)
        length = len(reactant.comp) + len(product.comp)
    except AttributeError:
        print('\033[0;31mError: Unable to balance the equation\033[0m')
        return None

    count = 0
    for i in it.product(range(1, attemps+1), repeat=length):
        reactant.comp[0].coef = i[0]
        product.comp[0].coef = i[1]
        if length >= 3:
            if len(reactant.comp) == 2:
                reactant.comp[1].coef = i[2]
            else:
                product.comp[1].coef = i[2]
        if length >= 4:
            product.comp[1].coef = i[3]
        if length >= 5:
            if len(reactant.comp) == 3:
                reactant.comp[2].coef = i[4]
            else:
                product.comp[2].coef = i[4]

        reactant.refresh()
        product.refresh()

        count += 1
        percent = count / attemps**length * 100
        print('\r', end='')
        print(f'Progress: {count}/{attemps**length} [{round(percent, 1)}%]', end="")
        sys.stdout.flush()

        if reactant.counter == product.counter:
            for i in range(len(reactant.comp)):
                reactant.comp[i].weigh()
            for i in range(len(product.comp)):
                product.comp[i].weigh()
            reactant.update()
            product.update()
            print(' | \033[0;32mCalculation Finished\033[0m')
            return reactant, product
    else:
        print(' | \033[0;31mCalculation Failed\033[0m')


class Stoichiometry():
    def __init__(self, equation):
        self.molar = {}
        self.temp ={}
        self.limit = {}
        self.excceed = None
        self.liname = None

        try:
            for a in range(2):
                for b in range(len(equation[a].comp)):
                    self.molar[equation[a].short[b]] = equation[a].comp[b].mass
        except TypeError:
            print('\033[0;31mError: Unable to load the calculation\033[0m')

    
    def match(self, data):
        for key, value in self.molar.items():
            if expand(data[0])[1:] == expand(key)[1:]:
                return data[1] / value
    

    def calculate(self, parameter):
        if len(parameter) == 1:
            times = self.match(parameter[0])
            for key, value in self.molar.items():
                self.temp[key] = round(value * times, 5)

        elif len(parameter) == 2:
            exratio = self.match(parameter[0])
            liratio = self.match(parameter[1])
            if exratio < liratio:
                liratio = exratio
                parameter[0], parameter[1] = parameter[1], parameter[0]

            for key, value in self.molar.items():
                self.limit[key] = round(value * liratio, 5)
                if expand(key)[1:] == expand(parameter[0][0])[1:]:
                    self.limit[key] = parameter[0][1]
                    self.excceed = (parameter[0][0], parameter[0][1] - round(value * liratio, 5))
            self.liname = parameter[1][0]
        else:
            print('\033[0;31mError: Too many or less input\033[0m')


class Console():
    def get(self):
        while True:
            content = input('>>> Input Value: ').split(' ')
            verify = input(f'>>> Confirmation: {content} (y/n)')
            if verify != 'n':
                return content


    def mass(self):
        parameter = self.get()
        molecule = Molecule(parameter[0])
        molecule.weigh()
        print(f'Mass: \033[1m{round(molecule.mass, 2)}\033[0m')


    def show(self):
        parameter = self.get()
        molecule = Equation(parameter)
        molecule.show()


    def solve(self):
        auto = input('>>> Auto mode? (y/n)')
        if auto != 'n':
            parameter = self.get()
            equation = balance(parameter)
        else:
            reactant = input('>>> Reactant: ').split(' ')
            product = input('>>> Product: ').split(' ')
            equation = balance(reactant, product)
            
        if equation is not None:
            print(f'Equation: \033[1m{equation[0].short} >>> {equation[1].short}\033[0m')
        

    def stoichi(self):
        auto = input('>>> Auto mode? (y/n)')
        if auto != 'n':
            parameter = self.get()
            equation = balance(parameter)
        else:
            reactant = Equation(input('>>> Reactant: ').split(' '))
            product = Equation(input('>>> Product: ').split(' '))
            for i in range(len(reactant.comp)):
                reactant.comp[i].weigh()
            for i in range(len(product.comp)):
                product.comp[i].weigh()
            equation = (reactant, product)

        if equation is not None:
            print(f'Equation: \033[1m{equation[0].short} >>> {equation[1].short}\033[0m')
            new = []
            data = input('>>> Knowns: ').split(' ')
            try:
                for i in range(len(data)//2):
                    new.append((data[i*2], int(data[i*2+1])))
            except ValueError:
                print('\033[0;31mError: Invalid value input\033[0m')

            operation = Stoichiometry(equation)
            operation.calculate(new)
            if len(new) == 1 and len(operation.temp) > 0:
                print(f'Mass Equation: \033[1m{operation.temp}\033[0m')
            elif len(new) == 2:
                print(f'Result: {operation.limit}')
                print(f'Limiting Reactant: \033[1m{operation.liname}\033[0m')
                print(f'[{operation.excceed[0]}] excceed: \033[1m{operation.excceed[1]}\033[0m')
            else:
                print('\033[0;31mCalculation Failed\033[0m')
            
    
def main():
    os.system('clear')
    print(TEXT)

    while True:
        mode = input('>>> Select Mode Code: ')
        console = Console()
        if mode == '1':
            console.stoichi()
        elif mode == '2':
            console.solve()
        elif mode == '3':
            console.show()
        elif mode == '4':
            console.mass()
        elif mode == '0':
            os.system('clear')
            break
        else:
            print('Command not found')

        clear = input('>>> Clear screen? (y/n)')
        if clear != 'n':
            os.system('clear')
            print(TEXT)
        else:
            print('')
            

if __name__ == '__main__':
    main()
    
