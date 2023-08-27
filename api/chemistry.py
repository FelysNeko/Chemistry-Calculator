from copy import deepcopy
from .database import Table
import itertools as it


def lookup(name):
    for atom, data in Table.periodic.items():
        if atom == name:
            return data
        

def expand(molecule):
    if not isinstance(molecule, str) or len(molecule) == 0:
        return None
    if not molecule[0].isdigit():
        molecule = '1' + molecule

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
            return None
        
        molecule = molecule[:begin] + molecule[i:]
        result.insert(0, [quantity, information])

    result.insert(0, int(molecule))
    return result


class Molecule:
    def __init__(self, name) -> None:
        self.atoms = []
        self.coef = 1
        
        expansion = expand(name)
        if expansion is not None:
            for each in expansion[1:]:
                each[1].qtt = each[0]
                self.atoms.append(each[1])
            self.coef = expansion[0]


    @property
    def info(self):
        return [self.coef, [i.info for i in self.atoms]]
    

    @property
    def short(self):
        result = str(self.coef)
        for each in self.atoms:
            atom = f'({each.sym}){str(each.qtt)}'
            result += atom
        return result
    

    @property
    def mass(self):
        return round(self.coef * sum([i.qtt * i.mass for i in self.atoms]), 4)
    

    @property
    def solubility(self):
        if len(self.atoms) == 2 and 'ionic' in self.bond:
            for rules in Table.solubility:
                neg = self.atoms[0].sym
                pos = 'Hg2' if self.atoms[1].sym == 'Hg' and self.atoms[1].qtt == 2 else self.atoms[1].sym
                if self.atoms[0].chrg[0] > 0:
                    neg, pos = pos, neg

                if neg in rules.neg and pos in rules.pos:
                    return True
            else:
                return False


    @property
    def bond(self):
        flag = 0
        for atom in self.atoms:
            if atom.sym == 'C' or atom.sym == 'H':
                flag += 1
        if flag >= 2:
            return 'organic'
        elif len(self.atoms) == 1:
            return 'single'
        elif len(self.atoms) == 2:
            charge = self.atoms[0].chrg[0] * self.atoms[1].chrg[0]
            if charge > 0:
                return 'covalent'
            else:
                total = sum([self.atoms[i].chrg[0] * self.atoms[i].qtt for i in range(2)])
                return '*ionic' if total == 0 else 'ionic'
        

    def rectify(self, flag=True):
        temp = deepcopy(self)
        if 'ionic' in temp.bond:
            if flag is True:
                if temp.atoms[0].chrg[0] > 0:
                    temp.atoms[0], temp.atoms[1] = temp.atoms[1], temp.atoms[0]
                charge = abs(temp.atoms[0].qtt * temp.atoms[0].chrg[0])
                for i in range(len(temp.atoms[1].chrg)):
                    if temp.atoms[1].chrg[i] * temp.atoms[1].qtt == charge:
                        temp.atoms[1].chrg = [temp.atoms[1].chrg[i]]
                        break
            elif flag is False:
                lcm = 1
                while lcm % temp.atoms[0].chrg[0] or lcm % temp.atoms[1].chrg[0]:
                    lcm += 1
                temp.atoms[0].qtt = int(abs(lcm / temp.atoms[0].chrg[0]))
                temp.atoms[1].qtt = int(abs(lcm / temp.atoms[1].chrg[0]))
        elif temp.bond == 'single' and temp.atoms[0].sym in Table.diatomic:
            temp.atoms[0].qtt = 2            
        return temp


class Equation:
    def __init__(self, equation) -> None:
        self.comp = []
        self.counter = {}

        for each in equation:
            parsed = Molecule(each).rectify(True)
            self.comp.append(parsed)


    @property
    def info(self):
        return [i.info for i in self.comp]
    

    @property
    def short(self):
        return [i.short for i in self.comp]
    

    @property
    def reaction(self):
        if len(self.comp) == 1 and len(self.comp[0].atoms) == 2:
            return 'decom'
        elif len(self.comp) == 2:
            if self.comp[0].bond == 'single' and self.comp[1].bond == 'single':
                charge = self.comp[0].atoms[0].chrg[0] * self.comp[1].atoms[0].chrg[0]
                return '*comb' if charge < 0 else 'comb'
            
            elif self.comp[0].bond == '*ionic' and self.comp[1].bond == '*ionic':
                for i in range(2):
                    if 'H' in self.short[0-i] and '(OH)' in self.short[1-i]:
                        return 'neutr'
                else:
                    return 'double'
                
            for i in range(2):
                if self.comp[0-i].bond == 'organic' and self.comp[1-i].atoms[0].sym == 'O':
                    return 'combust'
                elif self.comp[0-i].bond == 'single' and self.comp[1-i].bond == '*ionic':
                    return 'single'
                
                
    def refresh(self):
        self.counter.clear()
        for molecule in self.comp:
            for atom in molecule.atoms:
                if atom.structure is None:
                    polyatomic = expand(atom.sym)
                    for each in polyatomic[1:]:
                        count = each[0] * atom.qtt * molecule.coef
                        if each[1].sym not in self.counter:
                            self.counter[each[1].sym] = count
                        else:
                            self.counter[each[1].sym] += count
                else:
                    count = atom.qtt * molecule.coef
                    if atom.sym not in self.counter:
                        self.counter[atom.sym] = count
                    else:
                        self.counter[atom.sym] += count


    def predict(self):
        def decompose():
            a = f'{self.comp[0].atoms[0].sym}{self.comp[0].atoms[0].qtt}'
            b = f'{self.comp[0].atoms[1].sym}{self.comp[0].atoms[1].qtt}'
            product = Equation([a, b])
            return product


        def combine():
            product = Equation([f'({self.comp[0].atoms[0].sym})({self.comp[1].atoms[0].sym})'])
            return product
        
        
        def combust():
            for molecule in self.short:
                product = Equation(['H2O', 'CO2', 'SO2']) if 'S' in molecule else Equation(['H2O', 'CO2'])
            return product
                

        def single():
            if len(self.comp[0].atoms) > len(self.comp[1].atoms):
                self.comp[0], self.comp[1] = self.comp[1], self.comp[0]
            product = deepcopy(self)
            if Table.activity.index(self.comp[0].atoms[0].sym) < Table.activity.index(self.comp[1].atoms[1].sym):
                product.comp[0].atoms[0], product.comp[1].atoms[1] = product.comp[1].atoms[1], product.comp[0].atoms[0]
            return product
        

        def neutralize():
            if self.comp[0].atoms[0].sym == 'OH':
                salt = f'({self.comp[0].atoms[1].sym})({self.comp[1].atoms[0].sym})'
            else:
                salt = f'({self.comp[0].atoms[0].sym})({self.comp[1].atoms[1].sym})'
            product = Equation([salt, 'H2O'])
            return product
            

        def double():
            product = deepcopy(self)
            product.comp[0].atoms[0], product.comp[1].atoms[0] = product.comp[1].atoms[0], product.comp[0].atoms[0]
            return product
        
        if self.reaction == 'decom':
            product = decompose()
        elif self.reaction == '*comb':
            product = combine()
        elif self.reaction == 'combust':
            product = combust()
        elif self.reaction == 'single':
            product = single()
        elif self.reaction == 'neutr':
            product = neutralize()
        elif self.reaction == 'double':
            product = double()
        else:
            product = None

        if product is not None:
            for i in range(len(product.comp)):
                product.comp[i] = product.comp[i].rectify(False)

        return product


    def balance(self, manual=False, attemps=40):
        product =  self.predict() if manual == False else Equation(manual)
        reactant = deepcopy(self)
        length = len(reactant.comp) + len(product.comp)

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
        
            if reactant.counter == product.counter:
                return reactant, product
