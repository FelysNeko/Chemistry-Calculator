from .data import *
from copy import deepcopy
import itertools as it



def gcd(a:int, b:int) -> int:
    if not a and b:
        return b
    elif a and not b:
        return a
    else:
        return gcd(b, a%b)



def lcm(a:int, b:int) -> int:
    return int(a*b/gcd(a,b))



class Molecule:
    def __init__(self, molecule:str) -> None:
        coef, data = parse(molecule)
        self.coef = coef
        self.atoms = data

    
    def __eq__(self, __value: object) -> bool:
        return (
            self.coef == __value.coef and
            self.atoms == __value.atoms
        )
        

    @property
    def info(self) -> list:
        return [self.coef, [i.info for i in self.atoms]]
    

    @property
    def short(self) -> str:
        num = lambda x: '' if x==1 else str(x)
        par = lambda x: f'({x})' if x in Table.polyatomic else x
        return num(self.coef) + ''.join([''.join([par(i.symbol), num(i.quantity)]) for i in self.atoms])

    
    @property
    def mass(self) -> float:
        return round(self.coef * sum([i.quantity*i.mass for i in self.atoms]), 4)


    @property
    def bond(self) -> str:
        flag = sum([i.symbol=='C' or i.symbol=='H' for i in self.atoms])
        if flag >= 2:
            return 'organic'
        elif len(self.atoms) == 1:
            return 'single'
        elif len(self.atoms) == 2:
            charge = self.atoms[0].charge.head * self.atoms[1].charge.head
            if charge >= 0:
                return 'covalent'
            else:
                total = sum([self.atoms[i].charge.head * self.atoms[i].quantity for i in range(2)])
                return 'ionic' if total else '*ionic'
        else:
            return 'unknown'


    @property
    def solubility(self) -> bool:
        if len(self.atoms)==2 and 'ionic' in self.bond:
            pos, neg = self.atoms[0].symbol, self.atoms[1].symbol
            pos = 'Hg2' if self.atoms[0].symbol=='Hg' and self.atoms[0].quantity==2 else pos
            pos, neg = (neg, pos) if self.atoms[0].charge.head<0 else (pos, neg)

            for rule in Table.solubility:
                if pos in rule.positive and neg in rule.negative:
                    return True
            else:
                return False
        else:
            return False
        

    @staticmethod
    def null(value:str) -> object:
        return deepcopy(Molecule(value))
    

    @property
    def count(self):
        counter = {}
        for each in self.atoms:
            for key, value in each.count.items():
                counter[key] = counter[key]+value if key in counter else value
        result = {key:value*self.coef for key, value in counter.items()}
        return result


    def rectify(self, flag:bool=True) -> object:
        if 'ionic' in self.bond:
            if self.atoms[1].charge.head > 0:
                self.atoms[0], self.atoms[1] = self.atoms[1], self.atoms[0]
            if flag is True:
                charge = self.atoms[1].charge.head * self.atoms[1].quantity / self.atoms[0].quantity
                self.atoms[0].charge.match(charge)
            else:
                x = lcm(abs(self.atoms[0].charge.head), abs(self.atoms[1].charge.head))
                self.atoms[0].quantity = abs(int(x/self.atoms[0].charge.head))
                self.atoms[1].quantity = abs(int(x/self.atoms[1].charge.head))
        elif self.bond=='single' and self.atoms[0].symbol in Table.diatomic:
            self.atoms[0].quantity = 2
        return self



class Equation:
    def __init__(self, equation:list) -> None:
        self.molecule = [deepcopy(Molecule(i)) for i in equation]
        

    def rectify(self, flag:bool=True) -> object:
        for each in self.molecule:
            each.rectify(flag)
        return self


    @staticmethod
    def null(value:str) ->object:
        return deepcopy(Equation(value))


    @property
    def info(self) -> list:
        return [i.info for i in self.molecule]
    

    @property
    def short(self) -> list:
        return [i.short for i in self.molecule]
    

    @property
    def reaction(self) -> str:
        if len(self.molecule)==1 and len(self.molecule[0].atoms)==2:
            return 'demoleculeostion'
        elif len(self.molecule) == 2:
            if self.molecule[0].bond == self.molecule[1].bond == 'single':
                charge = self.molecule[0].atoms[0].charge.head * self.molecule[1].atoms[0].charge.head
                return '*combination' if charge<0 else 'combination'
            elif self.molecule[0].bond == self.molecule[1].bond == '*ionic':
                check = ['(OH)' in self.molecule[0-i].short and 'H' in self.molecule[1-i].short for i in range(2)]
                return 'neutralization' if True in check else 'double'
            elif sum(['O2' in self.molecule[0-i].short and self.molecule[1-i].bond=='organic' for i in range(2)]):
                return 'combustion'
            elif sum([self.molecule[0-i].bond=='single' and self.molecule[1-i].bond=='*ionic' for i in range(2)]):
                if self.molecule[1].bond=='single':
                    self.molecule[0], self.molecule[1] = self.molecule[1], self.molecule[0]
                try:
                    check = (
                        Table.activity.index(self.molecule[1].atoms[0].symbol) -
                        Table.activity.index(self.molecule[0].atoms[0].symbol)
                    )
                except Exception:
                    check = -1
                return '*single' if check > 0 else 'single'
        return 'unknown'
    

    @property
    def prediction(self) -> object:
        if self.reaction == 'demoleculeostion':
            product = [i.symbol for i in self.molecule[0].atoms]
            return deepcopy(Equation(product).rectify(False))
        elif 'combination' in self.reaction:
            product = [''.join([i.atoms[0].symbol for i in self.molecule])]
            return deepcopy(Equation(product).rectify(False))
        elif self.reaction == 'combustion':
            check = sum('S' in i for i in self.short)
            product = ['H2O', 'CO2', 'SO2'] if check else ['H2O', 'CO2']
            return deepcopy(Equation(product).rectify(False))
        elif self.reaction == 'neutralization':
            neg, pos = self.molecule[0].atoms[0].symbol, self.molecule[1].atoms[1].symbol
            temp = self.molecule[1].atoms[0].symbol+self.molecule[0].atoms[1].symbol if 'H' in neg else neg+pos
            product = ['H2O', temp]
            return deepcopy(Equation(product).rectify(False))
        elif self.reaction=='double' or self.reaction=='*single':
            product = deepcopy(self)
            product.molecule[0].atoms[0], product.molecule[1].atoms[0] = product.molecule[1].atoms[0], product.molecule[0].atoms[0]
            return product.rectify(False)
        else:
            return Equation.null('0')
        

    @property
    def count(self):
        counter = {}
        for each in self.molecule:
            for key, value in each.count.items():
                counter[key] = counter[key]+value if key in counter else value
        return counter
    

    def balance(self, manual=False, attemps=40):
        reactant = deepcopy(self.rectify(True))
        product = deepcopy(Equation(manual)) if manual else reactant.prediction
        length = len(reactant.molecule) + len(product.molecule)
        
        for i in it.product(range(1, attemps+1), repeat=length):
            reactant.molecule[0].coef = i[0]
            product.molecule[0].coef = i[1]

            if length >= 3:
                if len(reactant.molecule) == 2:
                    reactant.molecule[1].coef = i[2]
                else:
                    product.molecule[1].coef = i[2]
            if length >= 4:
                product.molecule[1].coef = i[3]
            if length >= 5:
                if len(reactant.molecule) == 3:
                    reactant.molecule[2].coef = i[4]
                else:
                    product.molecule[2].coef = i[4]

            if reactant.count == product.count:
                return reactant, product
            
        else:
            return self, Equation.null('0')
