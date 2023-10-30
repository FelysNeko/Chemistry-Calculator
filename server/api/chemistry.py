from .data import *
import sympy as sp



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

    
    def __eq__(self, __value:object) -> bool:
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
        return False
    

    @property
    def count(self):
        counter = {}
        for each in self.atoms:
            for key, value in each.count.items():
                counter[key] = counter[key]+value if key in counter else value
        result = {key:value*self.coef for key, value in counter.items()}
        return result
    

    def has(self, atom:str) -> bool:
        for each in self.atoms:
            if each.symbol == atom:
                return True
        else:
            return False


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
        elif self.bond=='single' and self.atoms[0].quantity==1 and self.atoms[0].symbol in Table.diatomic:
            self.atoms[0].quantity = 2
        return self
    

    @staticmethod
    def null(value:str) -> object:
        return Molecule(value)
    

    @staticmethod
    def isna(x:object) -> bool:
        return True if isinstance(x, Molecule) and not len(x.atoms) else False



class Equation:
    def __init__(self, equation:list) -> None:
        self.molecule = [deepcopy(Molecule(i)) for i in equation]
        

    def rectify(self, flag:bool=True) -> object:
        for each in self.molecule:
            each.rectify(flag)
        return self


    @property
    def info(self) -> list:
        return [i.info for i in self.molecule]
    

    @property
    def short(self) -> list:
        return [i.short for i in self.molecule]
    

    @property
    def reaction(self) -> str:
        if len(self.molecule)==1 and len(self.molecule[0].atoms)==2:
            return 'decomposition'
        elif len(self.molecule) == 2:
            if self.molecule[0].bond == self.molecule[1].bond == 'single':
                charge = self.molecule[0].atoms[0].charge.head * self.molecule[1].atoms[0].charge.head
                return '*combination' if charge<0 else 'combination'
            elif self.molecule[0].bond == self.molecule[1].bond == '*ionic':
                check = [self.molecule[0-i].has('OH') and self.molecule[1-i].has('H') for i in range(2)]
                return 'neutralization' if True in check else 'double'
            elif sum([self.molecule[0-i].has('O') and self.molecule[1-i].bond=='organic' for i in range(2)]):
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
                return '*single' if check>0 else 'single'
        return 'unknown'
    

    @property
    def prediction(self) -> object:
        if self.reaction == 'decomposition':
            product = [i.symbol for i in self.molecule[0].atoms]
            return Equation(product).rectify(False)
        elif 'combination' in self.reaction:
            product = [''.join([i.atoms[0].symbol for i in self.molecule])]
            return Equation(product).rectify(False)
        elif self.reaction == 'combustion':
            check = sum(['S' in i for i in self.short])
            product = ['H2O', 'CO2', 'SO2'] if check else ['H2O', 'CO2']
            return Equation(product).rectify(False)
        elif self.reaction == 'neutralization':
            neg, pos = self.molecule[0].atoms[0].symbol, self.molecule[1].atoms[1].symbol
            temp = self.molecule[1].atoms[0].symbol+self.molecule[0].atoms[1].symbol if 'H' in neg else neg+pos
            product = ['H2O', temp]
            return Equation(product).rectify(False)
        elif self.reaction=='double' or self.reaction=='*single':
            product = deepcopy(self)
            product.molecule[0].atoms[0], product.molecule[1].atoms[0] = product.molecule[1].atoms[0], product.molecule[0].atoms[0]
            return product.rectify(False)
        else:
            return Equation.null('0')
        

    @property
    def count(self) -> list: 
        return [i.count for i in self.molecule]
    

    def balance(self, manual:bool=False) -> tuple:
        reactant = deepcopy(self)
        product = deepcopy(Equation(manual)) if manual else reactant.prediction
        count = reactant.count + product.count

        key = {k for i in count for k in i.keys()}
        m = sp.Matrix([[x[s] if s in x else 0 for x in count] for s in key])
        result = sp.linsolve(m)

        if len(result) == 1:
            result = list(result)[0]
        else:
            return self, Equation.null('0')

        multiple = sp.lcm([i.denominator for i in result])
        coef = sp.Array(result) * multiple

        for i in range(len(reactant.molecule)):
            reactant.molecule[i].coef = coef[i]
        for i in range(len(product.molecule)-1):
            product.molecule[i].coef = abs(coef[i+len(reactant.molecule)])

        return reactant, product


    @staticmethod
    def null(value:str) -> object:
        return Equation([value])
    

    @staticmethod
    def isna(x:object) -> bool:
        return True if isinstance(x, Equation) and sum([Molecule.isna(i) for i in x.molecule]) else False
    