from __future__ import annotations
from copy import deepcopy
import re
import os



def lookup(molecule:str) -> Atom:
    if molecule in Table.periodic:
        return deepcopy(Table.periodic[molecule])
    else:
        return Atom.null()


def parse(molecule:str) -> tuple:
    poly = re.findall(r'\((.+?)\)(\d*)', molecule)
    molecule = re.sub(r'\((.+?)\)(\d*)', '', molecule)
    atom = re.findall(r'([A-Z][a-z]{0,1})(\d*)', molecule)
    coef = re.findall(r'(^\d*)', molecule)
    
    cvt = lambda x: 1 if x=='' else int(x)
    def shrink(symbol:str, coef:int) -> Atom:
        atom = lookup(symbol)
        atom.quantity = cvt(coef)
        return atom
    
    data = [shrink(i[0], i[1]) for i in atom+poly if i[0] in Table.periodic]
    coef = cvt(coef[0]) if len(data) else 0
    return coef, data



class Charge:
    def __init__(self, data:list) -> None:
        self.__data:list[int] = [int(i) for i in data]

    def __eq__(self, __value:Charge) -> bool:
        return not sum([self.__data[i]-__value.__data[i] for i in range(3)])
    
    def __repr__(self) -> str:
        return str(self.__data)

    @property
    def head(self) -> int:
        return self.__data[0]
    

    def match(self, charge:int) -> int:
        for i, each in enumerate(self.__data):
            if each+charge == 0:
                self.__data[0], self.__data[i] = self.__data[i], self.__data[0]
                return self.head
        else:
            return 0



class Atom:
    def __init__(self, data:list) -> None:
        self.quantity:int = 0
        self.symbol:str = str(data[0])
        self.number:int = int(data[1])
        self.mass:float = float(data[2])
        self.charge:Charge = Charge(data[3:6])
        self.structure:str = str(data[6])

    def __eq__(self, __value:Atom) -> bool:
        return (
            self.quantity == __value.quantity and
            self.symbol == __value.symbol and
            self.number == __value.number and
            self.mass == __value.mass and
            self.charge == __value.charge and
            self.structure == __value.structure
        )
    
    def __str__(self) -> str:
        symbol = f'({self.symbol})' if self.ispoly() else self.symbol
        quantity = str(self.quantity) if self.quantity>=2 else ''
        reuslt = 'null' if self.isnull() else symbol+quantity
        return reuslt
    
    def __repr__(self) -> str:
        return str([self.quantity, self.symbol, self.number, self.mass, self.charge, self.structure])
    

    @staticmethod
    def null() -> Atom:
        return Atom(['0'] * 7)


    @property
    def count(self) -> dict:
        counter:dict[str:int] = dict()
        if self.ispoly():
            for each in parse(self.symbol)[1]:
                counter[each.symbol] = each.quantity * self.quantity
        else:
            counter[self.symbol] = self.quantity
        return counter


    def ispoly(self) -> bool:
        return self.structure=='0'

    def isnull(self) -> bool:
        return self.symbol=='0'
    


class Rule:
    def __init__(self, data:list) -> None:
        if data[0] == 'anion':
            self.__negative = Table.anion
        else:
            self.__negative = set(data[0].split(','))

        if '~' in data[1]:
            temp = data[1][1:].split(',')
            self.__positive = {i for i in Table.cation if i not in temp}
        elif data[1] == 'cation':
            self.__positive = Table.cation
        else:
            self.__positive = set(data[1].split(','))


    def match(self, a, b) -> bool:
        return (
            (a in self.__positive and b in self.__negative) or
            (b in self.__positive and a in self.__negative)
        )



class Path:
    current = os.path.dirname(os.path.abspath(__file__))
    periodic = os.path.join(current, 'data', 'periodic.csv')
    others = os.path.join(current, 'data', 'others.txt')
    solubility = os.path.join(current, 'data', 'solubility.txt')



class Table:
    periodic:   dict[str:Atom]  = dict()
    activity:   list[str]       = list()
    diatomic:   set[str]        = set()
    cation:     set[str]        = set()
    anion:      set[str]        = set()
    solubility: set[Rule]       = set()


    @classmethod
    def initialize(cls) -> None:
        with open(Path.periodic) as file:
            for each in file:
                line = each.strip('\n').split(',')
                cls.periodic[line[0]] = Atom(line)

        cls.cation = {key for key, value in cls.periodic.items() if value.charge.head>0}
        cls.anion = {key for key, value in cls.periodic.items() if value.charge.head<0}

        with open(Path.solubility) as file:
            for each in file:
                line = each.strip('\n').split(':')
                cls.solubility.add(Rule(line))

        with open(Path.others) as file:
            cls.activity = file.readline().strip('\n').split(',')
            cls.diatomic = set(file.readline().strip('\n').split(','))
            