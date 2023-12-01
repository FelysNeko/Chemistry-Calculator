from copy import deepcopy
import re

def lookup(molecule:str) -> object:
    if molecule in Table.periodic:
        return deepcopy(Table.periodic[molecule])
    else:
        return Atom.null('0')



def parse(molecule:str) -> tuple:
    poly = re.findall(r'\((.+?)\)(\d*)', molecule)
    molecule = re.sub(r'\((.+?)\)(\d*)', '', molecule)
    atom = re.findall(r'([A-Z][a-z]{0,1})(\d*)', molecule)
    coef = re.findall(r'(^\d*)', molecule)
    
    cvt = lambda x:1 if x=='' else int(x)
    def shrink(name:str, coef:int) -> Atom:
        atom = lookup(name)
        atom.quantity = cvt(coef)
        return atom
    
    coef = cvt(coef[0])
    data = [shrink(i[0], i[1]) for i in atom+poly if i[0] in Table.periodic]
    return coef, data


class Path:
    periodic = './api/periodic.csv'
    others = './api/others.txt'
    solubility = './api/solubility.txt'



class Table:
    periodic = {}
    activity = []
    diatomic = []
    cation = []
    anion = []
    solubility = []
    polyatomic = []



class Charge:
    def __init__(self, data:list) -> None:
        self.data = [int(i) for i in data]


    def __eq__(self, __value:object) -> bool:
        return not sum([self.data[i]-__value.data[i] for i in range(3)])


    @property
    def head(self):
        return self.data[0]


    def match(self, charge:int) -> int:
        try:
            index = [i+charge for i in self.data].index(0)
            self.data[0], self.data[index] = self.data[index], self.data[0]
        except Exception:
            return 0
        else:
            return self.head



class Atom:
    def __init__(self, data:list) -> None:
        self.quantity = 0
        self.symbol = str(data[0])
        self.number = int(data[1])
        self.mass = float(data[2])
        self.charge = Charge(data[3:6])
        self.structure = str(data[6])


    def __eq__(self, __value:object) -> bool:
        return (
            self.quantity == __value.quantity and
            self.symbol == __value.symbol and
            self.number == __value.number and
            self.mass == __value.mass and
            self.charge == __value.charge and
            self.structure == __value.structure
        )


    @property
    def info(self):
        return [self.quantity, self.symbol, self.number, self.mass, self.charge.data, self.structure]


    @property
    def count(self):
        counter = {}
        if self.structure == '0':
            for each in parse(self.symbol)[1]:
                counter[each.symbol] = each.quantity * self.quantity
        else:
            counter[self.symbol] = self.quantity
        return counter
            

    @staticmethod
    def null(value:str) -> object:
        return Atom([value]*7)
    

    @staticmethod
    def isna(x:object) -> bool:
        return True if isinstance(x, Atom) and x.symbol.isdigit() else False



class Rule:
    def __init__(self, data:list) -> None:
        if data[0] == 'anion':
            self.negative = Table.anion
        else:
            self.negative = data[0].split(',')

        if '~' in data[1]:
            temp = data[1][1:].split(',')
            self.positive = [i for i in Table.cation if i not in temp]
        elif data[1] == 'cation':
            self.positive = Table.cation
        else:
            self.positive = data[1].split(',')



with open(Path.others) as file:
    lines = file.readlines()
    Table.activity = lines[0].strip('\n').split(',')
    Table.diatomic = lines[1].strip('\n').split(',')

with open(Path.periodic) as file:
    for each in file.readlines():
        line = each.strip('\n').split(',')
        Table.periodic[line[0]] = Atom(line)


Table.cation = [key for key, value in Table.periodic.items() if value.charge.head>0]
Table.anion = [key for key, value in Table.periodic.items() if value.charge.head<0]
Table.polyatomic = [key for key, value in Table.periodic.items() if value.structure=='0']

with open(Path.solubility) as file:
    for each in file.readlines():
        line = each.strip('\n').split(':')
        Table.solubility.append(Rule(line))