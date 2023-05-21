class Atom:
    def __init__(self, quantity, data) -> None:
        self.qtt = quantity
        self.sym = data[0]
        self.num = int(data[1])
        self.mass = float(data[2])
        self.chrg = [int(i) for i in data[3:-1] if i != '0']
        self.structure = data[-1] if data[-1] != '0' else None


    @property
    def info(self):
        return [self.qtt, self.sym, self.num, self.mass, self.chrg, self.structure]


class Rule:
    def __init__(self, para) -> None:
        if para[0] == 'anion':
            self.neg = Table.anion
        else:
            self.neg = para[0].split(',')

        if '~' in para[1]: 
            self.pos = [i for i in Table.cation if i not in para[1][1:].split(',')]
        elif para[1] == 'cation':
            self.pos = Table.cation
        else:
            self.pos = para[1].split(',')


class Table:
    periodic = {}
    diatomic = []
    activity = []
    cation = []
    anion = []
    solubility = []


with open('api/library/others.csv') as data:
    lines = data.readlines()
    Table.activity = lines[0].strip('\n').split(',')
    Table.diatomic = lines[1].strip('\n').split(',')

with open('api/library/periodic.csv') as data:
        for line in data.readlines():
            pre = line.strip('\n').split(',')
            Table.periodic[pre[0]] = Atom(None, pre)

Table.cation = [key for key, value in Table.periodic.items() if len(value.chrg) != 0 and value.chrg[0] > 0]
Table.anion = [key for key, value in Table.periodic.items() if len(value.chrg) != 0 and value.chrg[0] < 0]

with open('api/library/solubility.txt') as data:
    for each in data.readlines():
        line = each.strip('\n').split(':')
        Table.solubility.append(Rule(line))
