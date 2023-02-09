# PROGRAM IS UNDER RECONSTRUCTION
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
            information = lookup(molecule[begin + 1:end])
        else:
            for each in molecule:
                if each.isupper() is True:
                    break
            begin = molecule.find(each)
            end = begin+1 if molecule[-1] != each and molecule[begin + 1].islower() is True else begin
            information = lookup(molecule[begin:end + 1])

        if information is not None:
            quantity = 1
            for i in range(end + 1, len(molecule) + 1):
                if i == len(molecule) or molecule[i].isdigit() is False:
                    break
            if molecule[end + 1:i] != '':
                quantity = int(molecule[end + 1:i])
        else:
            return None

        molecule = molecule[:begin] + molecule[i:]
        result.insert(0, [quantity, information])

    result.insert(0, int(molecule))
    return result
