from __future__ import annotations
from .basis import Table, parse
from copy import deepcopy
import sympy as sp


class Molecule:
    def __init__(self, molecule:str) -> None:
        coef, data = parse(molecule)
        self.coef = coef
        self.data = data

    def __eq__(self, __value:Molecule) -> bool:
        return (
            self.coef == __value.coef and
            self.data == __value.data
        )
    
    def __str__(self) -> str:
        coef = '' if self.coef==1 else str(self.coef)
        molecule = ''.join([str(i) for i in self.data])
        return 'null' if self.isnull() else coef+molecule 

    def __repr__(self) -> str:
        return str(self.coef) + ' ' + str(self.data)
    

    @staticmethod
    def null() -> Molecule:
        return Molecule('')


    @property
    def mass(self) -> float:
        return round(self.coef * sum([i.mass * i.quantity for i in self.data]), 4)
    
    @property
    def bond(self) -> str:
        if self.isnull():
            return 'unknown'

        if len(self.data)==2 and self.data[0].charge.head*self.data[1].charge.head < 0:
            return 'ionic'
        elif sum([i.charge.head<0 for i in self.data]) == len(self.data):
            return 'covalent'
        elif sum([i.charge.head>0 for i in self.data]) == len(self.data):
            return 'metallic'
        else:
            return 'unknown'


    def isnull(self) -> bool:
        return not len(self.data)
    
