import random
import sqlite3
import sys
import time
import math


OXID = 'oxidizer'
RECO = 'reductant'
specs = ('id', 'name', 'name_latin', 'short_name', 'mass', 'radius', 'column', 'row', 'type', 'electronegativity',
         'density', 'cas_number', 'specs')


class Element:
    def __init__(self, element: str | int, oxidation: int = None):
        self.el = element
        self.info = self.get_info()
        self.neutrons = round(self.info['mass'], 0) - self.info['id']
        self.protons = self.info['id']
        self.electrons = self.protons

        if self.info['type'] == 'n':
            self.type = OXID
        else:
            self.type = RECO

        if oxidation is None:
            self.last_electrons = Chemistry().get_oxidation_degree(self)
        else:
            self.last_electrons = oxidation

    def __ne__(self, other):
        if not isinstance(self, type(other)):
            return True
        if other.el.lower() != self.el.lower():
            return True
        if other.info != self.info:
            return True
        if other.neutrons != self.neutrons:
            return True
        if other.protons != self.protons:
            return True
        if other.electrons != self.electrons:
            return True
        if other.type != self.type:
            return True
        return False

    def change_type(self, set_type: str = None):
        if set_type is not None:
            if self.type != set_type:
                self.type = set_type
        else:
            if self.type == RECO:
                self.type = OXID
            else:
                self.type = RECO

        self.last_electrons = Chemistry().get_oxidation_degree(self)

        return self

    def get_info(self):
        with sqlite3.connect("base.db") as sql:
            try:
                return dict(zip(specs, tuple(sql.cursor().execute(f"SELECT * FROM elements WHERE short_name = ?",
                                                                  (self.el.lower().capitalize(),)))[0]))
            except IndexError:
                try:
                    return dict(zip(specs, tuple(sql.cursor().execute(f"SELECT * FROM elements WHERE name = ?",
                                                                      (self.el.lower().capitalize(),)))[0]))
                except IndexError:
                    try:
                        return dict(zip(specs, tuple(sql.cursor().execute(f"SELECT * FROM elements WHERE name_latin = "
                                                                          f"?", (self.el.lower().capitalize(),)))[0]))
                    except IndexError:
                        try:
                            return dict(zip(specs, tuple(sql.cursor().execute(f"SELECT * FROM elements WHERE id = ?",
                                                                              (self.el,)))[0]))
                        except IndexError:
                            print(f"Элемент не найден - {self.el}")
                            time.sleep(5)
                            sys.exit()

    def __repr__(self):
        return self.info['short_name']


class Group:
    def __init__(self, *elements, ion: bool = False):
        self.elements = elements

        if ion:
            self.last_electrons = elements[0][1] * elements[0][0].last_electrons - elements[1][1] * \
                                  elements[1][0].last_electrons
        else:
            self.last_electrons = elements[1][1] * elements[1][0].last_electrons - elements[0][1] * \
                                  elements[0][0].last_electrons

    def get_name(self, brackets: bool = True):
        strr = ''
        for i in self.elements:
            if type(i[0]) is Group:
                strr = i[0].get_name()
                if i[1] != 1:
                    strr += str(i[1])
            elif type(i[0]) is Element:
                strr += i[0].info['short_name']
                if i[1] != 1:
                    strr += str(i[1])
        if brackets:
            return '(' + strr + ')'
        else:
            return strr

    def __repr__(self):
        return ''.join(f"{i[0]}{i[1] if i[1] != 1 else ''}" for i in self.elements)


class Substance:
    def __init__(self, *elss: Element):
        self.els = self.construct(*elss)
        self.raw = self.get_raw()

    def __repr__(self):
        return self.raw

    def __ne__(self, other):
        if not isinstance(self, type(other)):
            return True
        for el1, el2 in zip(other.els, self.els):
            if el1[1] != el2[1]:
                return True
            if el1[0] != el2[0]:
                return True
        return False

    @staticmethod
    def construct_group(element_1: Element, element_2: Element):
        if element_1.info['short_name'] == 'O' and element_2.info['short_name'] == 'H':
            return Group((element_1, 1), (element_2, 1), ion=True)
        else:
            element_1.change_type()

            el_1_ox, el_2_ox = 1, 1

            while el_2_ox * element_2.last_electrons <= element_1.last_electrons:
                el_2_ox += 1

            grpp = Group((element_1, el_1_ox), (element_2, el_2_ox))

            return grpp

    def construct(self, *els: Element | Group):
        els = tuple(set(els))
        metals = set()

        if len(els) == 1:
            if els[0].info['type'] == 'n':
                return (els[0], 2),

        for i in els:
            if i.info['type'] == 'm':
                metals.add(i)
                del i

        els = set(els) - metals
        nonmetal_sorted = sorted(els, key=lambda x: x.info['electronegativity'])
        el_sorted = list((*sorted(metals, key=lambda x: x.info['electronegativity']), *nonmetal_sorted))

        if len(el_sorted) > 2:
            for i in range(len(el_sorted)):
                if type(el_sorted[i]) is Element:
                    try:
                        if el_sorted[i].info['short_name'] == "H" and el_sorted[i + 1].info['short_name'] == "O":
                            el_sorted[i], el_sorted[i + 1] = el_sorted[i + 1], el_sorted[i]
                    except IndexError:
                        pass

        while len(el_sorted) > 2:
            constt = self.construct_group(el_sorted[-2], el_sorted[-1])
            el_sorted = el_sorted[0:-2]
            el_sorted.append(constt)

        oxids = tuple(i.last_electrons for i in el_sorted)

        zipped = tuple(zip((i for i in el_sorted), (int(math.lcm(*oxids) / i) for i in oxids)))

        return zipped

    def get_raw(self):
        itg = ""
        for i in self.els:
            if type(i[0]) is Element:
                itg += i[0].info['short_name']
                if i[1] != 1:
                    itg += str(i[1])
            elif type(i[0]) is Group:
                if i[1] != 1:
                    itg += i[0].get_name()
                    itg += str(i[1])
                else:
                    itg += i[0].get_name(brackets=False)
            else:
                print(f"Класс элемента не найден - {i}")
                time.sleep(5)
                sys.exit()
        return itg


class Chemistry:
    def __init__(self):
        ...

    @staticmethod
    def get_oxidation_degree(elem: Element):
        if elem.info['short_name'] == 'H':
            return 1
        if elem.type == OXID:
            return 8 - elem.info['column']
        elif elem.type == RECO:
            return elem.info['column']
        else:
            print(f"Тип элемента не найден - {elem.info['name']}")
            time.sleep(5)
            sys.exit()


class Require:
    def __init__(self, **requirements):
        self.reqs = requirements

    def get_requirements(self):
        return self.reqs

    def check(self, to_check: tuple[Element | Substance]):
        if len(to_check) < self.reqs['id']:
            return False
        cur_el = to_check[self.reqs['id'] - 1]

        for name, val in self.reqs.items():
            if name == 'needed':
                if not isinstance(cur_el, type(val)):
                    return False
                if cur_el != val:
                    return False
            elif name == 'specs':
                if not isinstance(cur_el, Element):
                    return False
                for val_f in val:
                    if val_f not in cur_el.info['specs'].split(';'):
                        return False
        return True


reactions = {
    'neutralization': {
        'reagents': {
            (
                Require(specs=('active_metal', ), id=1),
                Require(needed=Substance(Element('H'), Element('O')), id=2)
            )
        },
        'products': {
            (
                Substance(Element('Ca'), Element('O'), Element('H')),
                Substance(Element('H'))
            )
        }
    }
}


class Reaction:
    def __init__(self, *reagents: Substance | Element):
        self.reagents = reagents

    def process(self):
        reacts_avails = []

        for reac_name, reac_info in reactions.items():
            for reqs in reac_info['reagents']:
                if self.check_valid(reqs):
                    reacts_avails.append(reac_name)
        cur_reac = reactions[random.choice(reacts_avails)]

        return cur_reac['products']

    def check_valid(self, reqs: tuple[Require, Require]):
        for req in reqs:
            if not req.check(self.reagents):
                return False
        return True


Li = Element('Li')
H2O = Substance(Element("h"), Element("o"))

print(Reaction(Li, H2O).process())
input()
