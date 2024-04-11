import api.base.constants as cts
from catalog.eurocode.objet import Objet

    
class And(Objet):
    def __init__(self):
        super().__init__()
        pass


class For(Objet):
    SELECTION = ["Tous"]
    def __init__(self, liste: list, index: int=None, key: str=None, start: int=None, stop: int=None, **kwargs):
        """Boucle sur les valeurs d'une liste, en option on peux donner un index pour retourner une valeur spécifique d'une sous liste et une clé pour un dictionnaire.
        On peut aussi donner un start et un stop pour boucler que sur une partie de la liste.

        Args:
            liste (list): la liste sur laquelle on doit itérer.
            index (int|otional): l'index optionnel.
            key (str, optional): on utilise la clé spécifique si la liste transmise est un dictionnaire.

        Returns:
            Retourne l'item présent dans la liste.
        """
        super().__init__()
        self.liste = liste
        self.index = index
        self.key = key
        self.start = start
        self.stop = stop
        self.counter = -1

        for key, val in kwargs.items():
            setattr(self, key, val)
                        
        self._get_items()
        
    def _get_items(self):
        """Détermine le début d'une boucle
        """
        self.items = []
        if isinstance(self.liste, list):
            strt = 0
            stp = len(self.liste)+1
            if self.start:
                strt = self.start
            if self.stop:
                stp = self.stop
            self.liste = self.liste[strt:stp]
            for item in self.liste:
                if self.index:
                    self.items.append(item[self.index])
                else:
                    self.items.append(item)
        elif isinstance(self.liste, dict):
            for item in self.liste.values():
                if self.key:
                    self.items.append(item[self.key])
                else:
                    self.items.append(item)
        self.set_selection(*self.items)
        return self.items
    
    def set_selection(cls, *items):
        for item in items:
            cls.SELECTION.append(item)
        print(cls.SELECTION)

    def get_item(self, selection: str=SELECTION):
        if selection == "Tous":
            self.counter += 1
            if self.counter > len(self.items)-1:
                self.counter = 0
                return False
            print(self.counter, self.items[self.counter])
            return self.items[self.counter]
        else:
            selection
    
    def end_for(self):
        return f"j'arrète la boucle for, cela fait {self.counter+1} fois que je boucle."
                
