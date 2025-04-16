from ourocode.eurocode.objet import Objet
    
class And(Objet):
    def __init__(self):
        super().__init__()
        pass


class For(Objet):
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
        self.selected_item = None

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
        return self.items

    def _get_item(self, selection="Tous"):
        """Retourne un élément spécifique ou la sélection."""
        if selection == "Tous":
            self.counter += 1
            if self.counter >= len(self.items):
                self.counter = 0
            self.selected_item = self.items[self.counter]
        else:
            for item in self.items:
                if item == selection:
                    self.selected_item = item
        return self.selected_item
    
    def get_min_value_in_loop():
        """Retourne la valeur de sortie de la boucle pour laquelle le résultat (valeur numérique uniquement) est le plus faible.
        Le résultat peut ce trouver dans un dictionnaire, une liste, un tuple. Si c'est le cas alors on parcours toute les valeurs numérique et l'on récupère la plus faible.
        Cette méthode permet pour exemple de trouver le taux de travail le plus faible dans une 

        Args:
            result (float): le résultat à comparer avec la précédente valeur minimal

        """
                
