# coding in UTF-8 
# by Anthony PARISOT

import pickle
from tkinter import filedialog


class Objet(object):
    """Classe permetant la sauvegarde ou l'ouverture d'un objet ou de plusieur sous un fichier .ec
    """
    def __init__(self) -> None:
        pass


    @classmethod
    def _from_parent_class(cls, object, **kwargs):
        """Class méthode permetant l'intanciation des classe hérité de la classe parent, par une classe déjà instanciée.

        Args:
            object (class object): l'objet Element déjà créer par l'utilisateur
        """ 
        return cls(**object.__dict__, **kwargs)

    
    def _save_muliple_objects(self, object: list):
        with filedialog.asksaveasfile('wb', filetypes=(("savefile EC object", "*.ECobjt"), ('Text Document', '*.txt')), defaultextension='.ECobjt') as f:
            for ligne in object:
                pickle.dump(ligne, f)
    
    def _save_object(self):
        with filedialog.asksaveasfile('wb', filetypes=(("savefile EC object", "*.ECobjt"), ('Text Document', '*.txt')), defaultextension='.ECobjt') as f:
            pickle.dump(self, f)
            
            
    @classmethod
    def _open_multiple_objects(cls):
        data = []
        with filedialog.askopenfile('rb', filetypes=(("savefile EC object", "*.ECobjt"), ('Text Document', '*.txt')), defaultextension='.ECobjt') as f:
            while True:
                try:
                    data.append(pickle.load(f))
                except EOFError:
                    break
            return data
    
    @classmethod
    def _open_object(cls):
        with filedialog.askopenfile('rb', filetypes=(("savefile EC object", "*.ECobjt"), ('Text Document', '*.txt')), defaultextension='.ECobjt') as f:
            return pickle.load(f)