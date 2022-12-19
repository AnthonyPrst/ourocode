# coding in UTF-8 
# by Anthony PARISOT

import pickle
from tkinter import filedialog, messagebox


class Objet(object):
    """Classe permetant la sauvegarde ou l'ouverture d'un objet ou de plusieur sous un fichier .ec
    """
    def __init__(self) -> None:
        pass
    
    def save_muliple_objects(self, object: list):
        with filedialog.asksaveasfile('wb', filetypes=(("savefile EC object", "*.ECobjt"), ('Text Document', '*.txt')), defaultextension='.ECobjt') as f:
            for ligne in object:
                pickle.dump(ligne, f)
    
    def save_object(self):
        with filedialog.asksaveasfile('wb', filetypes=(("savefile EC object", "*.ECobjt"), ('Text Document', '*.txt')), defaultextension='.ECobjt') as f:
            pickle.dump(self, f)
            
            
    @classmethod
    def open_multiple_objects(cls):
        data = []
        with filedialog.askopenfile('rb', filetypes=(("savefile EC object", "*.ECobjt"), ('Text Document', '*.txt')), defaultextension='.ECobjt') as f:
            while True:
                try:
                    data.append(pickle.load(f))
                except EOFError:
                    break
            return data
    
    @classmethod
    def open_object(cls):
        with filedialog.askopenfile('rb', filetypes=(("savefile EC object", "*.ECobjt"), ('Text Document', '*.txt')), defaultextension='.ECobjt') as f:
            return pickle.load(f)