# ============================================================
# DNA Analyzer — BioSeq Lab
# Point d'entrée principal
# Lancer avec : python main.py
# ============================================================

import tkinter as tk
from gui.app import DNAAnalyzerApp

if __name__ == '__main__':
    root = tk.Tk()
    app = DNAAnalyzerApp(root)
    root.mainloop()
