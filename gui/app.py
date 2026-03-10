# ============================================================
# Fenêtre principale — Thème clair biologique
# ============================================================

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import math

from gui.input_frame     import InputFrame
from gui.results_frame   import ResultsFrame
from gui.visualize_frame import VisualizeFrame

from analysis.orf_finder   import find_orfs
from analysis.motif_finder import (find_promoters, find_shine_dalgarno,
                                   find_terminators, find_restriction_sites)
from analysis.statistics   import (nucleotide_composition, calc_tm_wallace,
                                   calc_tm_nearest_neighbor, molecular_weight)

from data.codon_table import RESTRICTION_SITES
from export.exporter  import (export_fasta, export_csv, export_json,
                               export_txt_report, export_excel)

# ── Palette biologique claire ────────────────────────────────
BG_MAIN    = '#F0F9F4'
BG_SURFACE = '#FFFFFF'
BG_PANEL   = '#E8F5EE'
BORDER     = '#B8DFC9'
FG_MAIN    = '#1A3A2A'
FG_DIM     = '#6B9E80'
FG_LABEL   = '#2D6A4F'
ACCENT     = '#1B7A4A'
ACCENT2    = '#2196A6'
ACCENT3    = '#E67E22'
ACCENT4    = '#8E44AD'


class DNAAnalyzerApp:
    def __init__(self, root):
        self.root = root
        self.root.title('🧬 DNA Analyzer — BioSeq Lab')
        self.root.geometry('1350x860')
        self.root.minsize(950, 650)
        self.root.configure(bg=BG_MAIN)
        self.results = {}
        self._apply_styles()
        self._build_menu()
        self._build_layout()

    def _apply_styles(self):
        style = ttk.Style()
        style.theme_use('clam')
        style.configure('TNotebook', background=BG_MAIN, borderwidth=0)
        style.configure('TNotebook.Tab', background=BG_PANEL, foreground=FG_DIM,
                         font=('Georgia', 9), padding=[12, 6], borderwidth=1)
        style.map('TNotebook.Tab',
                  background=[('selected', BG_SURFACE)],
                  foreground=[('selected', ACCENT)])
        style.configure('Treeview', background=BG_SURFACE, foreground=FG_MAIN,
                         fieldbackground=BG_SURFACE, rowheight=28,
                         font=('Courier New', 9), borderwidth=0)
        style.configure('Treeview.Heading', background=BG_PANEL, foreground=FG_LABEL,
                         font=('Georgia', 9, 'bold'), borderwidth=0, relief='flat')
        style.map('Treeview', background=[('selected', '#C8E6C9')],
                  foreground=[('selected', FG_MAIN)])
        style.configure('Vertical.TScrollbar', background=BG_PANEL,
                         troughcolor=BG_MAIN, arrowcolor=FG_DIM, borderwidth=0)
        style.configure('Horizontal.TScrollbar', background=BG_PANEL,
                         troughcolor=BG_MAIN, arrowcolor=FG_DIM, borderwidth=0)

    def _build_menu(self):
        menubar = tk.Menu(self.root, bg=BG_SURFACE, fg=FG_MAIN,
                          activebackground='#C8E6C9', activeforeground=FG_MAIN, tearoff=0)
        file_menu = tk.Menu(menubar, tearoff=0, bg=BG_SURFACE, fg=FG_MAIN,
                             activebackground='#C8E6C9')
        file_menu.add_command(label='📂  Ouvrir FASTA…',        command=self._open_file)
        file_menu.add_separator()
        file_menu.add_command(label='💾  Exporter CSV…',         command=lambda: self._export('csv'))
        file_menu.add_command(label='📊  Exporter Excel…',       command=lambda: self._export('excel'))
        file_menu.add_command(label='🗄️   Exporter JSON…',        command=lambda: self._export('json'))
        file_menu.add_command(label='📝  Exporter Rapport TXT…', command=lambda: self._export('txt'))
        file_menu.add_command(label='🧬  Exporter FASTA…',        command=lambda: self._export('fasta'))
        file_menu.add_separator()
        file_menu.add_command(label='Quitter', command=self.root.quit)
        menubar.add_cascade(label='Fichier', menu=file_menu)
        analysis_menu = tk.Menu(menubar, tearoff=0, bg=BG_SURFACE, fg=FG_MAIN,
                                 activebackground='#C8E6C9')
        analysis_menu.add_command(label='▶  Lancer toutes les analyses', command=self._run_analysis)
        analysis_menu.add_separator()
        analysis_menu.add_command(label='↔  Brin complémentaire inverse', command=self._reverse_complement)
        menubar.add_cascade(label='Analyse', menu=analysis_menu)
        help_menu = tk.Menu(menubar, tearoff=0, bg=BG_SURFACE, fg=FG_MAIN,
                             activebackground='#C8E6C9')
        help_menu.add_command(label='À propos', command=self._show_about)
        menubar.add_cascade(label='Aide', menu=help_menu)
        self.root.config(menu=menubar)

    def _build_layout(self):
        self._build_title_bar()
        main = tk.Frame(self.root, bg=BG_MAIN)
        main.pack(fill='both', expand=True, padx=8, pady=(0, 8))
        left = tk.Frame(main, bg=BG_PANEL, width=330,
                        highlightthickness=1, highlightbackground=BORDER)
        left.pack(side='left', fill='y', padx=(0, 6))
        left.pack_propagate(False)
        self.input_frame = InputFrame(left, self._run_analysis)
        self.input_frame.pack(fill='both', expand=True)
        right = tk.Frame(main, bg=BG_MAIN)
        right.pack(side='right', fill='both', expand=True)
        self.notebook = ttk.Notebook(right)
        self.notebook.pack(fill='both', expand=True)
        self.results_frame   = ResultsFrame(self.notebook)
        self.visualize_frame = VisualizeFrame(self.notebook)
        self.notebook.add(self.results_frame,   text='  🔬  Résultats  ')
        self.notebook.add(self.visualize_frame, text='  📊  Visualisation  ')

    def _build_title_bar(self):
        bar = tk.Frame(self.root, bg=ACCENT, height=58)
        bar.pack(fill='x', side='top')
        bar.pack_propagate(False)
        self.dna_canvas = tk.Canvas(bar, bg=ACCENT, height=58,
                                     bd=0, highlightthickness=0)
        self.dna_canvas.pack(fill='both', expand=True)
        tk.Label(self.dna_canvas, text='🧬  DNA Analyzer — BioSeq Lab',
                 bg=ACCENT, fg='#FFFFFF',
                 font=('Georgia', 14, 'bold')).place(x=16, y=13)
        self.status_label = tk.Label(self.dna_canvas, text='Prêt',
                                     bg=ACCENT, fg='#A8D5B8',
                                     font=('Georgia', 9, 'italic'))
        self.status_label.place(relx=1.0, x=-16, y=19, anchor='ne')
        self._dna_offset = 0.0
        self._animate_dna()

    def _animate_dna(self):
        c = self.dna_canvas
        c.delete('helix')
        try:
            w = c.winfo_width()
            h = c.winfo_height()
        except Exception:
            self.root.after(80, self._animate_dna)
            return
        if w < 10:
            self.root.after(80, self._animate_dna)
            return
        step = 16
        amp  = 9
        freq = 0.048
        off  = self._dna_offset
        pts1 = [(x, h/2 + amp * math.sin(freq*x + off))        for x in range(0, w+step, step)]
        pts2 = [(x, h/2 + amp * math.sin(freq*x + off + math.pi)) for x in range(0, w+step, step)]
        for i in range(len(pts1)-1):
            c.create_line(pts1[i][0], pts1[i][1], pts1[i+1][0], pts1[i+1][1],
                           fill='#A8D5B8', width=2, smooth=True, tags='helix')
        for i in range(len(pts2)-1):
            c.create_line(pts2[i][0], pts2[i][1], pts2[i+1][0], pts2[i+1][1],
                           fill='#C8E6D0', width=2, smooth=True, tags='helix')
        rung_cols = ['#74C69D', '#52B788', '#40916C']
        for i, (p1, p2) in enumerate(zip(pts1[::2], pts2[::2])):
            col = rung_cols[i % len(rung_cols)]
            c.create_line(p1[0], p1[1], p2[0], p2[1], fill=col, width=1.5, tags='helix')
            r = 3
            c.create_oval(p1[0]-r, p1[1]-r, p1[0]+r, p1[1]+r,
                           fill='#D8F3DC', outline='', tags='helix')
            c.create_oval(p2[0]-r, p2[1]-r, p2[0]+r, p2[1]+r,
                           fill='#B7E4C7', outline='', tags='helix')
        self._dna_offset -= 0.04
        self.root.after(40, self._animate_dna)

    def _run_analysis(self, seq=None, options=None):
        if seq is None:
            seq = self.input_frame.get_sequence()
        if not seq or len(seq) < 20:
            messagebox.showwarning('Séquence invalide',
                                   'Minimum 20 nucléotides requis.')
            return
        if options is None:
            options = self.input_frame.get_options()
        self._set_status('⏳  Analyse en cours…')
        self.root.update_idletasks()
        try:
            results = {}
            seq = seq.upper().replace(' ', '').replace('\n', '')
            min_orf = options.get('min_orf', 90)
            results['orfs'] = find_orfs(seq, min_length=min_orf) if options.get('ORFs', True) else []
            pd = find_promoters(seq) if options.get('Promoteurs', True) else {'promoters':[],'box10':[],'box35':[]}
            results['promoters'] = pd['promoters']
            results['box10']     = pd['box10']
            results['box35']     = pd['box35']
            results['sd_sites']    = find_shine_dalgarno(seq) if options.get('Shine-Dalgarno', True) else []
            results['terminators'] = find_terminators(seq)    if options.get('Terminateurs', True)   else []
            results['restriction'] = find_restriction_sites(seq, RESTRICTION_SITES) if options.get('Restriction', True) else {}
            stats = nucleotide_composition(seq)
            stats['tm_wallace'] = calc_tm_wallace(seq)
            stats['tm_nn']      = calc_tm_nearest_neighbor(seq)
            stats['mw_ss']      = molecular_weight(seq, double_stranded=False)
            stats['mw_ds']      = molecular_weight(seq, double_stranded=True)
            results['stats']      = stats
            results['seq']        = seq
            results['seq_length'] = len(seq)
            results['gc_pct']     = stats['gc']
            self.results = results
            self.results_frame.display(results)
            self.visualize_frame.draw(seq, results)
            total = (len(results['orfs']) + len(results['promoters']) +
                     len(results['sd_sites']) + len(results['terminators']))
            self._set_status(f'✓  Analyse terminée — {total} éléments détectés')
        except Exception as e:
            messagebox.showerror('Erreur', str(e))
            self._set_status('❌  Erreur')

    def _set_status(self, msg):
        self.status_label.config(text=msg)
        self.root.update_idletasks()

    def _open_file(self):
        path = filedialog.askopenfilename(filetypes=[('FASTA', '*.fasta *.fa *.txt *.seq'), ('Tous', '*.*')])
        if path:
            with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            lines = [l.strip() for l in content.splitlines() if not l.startswith('>')]
            self.input_frame.set_sequence(''.join(lines).upper())

    def _reverse_complement(self):
        from analysis.orf_finder import reverse_complement
        seq = self.input_frame.get_sequence()
        if seq:
            self.input_frame.set_sequence(reverse_complement(seq))

    def _export(self, fmt):
        if not self.results:
            messagebox.showinfo('Aucun résultat', 'Lancez d\'abord une analyse.')
            return
        exts = {'csv':'.csv','excel':'.xlsx','json':'.json','txt':'.txt','fasta':'.fasta'}
        path = filedialog.asksaveasfilename(defaultextension=exts.get(fmt,'.txt'))
        if not path:
            return
        try:
            seq = self.results.get('seq', '')
            if   fmt == 'csv':   export_csv(self.results, path)
            elif fmt == 'excel': export_excel(self.results, seq, path)
            elif fmt == 'json':  export_json(self.results, path)
            elif fmt == 'txt':   export_txt_report(self.results, seq, path)
            elif fmt == 'fasta': export_fasta(seq, path)
            messagebox.showinfo('Export réussi', f'Fichier enregistré :\n{path}')
        except Exception as e:
            messagebox.showerror('Erreur d\'export', str(e))

    def _show_about(self):
        messagebox.showinfo('À propos',
            'DNA Analyzer — BioSeq Lab\n\nPython / Tkinter\nThème biologique clair')