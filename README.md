# 🧬 DNA Analyzer — BioSeq Lab

> Application Python d'analyse bioinformatique de séquences ADN bactériennes.  
> Interface graphique complète, 6 types d'analyses, visualisations interactives et export multi-format.

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?style=flat-square&logo=python&logoColor=white)
![Tkinter](https://img.shields.io/badge/GUI-Tkinter-1B7A4A?style=flat-square)
![Matplotlib](https://img.shields.io/badge/Plots-Matplotlib-orange?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)

---

## 📋 Table des matières

- [Aperçu](#-aperçu)
- [Fonctionnalités](#-fonctionnalités)
- [Structure du projet](#-structure-du-projet)
- [Installation](#-installation)
- [Lancement](#-lancement)
- [Guide d'utilisation](#-guide-dutilisation)
- [Description des analyses](#-description-des-analyses)
- [Export des résultats](#-export-des-résultats)
- [Architecture technique](#-architecture-technique)
- [Dépendances](#-dépendances)

---

## 🔭 Aperçu

DNA Analyzer est une application de bureau développée en **Python / Tkinter** permettant d'analyser des séquences ADN bactériennes à partir d'une saisie directe ou d'un fichier FASTA.

Elle couvre les analyses bioinformatiques fondamentales demandées dans le cadre du projet :

- **Détection des ORFs** sur les 6 cadres de lecture (brin sens et antisens)
- **Identification de l'ORF codant** le plus probable avec traduction protéique
- **Recherche de motifs régulateurs** : promoteurs (-10 / -35), terminateurs Rho-indépendants, sites Shine-Dalgarno
- **Analyses physico-chimiques** : GC%, Tm, masse moléculaire
- **Visualisations graphiques** interactives via Matplotlib
- **Export** CSV, Excel, JSON, TXT, FASTA

---

## ✨ Fonctionnalités

### Analyses biologiques
| Analyse | Description |
|---|---|
| 📋 Cadres de lecture | Lecture codon par codon sur les 6 frames (+1/+2/+3/-1/-2/-3) |
| 🔍 ORFs ATG→Stop | Détection classique sur 6 cadres, classés par longueur |
| ⚙️ Traduction protéique | Traduction de l'ORF principal, colorée par propriété biochimique |
| 📍 Promoteurs | Pairs Box-35 / Box-10 avec tolérance aux mismatches |
| 🎯 Shine-Dalgarno | Consensus AGGAGG localisé 7–9 pb en amont d'ATG |
| 🔚 Terminateurs | Structures tige-boucle palindromiques + polyT |
| ✂️ Sites de restriction | 15 enzymes courantes (EcoRI, BamHI, HindIII…) |
| 📊 Statistiques | GC%, Tm Wallace, Tm Nearest-Neighbor, masse moléculaire |

### Interface & usabilité
- Interface graphique avec **thème biologique** (palette verte)
- Animation ADN double hélice en temps réel dans la barre de titre
- **Analyse threadée** — UI non bloquante sur les longues séquences FASTA
- Barre de statut avec progression étape par étape
- Import de fichiers FASTA (`.fasta`, `.fa`, `.txt`, `.seq`)
- Calcul du **brin complémentaire inverse** en un clic
- Seuil de longueur minimum des ORFs configurable

---

## 📁 Structure du projet

```
dna_analyzer/
│
├── main.py                    # Point d'entrée — lancer avec : python main.py
│
├── analysis/
│   ├── __init__.py
│   ├── orf_finder.py          # Détection ORFs (6 frames) + traduction
│   ├── motif_finder.py        # Promoteurs, SD, terminateurs, restriction
│   └── statistics.py          # Composition, Tm, masse moléculaire, GC glissant
│
├── data/
│   ├── __init__.py
│   └── codon_table.py         # Table des codons + sites de restriction
│
├── gui/
│   ├── __init__.py
│   ├── app.py                 # Fenêtre principale + orchestration analyses
│   ├── input_frame.py         # Panneau gauche : saisie & options
│   ├── results_frame.py       # Panneau droit : onglets de résultats
│   └── visualize_frame.py     # Onglet visualisation Matplotlib
│
├── export/
│   ├── __init__.py
│   └── exporter.py            # Export CSV, Excel, JSON, TXT, FASTA
│
└── requirements.txt           # Dépendances Python
```

---

## 🛠 Installation

### Prérequis

- Python **3.10 ou supérieur**
- `pip` à jour

### 1. Cloner le dépôt

```bash
git clone https://github.com/votre-utilisateur/dna_analyzer.git
cd dna_analyzer
```

### 2. Créer un environnement virtuel (recommandé)

```bash
python -m venv bioenv
source bioenv/bin/activate        # Linux / macOS
bioenv\Scripts\activate           # Windows
```

### 3. Installer les dépendances

```bash
pip install -r requirements.txt
```

#### Dépendances principales

```
matplotlib>=3.7.0
openpyxl>=3.1.0
```

> **Note** : `tkinter` est inclus avec Python standard.  
> Sur Ubuntu/Debian, si tkinter est manquant :
> ```bash
> sudo apt install python3-tk
> ```

---

## 🚀 Lancement

```bash
python main.py
```

L'interface s'ouvre avec une séquence de démonstration pré-chargée, prête à être analysée.

---

## 📖 Guide d'utilisation

### Panneau gauche — Saisie

```
┌─────────────────────────────────┐
│  Séquence ADN (5'→3')           │
│  ┌───────────────────────────┐  │
│  │ ATGCCCGGGTAA...           │  │
│  └───────────────────────────┘  │
│  300 pb  •  GC: 48.3%           │
│                                 │
│  [📁 Importer FASTA] [↔ RC] [✕] │
│                                 │
│  Analyses à effectuer :         │
│  ☑ Détection ORFs (6 cadres)    │
│  ☑ Promoteurs -10 / -35         │
│  ☑ Shine-Dalgarno               │
│  ☑ Terminateurs Rho-indép.      │
│  ☑ Sites de restriction         │
│                                 │
│  Longueur min. ORF (pb) : [30]  │
│                                 │
│  [▶  LANCER L'ANALYSE]          │
└─────────────────────────────────┘
```

| Élément | Rôle |
|---|---|
| Zone de texte | Coller la séquence ADN directement (ATGC uniquement) |
| `📁 Importer FASTA` | Charger un fichier FASTA — les lignes `>` sont ignorées |
| `↔ RC` | Remplace la séquence par son brin complémentaire inverse |
| `✕` | Efface la zone de saisie |
| Cases à cocher | Activer/désactiver chaque analyse individuellement |
| Longueur min. ORF | Seuil en pb — `30` pour tous les ORFs, `90+` pour les candidats biologiques |
| `▶ LANCER L'ANALYSE` | Démarre l'analyse (thread séparé, UI reste réactive) |

---

### Panneau droit — Onglets de résultats

#### 📋 Cadres de lecture
- Affiche les **6 frames** de lecture (+1, +2, +3, -1, -2, -3)
- Chaque frame est colorée distinctement
- Les codons sont affichés un par un : **ATG en vert**, **codons stop en rouge**
- Tableau récapitulatif des segments détectés par frame (début, fin, longueur, nb AA, protéine)

#### 🔍 ORFs (ATG → Stop)
- Liste de tous les ORFs classés par longueur décroissante
- Le meilleur candidat est marqué **★** en vert
- Colonnes : Frame, Début, Fin, Longueur, Nb AA, Brin, Codon Stop

#### ⚙️ Protéine
- Traduction de l'ORF principal (premier ORF du brin sens)
- Séquence colorée par propriété biochimique :
  - 🟠 **Orange** — Acides aminés hydrophobes (A, V, I, L, M, F, W, P)
  - 🔵 **Bleu** — Acides aminés polaires (S, T, N, Q, Y, C)
  - 🟣 **Violet** — Acides aminés chargés (K, R, H, D, E, P)
- Informations : position, longueur ADN, nb AA, poids moléculaire estimé

#### 📍 Promoteurs
- Tableau des paires Box-35 / Box-10 détectées
- Colonnes : position, séquence, mismatches, espacement, qualité
- Niveaux de qualité : `Consensus parfait` > `Fort` > `Modéré` > `Faible`

#### 🎯 Shine-Dalgarno
- Sites SD localisés avec leur ATG associé
- Colonnes : position SD, séquence, position ATG, espacement, mismatches, qualité

#### 🔚 Terminateurs
- Structures tige-boucle + polyT
- Colonnes : position, Bras1, Boucle, Bras2, PolyT, longueur tige

#### ✂️ Restriction
- Enzymes présentes dans la séquence, nombre de coupures et positions 1-based

#### 📊 Statistiques
- Longueur, GC%, AT%, ratio G/C, ratio A/T
- Tm Wallace (courtes séquences) et Tm Nearest-Neighbor
- Masse moléculaire simple brin et double brin (kDa)
- Barre de composition visuelle colorée

---

### Onglet Visualisation

Utiliser le menu déroulant pour choisir parmi 4 graphiques :

| Graphique | Description |
|---|---|
| Carte de séquence | Vue linéaire de la séquence avec ORFs, promoteurs, SD et terminateurs |
| GC% (fenêtre glissante) | Courbe du taux GC sur fenêtre de 100 pb avec ligne du GC% moyen |
| Composition nucléotidique | Histogramme + camembert A/T/G/C |
| Longueurs des ORFs | Barres horizontales des 20 premiers ORFs, colorées par brin |

> La barre d'outils Matplotlib permet de **zoomer**, **déplacer** et **enregistrer** en PNG.

---

## 💾 Export des résultats

Accessible via **Fichier** dans la barre de menu, après avoir lancé une analyse.

| Format | Contenu |
|---|---|
| 📊 Excel (`.xlsx`) | 6 onglets : ORFs, Promoteurs, SD, Terminateurs, Restriction, Statistiques |
| 📄 CSV (`.csv`) | Tous les résultats dans un fichier texte tabulé |
| 🗄️ JSON (`.json`) | Export complet structuré, idéal pour traitement programmatique |
| 📝 Rapport TXT (`.txt`) | Rapport lisible et formaté prêt à être partagé |
| 🧬 FASTA (`.fasta`) | Séquence analysée au format FASTA standard |

---

## 🏗 Architecture technique

### Flux de l'analyse

```
Saisie séquence
      │
      ▼
input_frame.py  ──►  app.py (_run_analysis)
                          │
                          ▼
                   Thread séparé (_worker)
                          │
              ┌───────────┼───────────────┐
              ▼           ▼               ▼
        orf_finder   motif_finder    statistics
              │           │               │
              └───────────┴───────────────┘
                          │
                          ▼
                  root.after(0, _finish_analysis)
                          │
              ┌───────────┴──────────────┐
              ▼                          ▼
       results_frame.display()    visualize_frame.draw()
```

### Threading

L'analyse est exécutée dans un **thread daemon** séparé du thread principal Tkinter. Les mises à jour de l'interface (barre de statut, affichage des résultats) sont toutes renvoyées dans le thread principal via `root.after(0, callback)`, respectant le modèle thread-safe de Tkinter.

```python
import threading

def _worker():
    # calculs lourds ici
    results = run_all_analyses(seq)
    # retour dans le thread principal
    self.root.after(0, lambda: self._finish_analysis(results))

threading.Thread(target=_worker, daemon=True).start()
```

### Algorithmes clés

**Recherche d'ORFs (6 frames)**
```
Pour chaque brin (sens, antisens) :
  Pour chaque décalage (0, 1, 2) :
    Lire les codons séquentiellement
    Identifier les segments entre codons stop
    Pour chaque segment : chercher le premier ATG
```

**Recherche de promoteurs**
```
Pour chaque position i dans la séquence :
  Calculer la distance de Hamming avec TATAAT (Box -10)
  Calculer la distance de Hamming avec TTGACA (Box -35)
  Si distance ≤ seuil :
    Chercher une paire valide avec espacement 14–22 pb
```

**Détection de terminateurs**
```
Pour chaque position i :
  Pour chaque longueur de tige (4–12) :
    Pour chaque longueur de boucle (3–8) :
      Vérifier si arm2 ≈ reverse_complement(arm1)
      Vérifier polyT après la tige-boucle
```

---

## 📦 Dépendances

```
matplotlib>=3.7.0       # Graphiques intégrés dans Tkinter
openpyxl>=3.1.0         # Export Excel (.xlsx)
```

> Les modules `tkinter`, `threading`, `csv`, `json`, `math`, `collections` sont inclus dans la bibliothèque standard Python.

Installation :
```bash
pip install matplotlib openpyxl
```

---

## 🧪 Séquence de test

La séquence suivante (issue du cahier des charges) est utilisée comme séquence de démonstration au lancement :

```
CATTTCTCTTAAGATTTATTCTATCTTAACACAACAACTTTTAATAAAAGATATGTAGAT
TACAATTTAAATAGATTGTAATATTTGTAACACTAACATTAATATAGTTGTTATTTTTGT
TACATAAACCACTAATAACTCATAATCTTTTAAAACTTATATTTGAGATAACATCAACTT
...
```

Elle contient plusieurs ORFs, des promoteurs, des sites Shine-Dalgarno et des terminateurs détectables, ce qui en fait un bon cas de démonstration.

---

## 👥 Auteurs

Développé dans le cadre du projet **Élaboration d'une application d'analyse d'ADN** — BioSeq Lab.

---

*DNA Analyzer — BioSeq Lab | Python / Tkinter / Matplotlib*
