# 🧬 DNA Analyzer — BioSeq Lab

> Application Python d'analyse bioinformatique de séquences ADN bactériennes.
> Interface graphique complète, analyse double brin, 7 types d'analyses, visualisations interactives et export multi-format.

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?style=flat-square&logo=python&logoColor=white)
![Tkinter](https://img.shields.io/badge/GUI-Tkinter-1B7A4A?style=flat-square)
![Matplotlib](https://img.shields.io/badge/Plots-Matplotlib-orange?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)

---

## 📋 Table des matières

- [Aperçu](#-aperçu)
- [Fonctionnalités](#-fonctionnalités)
- [Nouveautés & Améliorations](#-nouveautés--améliorations)
- [Structure du projet](#-structure-du-projet)
- [Installation](#-installation)
- [Lancement](#-lancement)
- [Guide d'utilisation](#-guide-dutilisation)
- [Analyse Double Brin](#-analyse-double-brin-nouvelle-fonctionnalité)
- [Description des analyses](#-description-des-analyses)
- [Export des résultats](#-export-des-résultats)
- [Architecture technique](#-architecture-technique)
- [Dépendances](#-dépendances)

---

## 🔭 Aperçu

DNA Analyzer est une application de bureau développée en **Python / Tkinter** permettant d'analyser des séquences ADN bactériennes à partir d'une saisie directe ou d'un fichier FASTA.

Elle couvre les analyses bioinformatiques fondamentales avec une architecture modulaire :

- **Détection des ORFs** sur les 6 cadres de lecture depuis le 1er nucléotide (brin sens et antisens)
- **Identification de l'ORF codant** le plus probable avec score /100 et traduction protéique
- **Analyse double brin simultanée** avec schéma visuel coloré et résultat optimal
- **Recherche de motifs régulateurs** : promoteurs (-10 / -35), terminateurs Rho-indépendants GC-riches, Shine-Dalgarno
- **Analyses physico-chimiques** : GC%, Tm, masse moléculaire
- **Visualisations graphiques** interactives via Matplotlib
- **Export** CSV, Excel, JSON, TXT, FASTA

![Interface principale — détection ORF](screenshots/screenshot_orf_main.png)
*Interface principale : détection des ORFs sur les 6 cadres de lecture, meilleur candidat identifié (Frame +1, score ★)*

---

## ✨ Fonctionnalités

### Analyses biologiques
| Analyse | Description |
|---|---|
| 📋 Cadres de lecture | Segmentation codon par codon sur les 6 frames depuis le 1er nucléotide |
| 🏆 ORF codant | Score /100 : longueur, ATG, taille protéique, brin, stop — meilleur candidat identifié |
| ⚙️ Traduction protéique | Traduction de l'ORF principal, colorée par propriété biochimique |
| 📍 Promoteurs | Pairs Box-35 / Box-10 avec tolérance aux mismatches |
| 🎯 Shine-Dalgarno | Consensus AGGAGG localisé 7–9 pb en amont d'ATG |
| 🔚 Terminateurs | Structures tige-boucle **riches en GC** (tige ≥ 50% GC) + polyT |
| ✂️ Sites de restriction | 15 enzymes courantes (EcoRI, BamHI, HindIII…) |
| 📊 Statistiques | GC%, Tm Wallace, Tm Nearest-Neighbor, masse moléculaire |

### Interface & usabilité
- Interface graphique avec **thème biologique** (palette verte émeraude)
- Animation ADN double hélice en temps réel dans la barre de titre
- **Validation des entrées** en temps réel (vide, trop court, caractères non-ATCG)
- **Analyse threadée** — UI non bloquante sur les longues séquences FASTA
- Import de fichiers FASTA (`.fasta`, `.fa`, `.txt`, `.seq`)
- Calcul du **brin complémentaire inverse** en un clic
- Seuil de longueur minimum des ORFs configurable (défaut : 30 pb)

---

## 🆕 Nouveautés & Améliorations

### Correctifs critiques

**Correctif A1 — Bouton d'analyse**
Le bouton `LANCER L'ANALYSE` manquait son argument `command=self._run` et ne répondait pas aux clics. Corrigé dans `gui/input_frame.py`.

**Correctif A2 — Seuil longueur minimale ORF**
La valeur par défaut (90 pb) était trop restrictive pour les séquences courtes. Réduite à **30 pb** dans `orf_finder.py` et `input_frame.py`.

### Nouvelles fonctionnalités

**Amélioration A3 — Validation des entrées utilisateur**
Trois vérifications avant toute analyse :
- Séquence vide → message d'erreur rouge inline
- Séquence < 15 nucléotides → avertissement avec la longueur
- Caractères non-ATCG → liste des caractères invalides

Validation aussi **en temps réel** pendant la saisie (bordure rouge si caractères invalides). Les messages disparaissent automatiquement après 5 secondes.

**Amélioration A4 — Refonte des onglets ORF**
Deux nouveaux onglets distincts remplacent l'onglet ORF unique :
- **📋 Cadres de lecture** : tous les segments entre codons stop sur les 6 cadres depuis le 1er nucléotide, avec indication de la présence d'un ATG
- **🏆 ORF Codant** : meilleur candidat avec score /100, justification détaillée, séquence ADN colorée et protéine traduite

**Amélioration A5 — Filtre biologique des terminateurs**
Les terminateurs Rho-indépendants actifs ont une tige riche en G-C. Seuls les terminateurs avec **GC tige ≥ 50%** sont retenus dans les analyses principales (les autres restent disponibles en fallback).

**Amélioration A6 — Fenêtre Analyse Double Brin** ⭐
Voir la section dédiée ci-dessous.

---

## 📁 Structure du projet

```
dna_analyzer/
│
├── main.py                    # Point d'entrée — lancer avec : python main.py
│
├── analysis/
│   ├── __init__.py
│   ├── orf_finder.py          # ORFs 6 frames + score codant [AMÉLIORÉ]
│   ├── motif_finder.py        # Promoteurs, SD, terminateurs GC-riche, restriction
│   └── statistics.py         # Composition, Tm, masse moléculaire, GC glissant
│
├── data/
│   ├── __init__.py
│   └── codon_table.py        # Table des codons + sites de restriction
│
├── gui/
│   ├── __init__.py
│   ├── app.py                # Fenêtre principale + orchestration analyses [AMÉLIORÉ]
│   ├── input_frame.py        # Saisie, validation, boutons action [AMÉLIORÉ]
│   ├── results_frame.py      # 7 onglets de résultats [AMÉLIORÉ]
│   ├── visualize_frame.py    # Graphiques Matplotlib
│   └── dual_strand_window.py # Fenêtre analyse double brin [NOUVEAU]
│
├── export/
│   ├── __init__.py
│   └── exporter.py           # Export CSV, Excel, JSON, TXT, FASTA
│
└── requirements.txt          # Dépendances Python
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
│  [🧬  SCHÉMA DOUBLE BRIN]       │
└─────────────────────────────────┘
```

| Élément | Rôle |
|---|---|
| Zone de texte | Coller la séquence ADN directement (ATCG uniquement) — validation en temps réel |
| `📁 Importer FASTA` | Charger un fichier FASTA — les lignes `>` sont ignorées |
| `↔ RC` | Remplace la séquence par son brin complémentaire inverse |
| `✕` | Efface la zone de saisie |
| Cases à cocher | Activer/désactiver chaque analyse individuellement |
| Longueur min. ORF | Seuil en pb — `30` pour tous les ORFs, `90+` pour les candidats biologiques |
| `▶ LANCER L'ANALYSE` | Démarre l'analyse avec validation des entrées |
| `🧬 SCHÉMA DOUBLE BRIN` | **[NOUVEAU]** Ouvre la fenêtre d'analyse visuelle des 2 brins |

> **Validation automatique :** si la séquence est vide, trop courte ou contient des caractères non-ATCG, un message d'erreur rouge s'affiche et l'analyse ne démarre pas.

---

### Panneau droit — Onglets de résultats

#### 📋 Cadres de lecture *(nouveau)*
Affiche tous les segments entre codons stop sur les **6 cadres de lecture** en partant du 1er nucléotide, sur les deux brins. Chaque ligne indique :
- Cadre (+1 à +3 sens, -1 à -3 antisens) et brin
- Position de début et de fin (1-based)
- Longueur en pb, nombre d'acides aminés
- Présence d'un ATG interne (★)
- Codon stop terminal

Le **meilleur segment** (plus long avec ATG) est surligné en vert.

#### 🏆 ORF Codant *(nouveau)*
Identifie le meilleur candidat codant avec :
- **Score /100** calculé sur 7 critères biologiques
- Justification détaillée point par point
- Cartes de résumé : cadre, brin, position, longueur, ATG, stop, masse moléculaire
- Séquence ADN de l'ORF (par triplets)
- Protéine traduite colorée (hydrophobe/polaire/chargé/autre)
- Classement des 10 meilleurs candidats

#### 📍 Promoteurs
Tableau des paires Box-35 / Box-10 détectées. Niveaux de qualité : `Consensus parfait` > `Fort` > `Modéré` > `Faible`.

#### 🎯 Shine-Dalgarno
Sites SD localisés avec leur ATG associé, espacement et nombre de mismatches.

#### 🔚 Terminateurs
Structures tige-boucle (GC-riches) + polyT. Colonnes : position, Bras1, Boucle, Bras2, PolyT, longueur tige.

#### ✂️ Restriction
Enzymes présentes dans la séquence, nombre de coupures et positions 1-based.

#### 📊 Statistiques
Longueur, GC%, AT%, ratio G/C, Tm Wallace, Tm Nearest-Neighbor, masse moléculaire ss/ds, barre de composition visuelle colorée.

---

## 🧬 Analyse Double Brin *(nouvelle fonctionnalité)*

Accessible via le bouton **"🧬 SCHÉMA DOUBLE BRIN"** dans le panneau gauche.

### Principe biologique
Une molécule d'ADN double brin possède deux orientations de lecture :
- **Brin sens 5'→3' (+)** — la séquence telle que saisie
- **Brin antisens 3'→5' (-)** — le brin complémentaire inverse

Un gène peut être codé sur l'un ou l'autre brin. Cette fenêtre analyse **simultanément** les deux orientations et identifie le meilleur candidat.

![Carte double brin — analyse simultanée des deux brins] ("2d137923-fcf8-4cee-995b-a141a202f439.png")
*Carte linéaire des deux brins : ORFs (vert/orange), promoteurs Box-35 (bleu), Box-10 (orange), terminateurs GC-riches. Le brin sens (+) est identifié comme meilleur candidat (score 90/100).*

### Onglets de la fenêtre

#### 📊 Carte Double Brin
Représentation linéaire des deux brins sur deux axes Matplotlib superposés, dans le même style que la carte de séquence principale :
- 🟩 **Rectangles verts/teal** — ORFs (or pour le meilleur)
- 🔵 **Lignes bleues** — promoteurs Box-35
- 🔴 **Lignes rouges** — promoteurs Box-10
- 🟠 **Zones orange** — terminateurs GC-riches
- 🔺 **Triangles violets** — sites Shine-Dalgarno
- Barre d'outils Matplotlib (zoom, sauvegarde PNG)

#### 🔤 Séquence Brin + et Brin −
Chaque nucléotide est coloré selon sa fonction biologique :

| Couleur | Élément |
|---|---|
| Vert foncé `#A9DFBF` | ORF meilleur candidat |
| Vert pâle `#D5F5E3` | ORF classique |
| Bleu pâle `#AED6F1` | Promoteur Box -35 |
| Rouge pâle `#F5B7B1` | Promoteur Box -10 |
| Orange pâle `#FDEBD0` | Terminateur Rho-indépendant (GC-riche) |
| Violet pâle `#E8DAEF` | Site Shine-Dalgarno |
| **Texte vert gras** | Codon ATG |
| **Texte rouge gras** | Codon Stop |

Chaque ligne de 60 bases affiche aussi les annotations de position à droite (`─35@pos`, `TERM@pos(GC:62%)`).

#### ★ Résultat Optimal

![Résultat optimal — analyse double brin](screenshots/screenshot_optimal_result.png)
*Résultat optimal : brin sens (+) sélectionné avec score 90/100. Tableau comparatif brin (+) vs brin (−), protéine traduite de 47 acides aminés colorée par propriété biochimique.*

- **Bannière or** avec les 12 métriques clés du brin gagnant
- **Mini-carte** du brin optimal avec légende
- **Protéine traduite** colorée par propriété biochimique
- **Tableau comparatif** brin (+) vs brin (-) avec cellules gagnantes surlignées en vert

---

### Onglet Visualisation (panneau principal)

![Carte de séquence — vue linéaire](screenshots/screenshot_sequence_map.png)
*Vue linéaire de la séquence (300 pb) : ORF détecté en vert (positions 202–234), promoteurs Box-35 (lignes bleues) et Box-10 (lignes orange) répartis sur la séquence.*

Utiliser le menu déroulant pour choisir parmi 4 graphiques :

| Graphique | Description |
|---|---|
| Carte de séquence | Vue linéaire avec ORFs, promoteurs, SD et terminateurs |
| GC% (fenêtre glissante) | Courbe du taux GC sur fenêtre de 100 pb avec ligne du GC% moyen |
| Composition nucléotidique | Histogramme + camembert A/T/G/C |
| Longueurs des ORFs | Barres horizontales des 20 premiers ORFs, colorées par brin |

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
input_frame.py ──► Validation (vide / trop court / non-ATCG)
      │
      ▼
app.py (_run_analysis)
      │
      ├──► find_reading_frames()  →  find_best_coding_orf()   [6 cadres + score]
      ├──► find_orfs()                                         [ATG→Stop classique]
      ├──► find_promoters()                                    [Box-35 / Box-10]
      ├──► find_shine_dalgarno()                               [AGGAGG]
      ├──► find_terminators() + filtre GC ≥ 50%               [tige-boucle]
      ├──► find_restriction_sites()                            [15 enzymes]
      └──► nucleotide_composition() + Tm + MW                 [statistiques]
               │
               ▼
      results_frame.display()    visualize_frame.draw()
```

### Score de potentiel codant

```python
score = 0
score += min(40, longueur_orf / longueur_seq * 200)  # longueur relative
score += 20 if has_atg else 0                         # codon start
score += 20 if num_aa >= 50 else 10 if num_aa >= 30 else 0  # taille
score += 5 if strand == '+' else 0                    # brin sens
score += 10 if stop_codon present else 0              # ORF complet
score += 5 if '*' not in protein else 0               # pas de stop interne
```

### Threading

```python
import threading

def _worker():
    results = run_all_analyses(seq)
    self.root.after(0, lambda: self._finish_analysis(results))

threading.Thread(target=_worker, daemon=True).start()
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

La séquence suivante (300 pb) est utilisée comme séquence de démonstration. Elle contient un ORF principal sur le **brin antisens** (score 90/100), des promoteurs, et des terminateurs GC-riches :

```
CATTTCTCTTAAGATTTATTCTATCTTAACACAACAACTTTTAATAAAAGATATGTAGAT
TACAATTTAAATAGATTGTAATATTTGTAACACTAACATTAATATAGTTGTTATTTTTGT
TACATAAACCACTAATAACTCATAATCTTTTAAAACTTATATTTGAGATAACATCAACTT
TACATTACAAGTTATAAAACAAAAGAAGTGGGACACAGAATTCGTCTTGAACACTGTGTC
CCACCTCGTCCCCAAAACTTGCTCTGTCCGTAGAAAAATAAAAAGGGGCCCCCTTTGTTG
```

**Résultat attendu :**

| Critère | Brin + | Brin − |
|---|---|---|
| Score codant | 76/100 | **90/100 ★** |
| Meilleur ORF | 54 pb / 18 aa | **141 pb / 47 aa** |
| Codon ATG | ✓ | ✓ |
| Terminateurs GC-riche | 2 | 2 |

---

## 👥 Auteurs

me xoxo

---

*DNA Analyzer — BioSeq Lab | Python / Tkinter / Matplotlib*
