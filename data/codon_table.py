# ============================================================
# Table des codons standard (code génétique universel)
# ============================================================

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

STOP_CODONS = ['TAA', 'TAG', 'TGA']
START_CODON = 'ATG'

# Noms complets des acides aminés
AA_NAMES = {
    'A': 'Alanine',    'R': 'Arginine',   'N': 'Asparagine',
    'D': 'Aspartate',  'C': 'Cystéine',   'Q': 'Glutamine',
    'E': 'Glutamate',  'G': 'Glycine',    'H': 'Histidine',
    'I': 'Isoleucine', 'L': 'Leucine',    'K': 'Lysine',
    'M': 'Méthionine', 'F': 'Phénylalanine', 'P': 'Proline',
    'S': 'Sérine',     'T': 'Thréonine',  'W': 'Tryptophane',
    'Y': 'Tyrosine',   'V': 'Valine',     '*': 'Stop',
}

# Propriétés des acides aminés
AA_HYDROPHOBIC = set('AVILMFWP')
AA_POLAR       = set('STNQYC')
AA_CHARGED     = set('KRHDEP')
AA_SPECIAL     = set('GPC')

# Sites de restriction courants
RESTRICTION_SITES = {
    'EcoRI':   ('GAATTC', 1),
    'BamHI':   ('GGATCC', 1),
    'HindIII': ('AAGCTT', 1),
    'SalI':    ('GTCGAC', 1),
    'XhoI':    ('CTCGAG', 1),
    'NcoI':    ('CCATGG', 1),
    'XbaI':    ('TCTAGA', 1),
    'SmaI':    ('CCCGGG', 3),
    'KpnI':    ('GGTACC', 5),
    'PstI':    ('CTGCAG', 5),
    'SphI':    ('GCATGC', 5),
    'ClaI':    ('ATCGAT', 2),
    'NotI':    ('GCGGCCGC', 2),
    'XmaI':    ('CCCGGG', 1),
    'MluI':    ('ACGCGT', 1),
}
