import pandas as pd
from os.path import dirname, join

# Biological constants
ROOT_RANK = 'root'
RANK_LIST = [
    'subspecies', 'species', 'species group', 'species subgroup', 'subgenus',
    'genus', 'subfamily', 'family', 'suborder', 'order', 'subclass', 'class',
    'subphylum', 'phylum', 'kingdom', 'superkingdom', 'no rank', 'varietas',
    'forma', 'tribe', ROOT_RANK
]
CANONICAL_RANKS = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
SPORE_FORMING_GENERA = [
    'Bacillus', 'Clostridium', 'Sporolactobacillus', 'Sporosarcina',
    'Cerasibacillus', 'Alkalibacillus', 'Amphibacillus', 'Anoxybacillus',
    'Filobacillus', 'Geobacillus', 'Gracilibacillus', 'Halobacillus',
    'Lentibacillus', 'Oceanobacillus', 'Paraliobacillus', 'Pontibacillus',
    'Tenuibacillus', 'Thalassobacillus', 'Virgibacillus', 'Jeotgalibacillus',
    'Marinibacillus', 'Planomicrobium', 'Ureibacillus', 'Sporobacterium',
    'Desulfitobacterium', 'Desulfonispora', 'Desulfosporosinus',
    'Desulfotomaculum', 'Sporotomaculum', 'Syntrophobotulus', 'Pelotomaculum',
    'Sporotomaculum', 'Syntrophobotulus', 'Thermincola', 'Filifactor',
    'Tepidibacter', 'Anaerotruncus', 'Oscillospira', 'Sporobacter',
]
NON_SPORE_FORMING_GENERA = [
    'Halolactibacillus', 'Marinococcus', 'Saccharococcus', 'Planococcus',
    'Caryophanon', 'Filibacter', 'Kurthia', 'Peptococcus', 'Cryptanaerobacter',
    'Dehalobacter', 'Anaerofilum', 'Acetivibrio', 'Acetanaerobacterium',
    'Fastidiosipila', 'Papillibacter', 'Subdoligranulum',
]
VIRUS = 'Virus'
BACTERIA = 'Bacteria'
FUNGI = 'Fungi'
ALLOWED_SUPERKINGDOMS = [VIRUS, BACTERIA, FUNGI]
DOMAINS = [el.lower() for el in ALLOWED_SUPERKINGDOMS]

# Programmatic
MICROBE_DIRECTORY_CSV = join(
    dirname(dirname(__file__)),
    'tables/microbe-directory.csv'
)
MICROBE_DIRECTORY = pd.read_csv(MICROBE_DIRECTORY_CSV, index_col=7)
