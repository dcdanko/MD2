import pandas as pd
from os.path import dirname, join

# Biological constants
ROOT_RANK = 'root'
RANK_LIST = [
    'forma', 'varietas', 'subspecies', 'species', 'species group', 'species subgroup', 'subgenus',
    'genus', 'tribe', 'subfamily', 'family', 'suborder', 'cohort', 'order', 'infraclass', 'subclass',
    'class', 'subphylum', 'phylum', 'subkingdom', 'kingdom', 'superkingdom', 'no rank', ROOT_RANK
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

NEG = 'Negative'
POS = 'Positive'
PHYLUM_GRAM_STAINS = {
    'Bacteroidetes': NEG,
    'Chlorobi': NEG,
    'Cyanobacteria': NEG,
    'Proteobacteria': NEG,
    'Firmicutes': POS,
    'Deinococcus-Thermus': POS,
    'Actinobacteria': POS,
    'Fusobacteria': NEG,
    'Nitrospirae': None,
    'Acidobacteria': None,
    'Fibrobacteres': NEG,
    'Thermodesulfobacteria': None,
    'Caldiserica': None,
    'Armatimonadetes': NEG,
    'Dictyoglomi': NEG,
    'Deferribacteres': NEG,
    'Verrucomicrobia': None,
    'Chrysiogenetes': None,
    'Kiritimatiellaeota': None,
    'Aquificae': NEG,
    'Thermotogae': NEG,
    'Chloroflexi': NEG,
    'Planctomycetes': NEG,
    'Spirochaetes': None,
    'Chlamydiae': NEG,
    'Lentisphaerae': None,
    'Synergistetes': NEG,
    'Tenericutes': None,
    'Ignavibacteriae': None,
    'Nitrospinae': None,
    'Rhodothermaeota': None,
    'Calditrichaeota': None,
    'Balneolaeota': None,
    'Abditibacteriota': None,
    'Coprothermobacterota': None,
}


VIRUS = 'Viruses'
BACTERIA = 'Bacteria'
FUNGI = 'Eukaryota'
ALLOWED_SUPERKINGDOMS = [VIRUS, BACTERIA, FUNGI]
DOMAINS = [el.lower() for el in ALLOWED_SUPERKINGDOMS]

# Programmatic
MICROBE_DIRECTORY_CSV = join(
    dirname(dirname(__file__)),
    'tables/microbe-directory.csv'
)
MICROBE_DIRECTORY = pd.read_csv(MICROBE_DIRECTORY_CSV, index_col=7)
