import pandas as pd
from os.path import dirname, join


MICROBE_DIRECTORY_CSV = join(
    dirname(dirname(__file__)),
    'tables/microbe-directory.csv'
)
MICROBE_DIRECTORY = pd.read_csv(MICROBE_DIRECTORY_CSV, index_col=7)
