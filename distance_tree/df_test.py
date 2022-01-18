#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import os
import csv
import sys
import pandas as pd

new_df = pd.DataFrame(
    [[1.0, .22352, 3.0], [4.12441, .22352, 8.0]], index=['A', 'BB'])

new_df.rename(index=lambda index_: index_.ljust(10), inplace=True)

# Write the remaining matrix, omit the column (header) names
new_df.to_csv(sys.stdout, sep='\t', header=False,
              index=True, index_label=False, float_format="%.4f")
sys.stdout.flush()
