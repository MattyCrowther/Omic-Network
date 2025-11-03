import numpy as np
import pandas as pd

def to_primitive(x):
    """Convert numpy/pandas scalars to Python primitives."""
    if pd.isna(x):
        return None
    if isinstance(x, (np.generic,)):  # catches np.int64, np.float64, np.bool_, etc.
        return x.item()
    return x