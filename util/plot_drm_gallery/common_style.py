#!/usr/bin/env python
"""
Common styles, functions, etc.

"""

import pandas as pd

# Helper function to create title from t_info
def plot_make_title(t_info):
    """Create plot title from metadata"""
    if t_info.empty:
        rv = ''
    elif 'experiment' in t_info:
        exp_name = str.strip(t_info['experiment'].iloc[0])
        if len(exp_name) < 50:
            chaser = f", Ensemble Size {t_info['ensemble_size'].iloc[0]}"
        else:
            chaser = ''
        rv = exp_name + chaser
    else:
        rv = ''
    return rv
    
