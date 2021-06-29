import re

def _human_key(key):
    parts = re.split('(\d*\.\d+|\d+)', key)
    return tuple((e.swapcase() if i % 2 == 0 else float(e))
            for i, e in enumerate(parts))