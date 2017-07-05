extinction = {
'lsc': {'U':0.46, 'u':0.46, 'B':0.27, 'g':0.20 ,'V':0.12, 'r':0.09, 'R':0.09, 'i':0.02, 'I':0.02, 'z':0.03, 'J':0.07, 'H':0.023,'K':0.045},
'coj': {'U':0.63, 'u':0.70, 'B':0.32, 'g':0.26 ,'V':0.18, 'r':0.15, 'R':0.13, 'i':0.08, 'I':0.07, 'z':0.06, 'J':0.00, 'H':0.00, 'K':0.00},
'ogg': {'U':0.45, 'u':0.48, 'B':0.21, 'g':0.16 ,'V':0.12, 'r':0.09, 'R':0.07, 'i':0.04, 'I':0.03, 'z':0.03, 'J':0.00, 'H':0.00, 'K':0.00},
'elp': {'U':0.51, 'u':0.51, 'B':0.23, 'g':0.20 ,'V':0.15, 'r':0.10, 'R':0.10, 'i':0.05, 'I':0.05, 'z':0.05, 'J':0.00, 'H':0.00, 'K':0.00},
'cpt': {'U':0.51, 'u':0.51, 'B':0.23, 'g':0.20 ,'V':0.15, 'r':0.10, 'R':0.10, 'i':0.05, 'I':0.05, 'z':0.05, 'J':0.00, 'H':0.00, 'K':0.00},
'tfn': {'U':0.46, 'u':0.46, 'B':0.22, 'g':0.16, 'V':0.12, 'r':0.08, 'R':0.18, 'i':0.04, 'I':0.04, 'z':0.06, 'J':0.12, 'H':0.06, 'K':0.09},
}
extinction['PS1'] = extinction['ogg']
extinction['SDSS'] = extinction['elp']

filterst = {
'U': ['U', 'Astrodon-U'],
'B': ['B', 'Bessell-B'],
'V': ['V', 'Bessell-V'],
'R': ['R', 'Bessell-R'],
'I': ['I', 'Bessell-I'],
'u': ['up', 'SDSS-U'],
'g': ['gp', 'SDSS-G'],
'r': ['rp', 'SDSS-R'],
'i': ['ip', 'SDSS-I'],
'z': ['zs', 'Pan-Starrs', 'Pan-Starrs-Z'],
'landolt': ['landolt'],
'sloan': ['sloan'],
'apass': ['apass']
}

filterst1 = dict()
for key, val in filterst.items():
    for subval in val:
        filterst1[subval] = key
