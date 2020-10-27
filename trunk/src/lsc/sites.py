extinction = {
'lsc': {'U':0.46, 'u':0.46, 'B':0.27, 'g':0.20 ,'V':0.12, 'r':0.09, 'R':0.09, 'i':0.02, 'I':0.02, 'z':0.03, 'J':0.07, 'H':0.023,'K':0.045, 'w':0.09},
'coj': {'U':0.63, 'u':0.70, 'B':0.32, 'g':0.26 ,'V':0.18, 'r':0.15, 'R':0.13, 'i':0.08, 'I':0.07, 'z':0.06, 'J':0.00, 'H':0.00, 'K':0.00, 'w':0.15},
'ogg': {'U':0.45, 'u':0.48, 'B':0.21, 'g':0.16 ,'V':0.12, 'r':0.09, 'R':0.07, 'i':0.04, 'I':0.03, 'z':0.03, 'J':0.00, 'H':0.00, 'K':0.00, 'w':0.09},
'elp': {'U':0.51, 'u':0.51, 'B':0.23, 'g':0.20 ,'V':0.15, 'r':0.10, 'R':0.10, 'i':0.05, 'I':0.05, 'z':0.05, 'J':0.00, 'H':0.00, 'K':0.00, 'w':0.10},
'cpt': {'U':0.51, 'u':0.51, 'B':0.23, 'g':0.20 ,'V':0.15, 'r':0.10, 'R':0.10, 'i':0.05, 'I':0.05, 'z':0.05, 'J':0.00, 'H':0.00, 'K':0.00, 'w':0.10},
'tfn': {'U':0.46, 'u':0.46, 'B':0.22, 'g':0.16, 'V':0.12, 'r':0.08, 'R':0.18, 'i':0.04, 'I':0.04, 'z':0.06, 'J':0.12, 'H':0.06, 'K':0.09, 'w':0.08},
None:  {'U':0.50, 'u':0.52, 'B':0.25, 'g':0.20, 'V':0.14, 'r':0.10, 'R':0.11, 'i':0.05, 'I':0.04, 'z':0.05, 'J':0.03, 'H':0.01, 'K':0.02, 'w':0.10}
} # use average extinction values if no telescope listed (probably merged images)
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
'w': ['w']
}

filterst1 = dict()
for key, val in filterst.items():
    for subval in val:
        filterst1[subval] = key

filterst['landolt'] = sum([filterst[f] for f in 'UBVRI'], [])
filterst['sloan'] = sum([filterst[f] for f in 'ugrizw'], [])
filterst['apass'] = sum([filterst[f] for f in 'BVgriw'], [])
filterst['gaia'] = sum([filterst[f] for f in filterst.keys()], [])
filterst[''] = filterst['landolt'] + filterst['sloan']

def chosecolor(allfilter, usegood=False):
    color = {filt: [] for filt in allfilter}
    for col in ['UB', 'BV', 'VR', 'RI', 'ug', 'gr', 'ri', 'iz']:
        if col[0] in allfilter and col[1] in allfilter:
            color[col[0]].append(col)
            color[col[1]].append(col)
    if usegood:
        goodcol = {'U': 'UB', 'B': 'BV', 'V': 'VR', 'R': 'VR', 'I': 'RI',
                   'u': 'ug', 'g': 'gr', 'r': 'ri', 'i': 'ri', 'z': 'iz'}
        for filt in color:
            if filt in goodcol and goodcol[filt] in color[filt]:
                color[filt] = [goodcol[filt]]
    return color
