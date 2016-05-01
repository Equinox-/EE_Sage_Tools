def tabular(data, cols):
    if cols == None:
        cpy = []
    else:
        cpy = [cols];
    cpy += data;
    return table(cpy, frame=true);

def tabular(data, cols, units, precise, smart=False):
    if smart:
        si_fn=SI_smart
    else:
        si_fn=SI
    if cols == None:
        cpy = []
    else:
        cpy = [cols];
    for row in data:
        rcpy = []
        for j in range(0, len(row)):
            if isinstance(row[j], str):
                rcpy.append(row[j])
            else:
                rcpy.append(si_fn(row[j], units[j], precise[j] if isinstance(precise, list) else precise))
        cpy.append(rcpy)
    return table(cpy, frame=true);

def tabularO(obj, precise, smart=False):
    return tabular(obj['data'], obj['cols'], obj['units'], precise, smart)