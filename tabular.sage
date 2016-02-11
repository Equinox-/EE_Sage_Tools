def tabular(data, cols):
    if cols == None:
        cpy = []
    else:
        cpy = [cols];
    cpy += data;
    return table(cpy, frame=true);

def tabular(data, cols, units, precise):
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
                rcpy.append(SI(row[j], units[j], precise[j]))
        cpy.append(rcpy)
    return table(cpy, frame=true);