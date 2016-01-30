def tabular(data, cols):
    cpy = [cols];
    cpy += data;
    return table(cpy, frame=true);

def tabular(data, cols, units, precise):
    cpy = [cols];
    for row in data:
        rcpy = []
        for j in range(0, len(row)):
            rcpy.append(SI(row[j], units[j], precise[j]))
        cpy.append(rcpy)
    return table(cpy, frame=true);