import os

libs = []
for root, d, files in os.walk("./"):
    for f in files:
        if f.endswith(".sage") and not(f == "lib.sage"):
            libs.append(os.path.join(root, f))


with open('lib.sage', 'w') as o:
    for f in libs:
        o.write("\n## begin " + f + " ##\n")
        with open(f) as i:
            for l in i:
                o.write(l)
        o.write("\n## end " + f + " ##\n")