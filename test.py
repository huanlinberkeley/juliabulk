import re

with open("param.txt", 'r') as f:
    lines = f.readlines()

len1, len2 = 25, 30

with open("out.jl", 'w') as f:
    for line in lines:
        if re.match(r'`BPRco', line):
            newline = line[line.find("(")+1 : line.rfind(")")]
            par, val, uni, lwr, upr, des = re.split(r"\s*,\s*", newline, 5)
            par = par.strip()
            val = val.strip('\`')
            des = des.strip().strip('\"')
            assign = f"{par} = {val}"
            if len(assign) < len1:
                f.write(assign.ljust(len1))
            else:
                f.write(assign.ljust(len2))
            f.write(f"# {des}, range ({lwr}, {upr}), unit {uni}\n")
        elif re.match(r'`BPRoc', line):
            newline = line[line.find("(")+1 : line.rfind(")")]
            par, val, uni, lwr, upr, des = re.split(r"\s*,\s*", newline, 5)
            par = par.strip()
            des = des.strip().strip('\"')
            assign = f"{par} = {val}"
            if len(assign) < len1:
                f.write(assign.ljust(len1))
            else:
                f.write(assign.ljust(len2))
            f.write(f"# {des}, range ({lwr}, {upr}), unit {uni}\n")
        elif re.match(r'`BPRoz', line):
            newline = line[line.find("(")+1 : line.rfind(")")]
            par, val, uni, des = re.split(r"\s*,\s*", newline, 3)
            par = par.strip()
            des = des.strip().strip('\"')
            assign = f"{par} = {val}"
            if len(assign) < len1:
                f.write(assign.ljust(len1))
            else:
                f.write(assign.ljust(len2))
            f.write(f"# {des}, unit {uni}\n")
        elif re.match(r'`BPRcz', line):
            newline = line[line.find("(")+1 : line.rfind(")")]
            par, val, uni, des = re.split(r"\s*,\s*", newline, 3)
            par = par.strip()
            des = des.strip().strip('\"')
            assign = f"{par} = {val}"
            if len(assign) < len1:
                f.write(assign.ljust(len1))
            else:
                f.write(assign.ljust(len2))
            f.write(f"# {des}, unit {uni}\n")
        elif re.match(r'`BPRnb', line):
            newline = line[line.find("(")+1 : line.rfind(")")]
            par, val, uni, des = re.split(r"\s*,\s*", newline, 3)
            par = par.strip()
            des = des.strip().strip('\"')
            assign = f"{par} = {val}"
            if len(assign) < len1:
                f.write(assign.ljust(len1))
            else:
                f.write(assign.ljust(len2))
            f.write(f"# {des}, unit {uni}\n")
        elif re.match(r'`BPRoo', line):
            newline = line[line.find("(")+1 : line.rfind(")")]
            par, val, uni, lwr, upr, des = re.split(r"\s*,\s*", newline, 5)
            par = par.strip()
            des = des.strip().strip('\"')
            assign = f"{par} = {val}"
            if len(assign) < len1:
                f.write(assign.ljust(len1))
            else:
                f.write(assign.ljust(len2))
            f.write(f"# {des}, range ({lwr}, {upr}), unit {uni}\n")
        elif re.match(r'`BPIcc', line):
            newline = line[line.find("(")+1 : line.rfind(")")]
            par, val, uni, lwr, upr, des = re.split(r"\s*,\s*", newline, 5)
            par = par.strip()
            des = des.strip().strip('\"')
            assign = f"{par} = {val}"
            if len(assign) < len1:
                f.write(assign.ljust(len1))
            else:
                f.write(assign.ljust(len2))
            f.write(f"# {des}, range ({lwr}, {upr}), unit {uni}\n")
        elif re.match(r'`IPIco', line):
            newline = line[line.find("(")+1 : line.rfind(")")]
            par, val, uni, lwr, upr, des = re.split(r"\s*,\s*", newline, 5)
            par = par.strip()
            des = des.strip().strip('\"')
            assign = f"{par} = {val}"
            if len(assign) < len1:
                f.write(assign.ljust(len1))
            else:
                f.write(assign.ljust(len2))
            f.write(f"# {des}, range ({lwr}, {upr}), unit {uni}\n")
        elif re.match(r'`MPRnb', line):
            newline = line[line.find("(")+1 : line.rfind(")")]
            par, val, uni, des = re.split(r"\s*,\s*", newline, 3)
            par = par.strip()
            des = des.strip().strip('\"')
            assign = f"{par} = {val}"
            if len(assign) < len1:
                f.write(assign.ljust(len1))
            else:
                f.write(assign.ljust(len2))
            f.write(f"# {des}, unit {uni}\n")
        elif re.match(r'`MPRex', line):
            newline = line[line.find("(")+1 : line.rfind(")")]
            par, val, uni, exc, des = re.split(r"\s*,\s*", newline, 4)
            par = par.strip()
            des = des.strip().strip('\"')
            assign = f"{par} = {val}"
            if len(assign) < len1:
                f.write(assign.ljust(len1))
            else:
                f.write(assign.ljust(len2))
            f.write(f"# {des}, except {exc}, unit {uni}\n")
        elif re.match(r'`MPRcc', line):
            newline = line[line.find("(")+1 : line.rfind(")")]
            par, val, uni, lwr, upr, des = re.split(r"\s*,\s*", newline, 5)
            par = par.strip()
            des = des.strip().strip('\"')
            assign = f"{par} = {val}"
            if len(assign) < len1:
                f.write(assign.ljust(len1))
            else:
                f.write(assign.ljust(len2))
            f.write(f"# {des}, range ({lwr}, {upr}), unit {uni}\n")
        elif re.match(r'`MPRoo', line):
            newline = line[line.find("(")+1 : line.rfind(")")]
            par, val, uni, lwr, upr, des = re.split(r"\s*,\s*", newline, 5)
            par = par.strip()
            des = des.strip().strip('\"')
            assign = f"{par} = {val}"
            if len(assign) < len1:
                f.write(assign.ljust(len1))
            else:
                f.write(assign.ljust(len2))
            f.write(f"# {des}, range ({lwr}, {upr}), unit {uni}\n")
        elif re.match(r'`MPRco', line):
            newline = line[line.find("(")+1 : line.rfind(")")]
            par, val, uni, lwr, upr, des = re.split(r"\s*,\s*", newline, 5)
            par = par.strip()
            des = des.strip().strip('\"')
            assign = f"{par} = {val}"
            if len(assign) < len1:
                f.write(assign.ljust(len1))
            else:
                f.write(assign.ljust(len2))
            f.write(f"# {des}, range ({lwr}, {upr}), unit {uni}\n")
        elif re.match(r'`MPRcz', line):
            newline = line[line.find("(")+1 : line.rfind(")")]
            par, val, uni, des = re.split(r"\s*,\s*", newline, 3)
            par = par.strip()
            des = des.strip().strip('\"')
            assign = f"{par} = {val}"
            if len(assign) < len1:
                f.write(assign.ljust(len1))
            else:
                f.write(assign.ljust(len2))
            f.write(f"# {des}, unit {uni}\n")
        elif re.match(r'`MPRoz', line):
            newline = line[line.find("(")+1 : line.rfind(")")]
            split = re.split(r"\s*,\s*", newline)
            if len(split) > 4:
                print(split)
            par = split[0]
            if '(' in split[1]:
                val = split[1] + split[2]
            # par, val, uni, des = [
            #     x.strip().strip('\"').strip('\`')
            #     for x in split
            # ]            
            assign = f"{par} = {val}"
            if len(assign) < len1:
                f.write(assign.ljust(len1))
            else:
                f.write(assign.ljust(len2))
            f.write(f"# {des}, unit {uni}\n")
        elif re.match(r'`MPIcc', line):
            newline = line[line.find("(")+1 : line.rfind(")")]
            split = re.split(r"\s*,\s*", newline, 5)
            par, val, uni, lwr, upr, des = [
                x.strip().strip('\"').strip('\`')
                for x in split
            ]
            assign = f"{par} = {val}"
            if len(assign) < len1:
                f.write(assign.ljust(len1))
            else:
                f.write(assign.ljust(len2))
            f.write(f"# {des}, range ({lwr}, {upr}), unit {uni}\n")
        elif re.match(r'//', line):
            f.write(line.replace("//", "#"))   
        else:
            f.write(line)
