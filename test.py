import re

with open("param.txt", 'r') as f:
    lines = f.readlines()

len1, len2 = 25, 30

with open("out_bk.jl", 'w') as f:
    for line in lines:
        if (re.match(r'`MPRnb', line) or 
            re.match(r'`MPRex', line) or 
            re.match(r'`MPRcc', line) or
            re.match(r'`MPRoo', line) or 
            re.match(r'`MPRco', line) or
            re.match(r'`MPRoc', line) or 
            re.match(r'`MPRcz', line) or
            re.match(r'`MPRoz', line) or 
            re.match(r'`IPRnb', line) or
            re.match(r'`IPRex', line) or 
            re.match(r'`IPRcc', line) or
            re.match(r'`IPRoo', line) or 
            re.match(r'`IPRco', line) or
            re.match(r'`IPRoc', line) or 
            re.match(r'`IPRcz', line) or
            re.match(r'`IPRoz', line) or 
            re.match(r'`BPRco', line) or
            re.match(r'`BPRoz', line) or 
            re.match(r'`BPRcz', line) or
            re.match(r'`BPRnb', line)):
            paran_1st = line.find("(")
            comma_1st = line.find(",")
            comma_2nd = line.find(",", comma_1st+1)
            par = line[paran_1st+1:comma_1st].strip()
            val = line[comma_1st+1:comma_2nd].strip()
            f.write(f"{par}::Float64 = get(param, \"{par}\", {val})\n")
        elif (re.match(r'`MPInb', line) or 
            re.match(r'`MPIex', line) or 
            re.match(r'`MPIcc', line) or
            re.match(r'`MPIoo', line) or 
            re.match(r'`MPIco', line) or
            re.match(r'`MPIoc', line) or 
            re.match(r'`MPIcz', line) or
            re.match(r'`MPIoz', line) or 
            re.match(r'`MPIsw', line) or
            re.match(r'`MPIty', line) or 
            re.match(r'`IPInb', line) or
            re.match(r'`IPIex', line) or 
            re.match(r'`IPIcc', line) or
            re.match(r'`IPIoo', line) or 
            re.match(r'`IPIco', line) or
            re.match(r'`IPIoc', line) or 
            re.match(r'`IPIcz', line) or
            re.match(r'`IPIoz', line) or 
            re.match(r'`BPIcc', line) or
            re.match(r'`BPInb', line)):
            paran_1st = line.find("(")
            comma_1st = line.find(",")
            comma_2nd = line.find(",", comma_1st+1)
            par = line[paran_1st+1:comma_1st].strip()
            val = line[comma_1st+1:comma_2nd].strip()
            f.write(f"{par}::UInt8 = get(param, \"{par}\", {val})\n")
        elif re.match(r'//', line):
            f.write(line.replace("//", "#"))   
        else:
            f.write(line)
