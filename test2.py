import re

with open("code.txt", 'r') as f:
    lines = f.readlines()

with open("code_out.jl", 'w') as f:
    for line in lines:
        new_line = line.replace("`","")
        new_line = new_line.replace("//","#")
        new_line = new_line.replace(";","")
        new_line = new_line.replace("begin","")
        new_line = new_line.replace("end else","else")
        new_line = new_line.replace("else if","elseif")
        new_line = new_line.replace("$finish(0)","")
        new_line = new_line.replace("$strobe","println")        
        f.write("\t" + new_line)

#        if re.match(r'//', line):
#            f.write(new_line.replace("//", "#") + "\n")   
#        else:
#            f.write(new_line.replace("`","") + "\n")
