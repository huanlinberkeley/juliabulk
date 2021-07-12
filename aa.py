import re
input = """`MPIcc(BULKMOD, 0, "", 0, 2, "0: SOI multi-gate; 1: Bulk multi-gate; 2: for decoupled bulk multi-gate")"""
bb = "0: SOI multi-gate; 1: Bulk multi-gate; 2: for decoupled bulk multi-gate"
aa = re.findall(r'\"(.*?)\"', input)
print(input.rfind(bb))
print(input[30:101])