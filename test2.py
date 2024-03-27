
combi = ["ELU 1.35 G", "ELU 1.35 G + 1.5 S"]
dictcombi = {}
for name in combi:
    flexion = "t"
    traction = "f"
    dictcombi[name] = {"flexion": flexion, "traction": traction}

print(dictcombi)
flex = dictcombi["ELU 1.35 G"]["flexion"]