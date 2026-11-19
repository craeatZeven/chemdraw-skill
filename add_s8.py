import io,sys,json
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding="utf-8")
p=r"assets/parity/p22-four-columns-wide.json"
d=json.load(open(p,encoding="utf-8"))
h=d["header"]
h["substrates"]=[h["substrates"][0], h["substrates"][1],
                 {"smiles":"S1SSSSSSS1","label":"S8"}]
h["above"]=[t for t in h["above"] if "S8" not in t and "S₈" not in t]
h.setdefault("above",[]).append("DIPEA (3.0 equiv)")
json.dump(d,open(p,"w",encoding="utf-8"),ensure_ascii=False,indent=1)
print("表头反应物:",[s.get("label") for s in h["substrates"]],"| 条件:",h["above"],h["below"])
