import os,sys
from fnmatch import fnmatch
import pickle

def parseLogFile(fn,smglu,smlsp):
    for i,l in enumerate(open(fn)):
        assert i<2
        idx = l.find("{")
        assert idx>=0
        hdr = l[:idx]
        lims = eval(l[idx:-1])
        fields = hdr.split()
        assert len(fields)==7 and fields[0]=="Result"
        assert fields[4]==smglu and fields[5]==smlsp
    return [ int(fields[4]), int(fields[5]), lims ]

assert len(sys.argv)==3

indir = sys.argv[1]
assert os.path.isdir(indir)

out = sys.argv[2]
assert not os.path.exists(out)

results = { }
for f in os.listdir(indir):
    if not fnmatch(f,"limit_[0-9]*_[0-9]*.log"):
        continue
    ffields = f.split(".")[0].split("_")
    assert len(ffields)==3
    assert ffields[1].isdigit() and ffields[2].isdigit()
    mglu,mlsp,lims = parseLogFile(indir+"/"+f,ffields[1],ffields[2])
    if not mglu in results:
        results[mglu] = { }
    assert not ( mlsp in results[mglu])
    results[mglu][mlsp] = lims

out = sys.argv[2]
assert not os.path.exists(out)
fout = open(out,"wb")
pickle.dump(results,fout)
fout.close()

