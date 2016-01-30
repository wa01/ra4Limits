import os,sys
from fnmatch import fnmatch
import pickle

def parseLogFile(fn,smglu,smlsp):
    for i,l in enumerate(open(fn)):
        if i>1:
            print "*** log file parsing failed for ",fn
            return None
        idx = l.find("{")
        if idx<0:
            print "*** log file parsing failed for ",fn
            return None
        hdr = l[:idx]
        lims = eval(l[idx:-1])
        fields = hdr.split()
        if not ( len(fields)==7 and fields[0]=="Result" ):
            print "*** log file parsing failed for ",fn
            return None
        if not ( fields[4]==smglu and fields[5]==smlsp ):
            print "*** log file parsing failed for ",fn
            return None
    return [ int(fields[4]), int(fields[5]), lims ]

assert len(sys.argv)==3

indir = sys.argv[1]
assert os.path.isdir(indir)

out = sys.argv[2]
assert not os.path.exists(out)

results = { }
nsucc = 0
nfail = 0
for f in os.listdir(indir):
    if not fnmatch(f,"limit_[0-9]*_[0-9]*.log"):
        continue
    ffields = f.split(".")[0].split("_")
    assert len(ffields)==3
    assert ffields[1].isdigit() and ffields[2].isdigit()
    parseRes = parseLogFile(indir+"/"+f,ffields[1],ffields[2])
    if parseRes==None:
        nfail += 1
        continue
    mglu,mlsp,lims = parseRes
    if not mglu in results:
        results[mglu] = { }
    assert not ( mlsp in results[mglu])
    results[mglu][mlsp] = lims
    nsucc += 1

print "Successfully parsed",nsucc,"files (",nfail,"files failed)"

out = sys.argv[2]
assert not os.path.exists(out)
fout = open(out,"wb")
pickle.dump(results,fout)
fout.close()

