import ROOT
import pickle
import sys

sigs = pickle.load(file(sys.argv[1]))

useSignalsKey = None

srs = [ ]
mglus = set()
mlsps = set()
for nj in sigs:
    for lt in sigs[nj]:
        for ht in sigs[nj][lt]:
            sr = ( nj, lt, ht )
            assert not sr in srs
            srs.append(sr)
            subdict = sigs[nj][lt][ht]
            if useSignalsKey==None:
                useSignalsKey = "signals" in subdict
                print "Set useSignalsKey to ",useSignalsKey
            if useSignalsKey:
                subdict = subdict["signals"]
            for mg in subdict:
                if not mg in mglus:
                    mglus.add(mg)
                for ml in subdict[mg]:
                    if not ml in mlsps:
                        mlsps.add(ml)

