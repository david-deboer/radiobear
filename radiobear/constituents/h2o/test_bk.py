import h2o_bk
import h2o_ddb
import rdh2o as r
otherpar=[]
iFig = 0
setNo = 2
for setNo in range(len(r.setInd)):
    fmin = min(r.dset[setNo][1][0])
    fmin*=0.9
    fmax = max(r.dset[setNo][1][0])
    fmax*=1.1
    fstep = (fmax-fmin)/100.0
    f = []
    for i in range(100):
        f.append(fmin + i*fstep)

    P = r.dset[setNo][0][0]
    T = r.dset[setNo][0][1]
    X_h2o = r.dset[setNo][0][2]
    X_he = r.dset[setNo][0][3]
    X_h2 = r.dset[setNo][0][4]
    X = [X_h2,X_he,X_h2o]
    P_dict = {'H2':0,'HE':1,'H2O':2}
    a = h2o_bk.alpha(f,T,P,X,P_dict,otherpar)
    b = h2o_ddb.alpha(f,T,P,X,P_dict,otherpar)
    r.plt.figure(iFig)
    iFig+=1
    r.plt.plot(r.dset[setNo][1][0],r.dset[setNo][1][1],'*')
    r.plt.errorbar(r.dset[setNo][1][0],r.dset[setNo][1][1],r.dset[setNo][1][2],fmt='*')
    r.plt.plot(f,a,'r')
    r.plt.plot(f,b,'b')
    r.plt.title(str(r.dset[setNo][0]))
