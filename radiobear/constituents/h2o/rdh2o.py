import matplotlib.pyplot as plt

fp = open('h2o_data-1.csv','r')

h2o = []
columns = []
i=0
for line in fp:
    dataStr = line.split(',')
    if i==0:
        columns = dataStr
    else:
        data = []
        for v in dataStr:
            v=v.strip()
            data.append(float(v))
        h2o.append(data)
    i+=1
print str(i)+' lines of data'

cols={}

i=0
for v in columns:
    v = v.strip()
    cols[v]=i
    i+=1

f = []
P = []
T = []
Xh2o = []
Xh2 = []
Xhe = []
ah2o = []
aerr = []
for d in h2o:
    f.append(d[cols['Frequency (Hz)']]/1.0E9)
    P.append(d[cols['Pressure (bars)']])
    T.append(d[cols['Temperature K']])
    Xh2o.append(d[cols['X H2O']])
    Xh2.append(d[cols['X H2 mole fraction']])
    Xhe.append(d[cols['X He mole fraction']])
    ah2o.append(d[cols['Absorption (dB/km)']])
    aerr.append(d[cols['2_ Absorption']])

dset = []
print 'Find sets of measurements specified by P, T, Xh2o, Xh2, Xhe'
def approx(a,b,c=0.01):
    m = abs(a-b)/a
    if m > c:
        return False
    else:
        return True
def findSet(PSetVal,TSetVal,Xh2oSetVal,XheSetVal,Xh2SetVal):
    print 'PsetVal:  ',PSetVal
    print 'TsetVal:  ',TSetVal
    print 'Xh2osetVal:  ',Xh2oSetVal
    print 'XhesetVal:  ',XheSetVal
    print 'Xh2setVal:  ',Xh2SetVal
    dentry = [[PSetVal,TSetVal,Xh2oSetVal,XheSetVal,Xh2SetVal]]
    j = 0
    fSet = []
    aSet = []
    eSet = []
    for i,pval in enumerate(P):
        if approx(pval,PSetVal) and approx(T[i],TSetVal) and approx(Xh2o[i],Xh2oSetVal) and approx(Xhe[i],XheSetVal) and approx(Xh2[i],Xh2SetVal):
            j+=1
            fSet.append(f[i])
            aSet.append(ah2o[i])
            eSet.append(aerr[i])
    val = [fSet,aSet,eSet]
    dentry.append(val)
    print 'Found '+str(j)+' with  ',dentry[0]
    return dentry

# index numbers for different measurement sets
setInd = [7, 17, 33, 43, 52, 61, 71, 81, 98, 107, 116, 125, 135, 153, 161, 327]
for si in setInd:
    dentry = findSet(PSetVal = P[si], TSetVal = T[si], Xh2oSetVal = Xh2o[si], XheSetVal = Xhe[si], Xh2SetVal = Xh2[si])
    dset.append(dentry)
    #plt.plot(dentry[1][0],dentry[1][1],'*')
    #plt.errorbar(dentry[1][0],dentry[1][1],dentry[1][2],fmt='*')


