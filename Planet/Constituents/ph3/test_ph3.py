import ph3_jh
import matplotlib.pyplot as plt

f = []
for i in range(500):
    f.append(float(i)/1.0)

P_dict = {'H2':0,'HE':1,'PH3':2}
P = 6.3
T = 135.0
X_ph3 = 4.3E-7
X_he = 0.12
X_h2 = 1.0 - (X_ph3+X_he)
X=[X_h2,X_he,X_ph3]
a = ph3_jh.alpha(f,T,P,X,P_dict)
#plt.plot(f,a)
plt.semilogy(f,a)

#Data from Icarus152
P = 3.15
T = 175.
X_ph3 = 0.022
X_he = 0.098
X_h2 = 0.88
X=[X_h2,X_he,X_ph3]
f = [1.51,2.25,8.3,13.3,21.6,27.0]
am= [2.15,3.92,17.54,21.13,21.03,16.44]
e = [0.25,0.47,0.81,4.66,4.31,14.7]
plt.figure(2)
plt.plot(f,am)
plt.errorbar(f,am,e)
ac = ph3_jh.alpha(f,T,P,X,P_dict)
plt.plot(f,ac)
