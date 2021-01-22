
"""
###TRABALHO FISICA NUCLEAR###
   
NOME: Lhwan Silva 
RA:182100
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#-----CONSTANTES OBTIDAS DO LIVRO DO WILLIAMS-----#

Av = 15.56  # constante de volume, em MeV
As = 17.23  # constante de superfície, em MeV 
Ac = .697   # constante de Coulomb, em MeV
Aa = 23.285 # constante de assimetria, em MeV
Ap = 12     # constante de paridade, em MeV
mh = 938.791  #*na verdade é mp*c^2*
mn = 939.573   #(*na verdade é mn*c^2*)

#-----FUNÇÕES MASSA SEMI-EMPÍRICA-----#

def delta(a,z):
  a=int(a)
  z=int(z)
  if (z%2) and ((a-z)%2):
    return -Ap/(a**.5)
  if (z%2 and not ((a-z)%2)) or (not (z%2)and ((a-z)%2)):
    return 0
  if (not (z%2)) and (not ((a-z)%2)):
    return Ap/(a**.5)

def be(a, z):
 aux = Av*a - As*(a**(2/3)) - Ac*(z**2)/(a**(1/3)) - Aa*((a - 2*z)**2)/a + delta(a, z)
 return aux

def sem(a, z):
  aux = z*mh + (a-z)*mn - be(a, z)
  return aux

#-----A=101 (COMPARANDO COM O WILLIAMS)-----#
isobar101=[[],[]]

# eixo X, numero Z
for i in range(7) :
    isobar101[0].append(41+i) 
# eixo y, massas semi empirica
for i in range(7) :  
    isobar101[1].append(sem(101, 41+i))

x101 = np.arange(41,47.1,.2)
for i in range(len(x101)):
  x101[i]=round(x101[i],2)

y101=[sem(101, i) for i in x101]

plt.plot(isobar101[0],isobar101[1],'ro')
plt.plot(x101,y101)
plt.ticklabel_format(useOffset=False)
plt.ylabel('Massa Semi-Empírica [MeV/c^2]')
plt.xlabel('Z')
plt.show()

#-----A=55-----#
iso55Ini=22
iso55Fim=28

isobar55=[[],[]]
# eixo X, numero Z
for i in range(7) :
    isobar55[0].append(iso55Ini+i) 
# eixo y, massas semi empirica
for i in range(7) :  
    isobar55[1].append(sem(55, iso55Ini+i))

x55 = np.arange(iso55Ini,iso55Fim+.2,.2)
for i in range(len(x55)):
  x55[i]=round(x55[i],2)

y55=[sem(55, i) for i in x55]

plt.plot(isobar55[0],isobar55[1],'ro')
plt.plot(x55,y55)
plt.ticklabel_format(useOffset=False)
plt.xlabel('Z')
plt.ylabel('Massa Semi-Empírica [MeV/c^2]')
plt.show()

#-----A=97-----#
iso97Ini=38.9
iso97Fim=46

isobar97=[[],[]]
# eixo X, numero Z
for i in range(int(iso97Fim-iso97Ini)+1) :
    isobar97[0].append(iso97Ini+i) 
# eixo y, massas semi empirica
for i in range(int(iso97Fim-iso97Ini)+1):  
    isobar97[1].append(sem(97, iso97Ini+i))

x97 = np.arange(iso97Ini,iso97Fim+.2,.2)
for i in range(len(x97)):
  x97[i]=round(x97[i],2)

y97=[sem(97, i) for i in x97]

plt.plot(isobar97[0],isobar97[1],'ro')
plt.plot(x97,y97)
plt.ticklabel_format(useOffset=False)
plt.xlabel('Z')
plt.ylabel('Massa Semi-Empírica [MeV/c^2]')
plt.show()

#-------DECAIMENTO DOS ISÓBAROS, A=101-------#

# meias vidas, em minutos, dos isótopos Nb, Mo e Tc #
mv=[7.1/60, 14.61, 14.22]

n0=1
w1=round(np.log(2)/mv[0],5)
w2=round(np.log(2)/mv[1],5)
w3=round(np.log(2)/mv[2],5)

# Equações de decaimento para o Nb, Mo e Tc com t em minutos #
# População do Nb
def n1(t): 
  return (n0*np.exp(-w1*t))
# População do Mo
def n2(t):
  return ((n0*w1*np.exp(-w1*t)-n0*w2*np.exp(-w2*t))/(w1 - w2))
# População do Nb
def n3(t):
  return ((n0*np.exp(-(w1 + w2 - w3)*t)*(-w2**2*(w1-w3)*np.exp((w1-w3)*t) + w1**2*(w2-w3)*np.exp((w2-w3)*t) + (w1-w2)*w3**2*np.exp((w1+w2-2*w3)*t)))/((w1-w2)*(w1-w3)*(w2-w3)))

tIni=0.2
tFim=0.6
tQuant=100
tPasso=(tFim-tIni)/tQuant

porc=100/n0
tempo=[(tIni+round(x*tPasso,2)) for x in range(tQuant+1)]
quant1=[(round(n1(tIni+round(x*tPasso,2)),6)*porc) for x in range(tQuant+1)]
quant2=[(round(n2(tIni+round(x*tPasso,2)),6)*porc) for x in range(tQuant+1)]
quant3=[(round(n3(tIni+round(x*tPasso,2)),6)*porc) for x in range(tQuant+1)]

plt.plot(tempo,quant1,'r',label='Nb')
plt.plot(tempo,quant2,'b',label='Mo')
plt.plot(tempo,quant3,'g',label='T c')
plt.plot(tempo,[10 for x in range(tQuant+1)],'k')
plt.ticklabel_format(useOffset=False)
plt.ylabel('% do elemento')
plt.xlabel('Tempo [min]')
plt.legend()
plt.show()
