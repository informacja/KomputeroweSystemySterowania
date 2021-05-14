# Program Grzanie.m

Informacje o pomieszczeniu:
 - 120 m2  powierzchni
 - 2,8 m wysokość

  Wentylacja 30% wymiany ciepła,
  druga strata ciepła to dyfuzja przez scianę 
  
  Jednostki:
 - energia kiloJoule, 
 - godzina jest jednostą czasu
 - obj w m3
  
Przebieg symulacji regulacji, zima mroźna do -15 *C 
![[./docs/Pasted image 20210514100906.png]]
![[./docs/Pasted image 20210514100954.png]]

Model wymennika ciepła 
![[docs/Pasted image 20210514101043.png]]

Charakterystyki
![[docs/Pasted image 20210514101232.png]]

> Symulacja grzania
> ![[docs/Pasted image 20210514101808.png]]

>
>![[docs/Pasted image 20210514102005.png]]

> Regulacja PID przepływu powietrza
> ![](https://i.imgur.com/w4rE5JA.gif)


> Regulacja PID 2 (zakłucenie nie odpowada ZeroOverHold (założenia Wienera))
> ![[docs/7FTDjT5Ftv.gif]]

> Synteza Bezpośrednia
> ![[docs/mW87XO2Vj1.gif]] 

> Przeregulowanie parametryczne 5%
> 12,15 profil wartości zadanej (lepsze od głowicy termostatycznej)
> ![[docs/V2eiJPBPZA.gif]]


> Reguloator przyrostowy/pozycyjny
> ![[docs/TXi01Oljzp.gif]]

> Regulator optymalny ( dewastuje układ wykonawczy)
> ![[docs/c1XV2AW0uL.gif]]
> układ jest niewłaściwy,
> ![[docs/zdygsSFYmc.gif]]

> Synteza bezpośrednia (ładnie reaguje)
> ![[docs/TEQxi7edy7.gif]]

> Ukłąd z przeregulowanie (tablicowe z podręcznika)
> ![[docs/Kp6DV7vy82.gif]]
> Model analityczny i regresyjny
> ![[docs/cAFZHkT2PD.gif]]

W optymalnej syntezie liczy się stosunek opóźnienia do stałej czasowej.
PS Zakłócenia są identyczne, raz generowane by można porównywać.

i - intuicyjna
o - optymalna stała czasowa/ del(opóźnienie)
b - bespośrednia
p - przeregulowanie

> Analityczna aproksymacja
> ![](https://i.imgur.com/PxVytF1.gif)

> Model regresyjny z tą samą małą wagą.
> ![](https://i.imgur.com/0ZGgRWR.gif)
> horyzont  liczenia  opytmalizacji (10) 4,5 powinno wystarczyć, wydłuża czas liczenia

> Liczenie kar za przyrosty sterowań, ważne bo pomiarowy
![[docs/xrqvdiUfkn.gif]]

80 wartości sterowania i 80 odpowiedzi wartości skokowej, dobrze by było używać zmiennoprzecinkoe, w automatyce można zastosować stało przecinkowy regulator PID

![[docs/Pasted image 20210514110713.png]]

wzór funkcji (sgnalizacyjny)
przepis na obiekt (obserwowalny)

Profil temperaturowy
![[docs/MATLAB_mLJJEISVzy.png]]