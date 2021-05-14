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
![[Pasted image 20210514100906.png]]
![[Pasted image 20210514100954.png]]

Model wymennika ciepła 
![[Pasted image 20210514101043.png]]

Charakterystyki
![[Pasted image 20210514101232.png]]

> Symulacja grzania
> ![[Pasted image 20210514101808.png]]

>
>![[Pasted image 20210514102005.png]]

> Regulacja PID przepływu powietrza
> ![](https://i.imgur.com/w4rE5JA.gif)


> Regulacja PID 2 (zakłucenie nie odpowada ZeroOverHold (założenia Wienera))
> ![[7FTDjT5Ftv.gif]]

> Synteza Bezpośrednia
> ![[mW87XO2Vj1.gif]] 

> Przeregulowanie parametryczne 5%
> 12,15 profil wartości zadanej (lepsze od głowicy termostatycznej)
> ![[V2eiJPBPZA.gif]]


> Reguloator przyrostowy/pozycyjny
> ![[TXi01Oljzp.gif]]

> Regulator optymalny ( dewastuje układ wykonawczy)
> ![[c1XV2AW0uL.gif]]
> układ jest niewłaściwy,
> ![[zdygsSFYmc.gif]]

> Synteza bezpośrednia (ładnie reaguje)
> ![[TEQxi7edy7.gif]]

> Ukłąd z przeregulowanie (tablicowe z podręcznika)
> ![[Kp6DV7vy82.gif]]
> Model analityczny i regresyjny
> ![[cAFZHkT2PD.gif]]

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
![[xrqvdiUfkn.gif]]

80 wartości sterowania i 80 odpowiedzi wartości skokowej, dobrze by było używać zmiennoprzecinkoe, w automatyce można zastosować stało przecinkowy regulator PID

![[Pasted image 20210514110713.png]]

wzór funkcji (sgnalizacyjny)
przepis na obiekt (obserwowalny)

Profil temperaturowy
![[MATLAB_mLJJEISVzy.png]]