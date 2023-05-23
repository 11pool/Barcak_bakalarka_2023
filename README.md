# Barcak_bakalarka_2023
Tento repozitár obsahuje kód a dáta k bakalárskej práci: Určovanie počtu zhlukov v dátach, Adam Barčák, Univerzita Komenského FMFI, 2023<br>
Tento readme obsahuje popis štruktúri kódu a najmä štruktúri a pomenovania dát

Kód aj dáta sú v programovacom jazyku R.

upratovanie kódu in progress. Hotové bude 23.5. pred polnocou najneskôr

############<br>
**Štruktúra kódu:**<br>
R scripty sú písané ako zoznam funkcií, ktoré sú spúšťané manuálne podľa potreby, nie ako script ktorý sa spustí celý.<br>
Script generate.R obsahuje funkcie pre generovanie modelových dát.<br>
Script stat.R obsahuje funkcie pre realizáciu štatistiky na vygenerovaných dátach.<br>
Konkrétnejšie informácie sú uvedené ako komenty v kóde.<br>
###########<br>

###########<br>
**Štruktúra dát:**<br>
Štatistika prebiehala na 1000 realizáciach každej modelovej situácie, okrem manuálnych ako lakťový graf a siluetový diagram ktoré prebehli na realizácií prvých 50 dát.<br>
Voľne povedané, premenné majú nasledovné názvoslovie: typdát_situácia_špecifikácianavyše

x_situácia:
+ realizácia vygenerovaných modelových dát ktoré boli použité v bakalárskej práci
+ list dĺžky 1000 ktorého zložky sú matice, kde každá matica reprezentuje jednu realizáciu modelovej situácie a každý riadok reprezentuje vektor
+ pamäťovo najväčšie
              
stat_situácia:
+ výsledky štatistiky pre metódy zhlukovania dát
+ list dĺžky 5 ktorého zložky sú vektory dĺžky 1000, kde prvá zložka listu je hodnotenie metódy kmeans, druhá kmedoid, tretia MClust, štvrtá DBScan a každý člen vektoru reprezentuje zlomok                        správne zaradených bodov pre danú realizáciu modelových dát. Piata zložka listu je vektor najlepšieho DBScan epsilon parametru pre danú realizáciu dát.
                 
stat_urcenie_situácia:
+ výsledky štatistiky pre metódy určovania počtu zhlukov v dátach
+ list DOKONCI

názvy modelových situácií v kóde:<br>
+ Benchmark -> bench
+ Málo zhluku -> lownoise
+ Stredne zhluku -> midnoise
+ Veľa zhluku -> highnoise
+ T rozdelenie -> distr_t
+ Rovnomerné rozdelenie -> distr_unif
+ Variabilné kovariancie -> covar
+ Blízke zhluky -> near
+ Vnorené zhluky -> nearnear
+ Rôzne veľké zhluky -> weight<br>
###########
