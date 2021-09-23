# stress bicycle map italy
Mappa dello stress da traffico in bicicletta in Italia- basato sui modelli definiti da [BikeOttawa](https://github.com/BikeOttawa)

Ne parliamo anche qui 

https://medium.com/tantotanto/lo-stress-dei-percorsi-ciclabili-ee7573ec8284

## Significato dei livelli

Il livello di stress da traffico (LST) è una categorizzazione applicata a un segmento di strada o un incrocio e indica lo stress da traffico che grava sui ciclisti. 

I livelli di stress da traffico vanno da 1 a 4 come segue:

* LST 1: separazione marcata da tutti i tipi di traffico, fatta eccezione per quello a bassa velocità e bassa intensità. Incroci semplici. Adatto ai bambini
* LST 2: Tranne che nei casi di traffico a bassa velocità/intensità, i ciclisti hanno uno spazio specifico per muoversi che evita loro le interazioni con il traffico, eccetto che negli incroci. C'è separazione fisica dal traffico a velocità più alta e a più corsie. Gli incroci sono facili da affrontare per un adulto. Corrisponde ai criteri per le ciclabili olandesi. Si tratta di un livello di stress da traffico che la maggior parte degli adulti riesce a tollerare, in particolare coloro che sono a volte classificati come "interessati, ma preoccupati"
* LST 3: Comporta l'interazione con traffico a velocità moderata o a più corsie, o vicinanza con traffico a velocità più alta. Un livello di stress da traffico accettabile da coloro classificati come "entusiasti e sicuri di sé"
* LST 4: Comporta l'interazione con traffico a velocità più alta o vicinanza con traffico a velocità elevata. Un livello di stress accettabile soltanto da coloro che sono classificati come "forti e impavidi"

Ci sono dei criteri che determinano l'LST per i tratti di strada e gli incroci. I criteri per l'LST sono stati pubblicati la prima volta nel 2012 in un report di Mekuria, Furth e Nixon, edito dal Mineta Transportation Institute.

I livelli di stress da traffico per un percorso vengono composti a partire dai suoi segmenti usando una logica di collegamento più debole. Ciò significa che se la maggior parte dei segmenti di un percorso hanno LST 1 o 2, ma uno o un piccolo numero di essi hanno LST 3, il percorso ha complessivamente LST 3.

## Codice sorgente originario
- mappa<br/>https://github.com/rcmc2020/stressmap/
- calcolo dello stress<br/>https://github.com/BikeItaly/stressmodel (fork italiano https://github.com/BikeOttawa/stressmodel)
- criteri con cui è stato fatto il calcolo<br/>https://peterfurth.sites.northeastern.edu/2014/05/21/criteria-for-level-of-traffic-stress/
