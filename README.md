# Genexpressie verschil tussen gezonde en patiënten met Reumatoïde Artritis

## 📁 Inhoud/structuur

- `data/raw/` – fictionele datasets voor de analyse van spreuk effectiviteit, gevaar en welke spreuken het beste samengaan met verschillende types staf.  
- `data/processed` - verwerkte datasets gegenereerd met scripts 
- `scripts/` – scripts om prachtige onzin te genereren
- `resultaten/` - grafieken en tabellen
- `bronnen/` - gebruikte bronnen 
- `README.md` - het document om de tekst hier te genereren
- `assets/` - overige documenten voor de opmaak van deze pagina
- `data_stewardship/` - Voor de competentie beheren ga je aantonen dat je projectgegevens kunt beheren met behulp van GitHub. In deze folder kan je hulpvragen terugvinden om je op gang te helpen met de uitleg van data stewardship. 

---

## Inleiding

Reumatoïde Artritis (RA) is de meest voorkomende chronische autoimmuunziekte van gewrichten. RA is een systemische auto-immuunziekte dit betekent dat het immuunsysteem het eigen lichaam aanvalt, niet op één specifieke plek maar verspreidt over het hele lichaam. Het veroorzaakt blijvende inflammatie wat zorgt voor zwellingen van gewrichten, vervorming, een verminderde dagelijkse functionaliteit en levenskwaliteit [(Sharif et al., 2018)](Bronnen/Clinical Anatomy - 2017 - Sharif - Rheumatoid arthritis in review  Clinical  anatomical  cellular and molecular points of.pdf). RA beïnvloed vooral de gewrichten en kan ook organen beïnvloeden, wat kan leiden tot permanente schade en beperking [(Bullock et al., 2019)](Bronnen/Rheumatoid Arthritis A Brief Overview.pdf). De exacte oorzaak van de aandoening is nog niet bekend, zowel genetica als milieu dragen bij aan de ontwikkeling van de ziekte [(Aho & Heliövaara, 2004)](Bronnen/Risk factors for rheumatoid arthritis.pdf). Er is tot op heden nog geen geneesmiddel voor de aandoening, wat een aanleiding is voor het uitvoeren van het onderzoek [(Bullock et al., 2019)](Bronnen/Rheumatoid Arthritis A Brief Overview.pdf). Het doel van dit onderzoek is om te bepalen of er een verschil is in genexpressie tussen gezonde patiënten en patiënten met Reumatoïde Artritis. Dit doel wordt bereikt door te kijken naar welke genen differentieel tot expressie gebracht zijn bij RA-patiënten, door een DESeq2 analyse uit te voeren. Ook wordt er met behulp van een Volcanoplot gekeken welke genen de grootste verschillen in expressie tussen patiënten vertonen. Als laatst wordt er gekeken naar welke pathways betrokken zijn bij ontwikkeling van RA, door een pathway analyse uit te voeren.


## Methoden

De data is verstrekt door de Afdeling Magische Wetshandhaving en Ollivanders Wandwinkel Archieven. 

De ruwe data van spreuken is eerst bewerkt voor analyse met behulp van [scripts/01_clean_spell_data.R](scripts/01_clean_spell_data.R). Vervolgens zijn de spreuken geanalyseerd op kracht en nauwkeurigheid met [scripts/02_spell_analysis.R](scripts/02_spell_analysis.R).

## Resultaten

Om inzicht te krijgen in eigenschappen van de te gebruiken spreuken is er een overzicht gemaakt, te vinden in [deze tabel](resultaten/top_10_spells.csv). Onvergeeflijke vloeken zijn niet meegenomen in dit overzicht. 

Om een afweging te maken welke spreuken het meest effectief zijn, is er onderzocht of er een verband te vinden is tussen kracht en accuraatheid. In [het resultaat hiervan](resultaten/spell_power_vs_accuracy.png) is te zien dat er een negatieve daling lijkt te zijn in kracht als de accuraatheid toeneemt. Een uitschieter is de onvergeeflijke vloek *Avada Kedavra*, met zowel hoge kracht als accuraatheid. 

## Conclusie

Spreuken met meer accuraatheid lijken minder krachtig te zijn. Een uitzondering op deze trend is de onvergeeflijke vloek *Avada Kedavra*, welke beter niet gebruikt kan worden. 

