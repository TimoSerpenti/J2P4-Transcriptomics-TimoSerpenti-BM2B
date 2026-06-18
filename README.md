## 📁 Inhoud/structuur

- `Assets/` - Opmaak
- `Bronnen/` - Gebruikte bronnen 
- `Data/Raw/` – Onbewerkte sample data
- `Data/Verwerkt` - verwerkte datasets gegenereerd met scripts 
- `Data_stewardship/` - Voor de competentie beheren ga je aantonen dat je projectgegevens kunt beheren met behulp van GitHub. In deze folder kan je hulpvragen terugvinden om je op gang te helpen met de uitleg van data stewardship. 
- `Resultaten/` - Resultaten
- `Scripts/` – Gebruikte script om analyzes uit te voeren
- `.gitattributes` - Voor het uploaden van grote bestanden
- `.gitignore` - Bestanden negeren
- `README.md` - Tekstbestand

---

<p align="center">
  <img src="Assets/Voorblad.png" alt="Flowchart" />
</p>

# Genexpressie verschil tussen gezonde en patiënten met Reumatoïde Artritis

## Inleiding

Reumatoïde Artritis (RA) is de meest voorkomende chronische autoimmuunziekte van gewrichten. RA is een systemische auto-immuunziekte dit betekent dat het immuunsysteem het eigen lichaam aanvalt, niet op één specifieke plek maar verspreidt over het hele lichaam. Het veroorzaakt blijvende inflammatie wat zorgt voor zwellingen van gewrichten, vervorming, een verminderde dagelijkse functionaliteit en levenskwaliteit [(Sharif et al., 2018)](Bronnen/Clinical%20Anatomy%20-%202017%20-%20Sharif%20-%20Rheumatoid%20arthritis%20in%20review%20%20Clinical%20%20anatomical%20%20cellular%20and%20molecular%20points%20of.pdf). RA beïnvloed vooral de gewrichten en kan ook organen beïnvloeden, wat kan leiden tot permanente schade en beperking [(Bullock et al., 2019)](Bronnen/Rheumatoid%20Arthritis%20A%20Brief%20Overview.pdf).

De exacte oorzaak van de aandoening is nog niet bekend, zowel genetica als milieu dragen bij aan de ontwikkeling van de ziekte [(Aho & Heliövaara, 2004)](Bronnen/Risk%20factors%20for%20rheumatoid%20arthritis.pdf). Er is tot op heden nog geen geneesmiddel voor de aandoening, wat een aanleiding is voor het uitvoeren van het onderzoek [(Bullock et al., 2019)](Bronnen/Rheumatoid%20Arthritis%20A%20Brief%20Overview.pdf). 

Het doel van dit onderzoek is om te bepalen of er een verschil is in genexpressie tussen gezonde patiënten en patiënten met Reumatoïde Artritis. Dit doel wordt bereikt door te kijken naar welke genen differentieel tot expressie gebracht zijn bij RA-patiënten, door een DESeq2 analyse uit te voeren. Ook wordt er met behulp van een Volcanoplot gekeken welke genen de grootste verschillen in expressie tussen patiënten vertonen. Als laatst wordt er gekeken naar welke pathways betrokken zijn bij ontwikkeling van RA, door een pathway analyse uit te voeren.


## Methoden

Voor dit onderzoek is data gebruikt uit een eerdere studie van [(Platzer et al., 2019)](Bronnen/samples.pdf), deze data is verkregen van 4 samples van gezonde patiënten en 4 van patiënten met RA te zien in *Tabel 1*. RNA is geïsoleerd uit gewrichtsslijmvlies weefsel van patiënten met RA. Alle analyses zijn uitgevoerd met R (versie 4.5.3) in RStudio. Er is gebruik gemaakt van een script om reproduceerbaarheid te garanderen. De methode is te zien in de flowchart hieronder.
<p align="center">
  <img src="Assets/Flowchart.png" alt="Flowchart" width="400"/>
</p>

Alle packages die zijn gebruikt voor het uitvoeren van het script zijn Rsubread(2.24.0) <a href="Bronnen/RSubread.pdf"> (Liao et al., 2019) </a>, Rsamtools(2.26.0) <a href="Bronnen/RSamtools.pdf"> (Morgan et al., 2013) </a>, DESeq2(1.50.2) <a href="Bronnen/DESeq2.pdf"> (Love et al., 2014) </a>, EnhancedVolcano(1.28.2) <a href="Bronnen/EnhancedVolcano.pdf"> (O’Connell, 2025) </a>, goseq(1.62.0) <a href="Bronnen/goseq.pdf"> (young et al., 2012) </a>, geneLenDataBase(1.46.0) <a href="Bronnen/goseq.pdf"> (young et al., 2012) </a>, org.Hs.eg.db(3.22.0) <a href="Bronnen/OrganismDbi.pdf"> (Carlson et al., 2015) </a>, AnnotationDbi(1.72.0) <a href="Bronnen/AnnotationDBI.pdf"> (Pages et al., 2017)  </a>, GO.db(3.22.0) <a href="Bronnen/AnnotationDBI.pdf"> (Pages et al., 2017)  </a>, pathview(1.50.0) <a href="Bronnen/Pathview.pdf"> (Luo & Brouwer, 2013) </a>, dplyr(1.2.1) <a href="Bronnen/dplyr.pdf"> (Broatch et al., 2019) </a> en KEGGREST(1.50.0) <a href="Bronnen/KEGGREST.pdf"> (Tenenbaum et al., 2019)</a>.

**Tabel 1:** Data samples. In de tabel zijn de accension nummer, de leeftijden, het geslacht en of de patiënten reuma hebben te zien.

| Accession nummer | Leeftijd | Geslacht | Diagnose |
| :--- | :---: | :---: | :--- |
| SRR4785819 | 31 | female | Normal |
| SRR4785820 | 15 | female | Normal |
| SRR4785828 | 31 | female | Normal |
| SRR4785831 | 42 | female | Normal |
| SRR4785979 | 54 | female | Rheumatoid arthritis (established) |
| SRR4785980 | 66 | female | Rheumatoid arthritis (established) |
| SRR4785986 | 60 | female | Rheumatoid arthritis (established) |
| SRR4785988 | 59 | female | Rheumatoid arthritis (established) |

Eerst is er een index van het referentiegenoom gemaakt, om uitlijnen sneller te maken. Reads werden uitgelijnd tegen het humane referentiegenoom GRCh38.p14 met Rsubread (2.24.0), waarna BAM-bestanden gesorteerd zen geïndexeerd zijn met Rsamtools (2.26.0). Om te bepalen hoe vaak een read van een bepaald gen voorkomt. Is er vervolgens een count matrix gemaakt met featureCounts op basis van de bijbehorende GTF-annotatie.

Differentiële genexpressieanalyse werd uitgevoerd met het package DESeq2 (1.50.2).. Hiermee is vervolgens getest op statistisch significante verschillen tussen samples. Genen werden als significant beschouwd wanneer de p-waarde kleiner was dan 0,05 en wanneer de log2 fold change groter was dan 1 of -1, wat betekent dat er sprake is van duidelijke op- of neerregulatie.

Om de significante genen te visualiseren is een Volcano plot gemaakt met behulp van de EnhancedVolcano (1.28.2) package. Voor het vinden van een pathway met veranderde expressie is er een GO analyse uitgevoerd met goseq (1.62.0) en hulp van geneLenDataBase (1.46.0), org.Hs.eg.db (3.22.0) en AnnotationDbi (1.72.0).

Met GO.db zijn de Go termen en beschrijvingen verkregen. Als laatst is er een pathway analyse gedaan voor de TNF signaling pathway met behulp van de package pathview(1.50.0) en KEGGREST(1.50.0).

## Resultaten

Om inzicht te krijgen in eigenschappen van de te gebruiken spreuken is er een overzicht gemaakt, te vinden in [deze tabel](resultaten/top_10_spells.csv). Onvergeeflijke vloeken zijn niet meegenomen in dit overzicht. 

Om een afweging te maken welke spreuken het meest effectief zijn, is er onderzocht of er een verband te vinden is tussen kracht en accuraatheid. In [het resultaat hiervan](resultaten/spell_power_vs_accuracy.png) is te zien dat er een negatieve daling lijkt te zijn in kracht als de accuraatheid toeneemt. Een uitschieter is de onvergeeflijke vloek *Avada Kedavra*, met zowel hoge kracht als accuraatheid. 

## Conclusie

Spreuken met meer accuraatheid lijken minder krachtig te zijn. Een uitzondering op deze trend is de onvergeeflijke vloek *Avada Kedavra*, welke beter niet gebruikt kan worden. 

