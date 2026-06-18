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

`Reumatoïde Artritis (RA)` is de meest voorkomende chronische autoimmuunziekte van gewrichten. RA is een systemische auto-immuunziekte dit betekent dat het immuunsysteem het eigen lichaam aanvalt, niet op één specifieke plek maar verspreidt over het hele lichaam. Het veroorzaakt blijvende inflammatie wat zorgt voor zwellingen van gewrichten, vervorming, een verminderde dagelijkse functionaliteit en levenskwaliteit [(Sharif et al., 2018)](Bronnen/Clinical%20Anatomy%20-%202017%20-%20Sharif%20-%20Rheumatoid%20arthritis%20in%20review%20%20Clinical%20%20anatomical%20%20cellular%20and%20molecular%20points%20of.pdf). RA beïnvloed vooral de gewrichten en kan ook organen beïnvloeden, wat kan leiden tot permanente schade en beperking [(Bullock et al., 2019)](Bronnen/Rheumatoid%20Arthritis%20A%20Brief%20Overview.pdf).

De exacte oorzaak van de aandoening is nog niet bekend, zowel `genetica` als `milieu` dragen bij aan de ontwikkeling van de ziekte [(Aho & Heliövaara, 2004)](Bronnen/Risk%20factors%20for%20rheumatoid%20arthritis.pdf). Er is tot op heden nog geen geneesmiddel voor de aandoening, wat een aanleiding is voor het uitvoeren van het onderzoek [(Bullock et al., 2019)](Bronnen/Rheumatoid%20Arthritis%20A%20Brief%20Overview.pdf). 

Het doel van dit onderzoek is om te bepalen of er een verschil is in `genexpressie` tussen gezonde patiënten en patiënten met Reumatoïde Artritis. Dit doel wordt bereikt door te kijken naar welke genen `differentieel` tot expressie gebracht zijn bij RA-patiënten, door een `DESeq2` analyse uit te voeren. Ook wordt er met behulp van een `Volcanoplot` gekeken welke genen de grootste verschillen in expressie tussen patiënten vertonen. Als laatst wordt er gekeken naar welke `pathways` betrokken zijn bij ontwikkeling van RA, door een pathway analyse uit te voeren.


## Methoden

Voor dit onderzoek is data gebruikt uit een eerdere studie van [(Platzer et al., 2019)](Bronnen/samples.pdf), deze data is verkregen van 4 samples van gezonde patiënten en 4 van patiënten met RA te zien in **Tabel 1**. RNA is geïsoleerd uit gewrichtsslijmvlies weefsel van patiënten met RA. Alle analyses zijn uitgevoerd met `R (versie 4.5.3)` in `RStudio`. Er is gebruik gemaakt van een script om reproduceerbaarheid te garanderen. De methode is te zien in de `flowchart` hieronder.
<p align="center">
  <img src="Assets/Flowchart.png" alt="Flowchart" width="400"/>
</p>

Alle packages die zijn gebruikt voor het uitvoeren van het script zijn `Rsubread(2.24.0)` <a href="Bronnen/RSubread.pdf"> (Liao et al., 2019) </a>, `Rsamtools(2.26.0)` <a href="Bronnen/RSamtools.pdf"> (Morgan et al., 2013) </a>, `DESeq2(1.50.2)` <a href="Bronnen/DESeq2.pdf"> (Love et al., 2014) </a>, `EnhancedVolcano(1.28.2)` <a href="Bronnen/EnhancedVolcano.pdf"> (O’Connell, 2025) </a>, `goseq(1.62.0)` <a href="Bronnen/goseq.pdf"> (young et al., 2012) </a>, `geneLenDataBase(1.46.0)` <a href="Bronnen/goseq.pdf"> (young et al., 2012) </a>, `org.Hs.eg.db(3.22.0)` <a href="Bronnen/OrganismDbi.pdf"> (Carlson et al., 2015) </a>, `AnnotationDbi(1.72.0)` <a href="Bronnen/AnnotationDBI.pdf"> (Pages et al., 2017)  </a>, `GO.db(3.22.0)` <a href="Bronnen/AnnotationDBI.pdf"> (Pages et al., 2017)  </a>, `pathview(1.50.0)` <a href="Bronnen/Pathview.pdf"> (Luo & Brouwer, 2013) </a>, `ggplot2(4.0.3)`  <a href="Bronnen/ggplot2.pdf"> (Wickham, 2011) </a>, `dplyr(1.2.1)` <a href="Bronnen/dplyr.pdf"> (Broatch et al., 2019) </a> en `KEGGREST(1.50.0)` <a href="Bronnen/KEGGREST.pdf"> (Tenenbaum et al., 2019)</a>.

**Tabel 1: Data samples** In de tabel zijn de `accension nummers`, `de leeftijden`, `het geslacht` en of de `patiënten reuma` hebben te zien.

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

Eerst is er een `index` van het `referentiegenoom` gemaakt, om uitlijnen sneller te maken. Reads werden uitgelijnd tegen het humane referentiegenoom `GRCh38.p14` met `Rsubread (2.24.0)`, waarna `BAM-bestanden` gesorteerd zen geïndexeerd zijn met `Rsamtools (2.26.0)`. Om te bepalen hoe vaak een read van een bepaald gen voorkomt. Is er vervolgens een `count matrix` gemaakt met `featureCounts` op basis van de bijbehorende `GTF-annotatie`.

`Differentiële genexpressieanalyse` werd uitgevoerd met het package `DESeq2 (1.50.2)`. Hiermee is vervolgens getest op statistisch significante verschillen tussen samples. Genen werden als significant beschouwd wanneer de `p-waarde` kleiner was dan `0,05` en wanneer de `log2 fold change` groter was dan `1` of `-1`, wat betekent dat er sprake is van duidelijke op- of neerregulatie.

Om de significante genen te visualiseren is een `Volcano plot` gemaakt met behulp van de `EnhancedVolcano (1.28.2)` package. Voor het vinden van een pathway met veranderde expressie is er een `GO analyse` uitgevoerd met `goseq (1.62.0)` en hulp van `geneLenDataBase (1.46.0)`, `org.Hs.eg.db (3.22.0)` en `AnnotationDbi (1.72.0)`.

Met `GO.db` zijn de Go termen en beschrijvingen verkregen. Als laatst is er een pathway analyse gedaan voor de `TNF signaling pathway` met behulp van de package `pathview(1.50.0)` en `KEGGREST(1.50.0)`.

## Resultaten

## Differentiele expressie van genen

De RNA-seq analyse heeft `2085` significant op gereguleerde genen gevonden bij patiënten met Reuma in vergelijking met de gezonde controles, en `2487` genen die significant neer gereguleerd zijn. In totaal werden er `4572` significant differentieel tot expressie gebracht, van de totaal `29407` variabele genen die zijn meegenomen in de analyse.

## Visualisatie genexpressie

In de `volcano plot` in **Figuur 1** is te zien er is ook een [tweede volcano plot]( Resultaten/VolcanoplotREUMA2.png) gemaakt waar meer genen in te zien zijn doordat er minder strikte significantie gebruikt is.
<p align="center">
  <img src="Resultaten/VolcanoplotREUMA.png" alt="VolcanoPlot" width="400"/>
</p>
**Figuur 1: VolcanoPlot** Op de  x-as is de `log2foldchange` te zien en op de y-as de `p-waarde`. Genen zijn gefilterd op basis van een gecorrigeerde p-waarde `($p_{\text{adj}}$) $< 0.05$` en een absolute `$\log_2 \text{Fold Change } (\text{FC}) > 1$ (of $< -1$)`. Genen die aan beide criteria voldoen, zijn in de plot `rood` gekleurd. Genen die wel een sterke expressieverandering ondergingen maar niet statistisch significant bleken, zijn `groen` gemarkeerd, terwijl niet-significante genen `blauw of grijs` zijn weergegeven.
De volcano plots lieten een duidelijke spreiding van genen zien die significant omhoog of omlaag gereguleerd zijn. De plots toonde bovendien enkele zeer sterke uitschieters. Zoals `SRGN`, `ANKRD30BL` en `MT-ND6 `op door een erg hoge statistische significantie.

## GO-termen

In **Figuur 2** is een staafdiagram te zien waarin de meest verrijkte `GO-termen` te zien zijn. Er werd aangetoond dat immune response, adaptive immune response, immune system process, immunoglobulin complex en antigen binding een oververtegenwoordiging hadden bij mensen met Reuma. 
<p align="center">
  <img src="Resultaten/Rplot.png" alt="GO-termen" width="400"/>
</p>
**Figuur 2: GO-termen vertegenwoordiging** In het figuur is te zien hoe de vertegenwoordiging van de GO-termen verdeeld is tussen de 20 meest vertegenwoordigde. Op de x-as is de p waarde te zien en op de y-as de GO-termen.

### Moleculaire pathway

De gekozen pathway is [TNF signaling]( Resultaten/hsa04668.png), in **Figuur 3** is de pathway analyse te zien. Deze toont de opregulatie van genen betrokken bij de instroom en activering van leukocyten. Ook is er een opregulatie van inflammatoire cytokines en signalering.
<p align="center">
  <img src=" Resultaten/hsa04668.pathview.png" alt="Pathview " width="400"/>
</p>
**Figuur 3: Pathview analyse** In het figuur is te zien welke genen binnen de pathway op- en neergereguleerd waren. Rood staat voor opregulatie van het gen en groen neerregulatie.


## Conclusie

Spreuken met meer accuraatheid lijken minder krachtig te zijn. Een uitzondering op deze trend is de onvergeeflijke vloek *Avada Kedavra*, welke beter niet gebruikt kan worden. 

