## Electronic Supplementary Information for "The biogeography of population resilience in lowland South America" chapter

This repository holds the code and data for reproducing the analysis in the chapter "The biogeography of population resilience in lowland South America", forthcoming in the _Oxford Handbook of Resilience in Climate History_. It also contains supplementary information on the statistical modelling undertaken as part of the aforementioned work. 

The code comprises five principal sections, plus setup:

0. Setup & data loading
1. Data processing and display
2. Radiocarbon analysis
3. Statistical modelling
4. Output
5. Supplementary output

The code features an extended version of the [p2pPerm function](https://github.com/philriris/p2pPerm) that was first introduced in [Riris and De Souza (2021)](https://doi.org/10.3389/fevo.2021.740629), here called the `resmet` (**res**ilience **met**rics) function. 

In addition, the code is accompanied by three datasets:

- A [table](https://github.com/philriris/resilience-handbook-chapter/blob/main/rhdata.csv) containing archaeological radiocarbon dates from lowland tropical South America 
- A [shapefile](https://github.com/philriris/resilience-handbook-chapter/tree/main/sa_eco) of South American ecoregions, original data available [here](http://ecologicalregions.info/data/sa/)
- A [table](https://github.com/philriris/resilience-handbook-chapter/blob/main/domesticates.csv) of domesticated Neotropical plants resolved in Amazonian palaeoecological records, after [Iriarte et al. (2020)](https://doi.org/10.1016/j.quascirev.2020.106582)

Briefly, the georeferenced radiocarbon data used in the paper have been compiled from a wide range of sources, including Goldberg et al. (2016), Riris & Arroyo-Kalin (2019), Napolitano et al. (2019) Arroyo-Kalin & Riris (2021), De Souza & Riris (2021), and Bird et al. (2022). These sources have been extensively cross-checked for duplicate lab codes, variation in site naming conventions, and reported locations, in order to minimise errors arising from these variables. It does not purport to be error-free, although it is adequate for the current analysis. 

As well as its use in the production of Figure 1, the ecoregions shapefile has been intersected with the radiocarbon date locations to append this information to `rhdata.csv`. An extended description and rationale for its use can be found in the main text. 

For the convenience of the end-user, an additional file containing the main results ([metrics_regular.csv](https://github.com/philriris/resilience-handbook-chapter/blob/main/metrics_regular.csv)) is included. The data contained in this table is the subject of section 3 of the code, **Statistical Modelling**. It forms the basis of the discussion in the chapter. 

Data cleaning was carried out manually on the raw output of the `resmet` function to remove false positives from the table. These "events" are either: _a)_ statistically significant downturns present in periods where, logically, no humans should be present, e.g. in the Greater Antilles before ~6000 cal BP, or: _b)_ downturns where there are no minima in the summed probability distributions of calibrated radiocarbon dates, returning nonsensical resilience metrics. Removing these data rows introduces errors to the variable Cumulative, which counts the cumulative number of downturns detected by the `permTest` function in `rcarbon`. The file version of the output in this repository should be considered authoritative for present purposes, as these counting errors in Cumulative have been manually fixed too. 

### References

- Arroyo-Kalin, M. and Riris, P. 2021. Did pre-Columbian populations of the Amazonian biome reach carrying capacity during the Late Holocene? *Phil. Trans. R. Soc. B* 376: 20190715 http://doi.org/10.1098/rstb.2019.0715

- Bird, D., Miranda, L., Vander Linden, M., Robinson, E., Bocinsky, R.K., Nicholson, C., Capriles, J.M., Finley, J.B., Gayo, E.M., Gil, A. and d’Alpoim Guedes, J. 2022. p3k14c, a synthetic global database of archaeological radiocarbon dates. *Scientific Data* 9: 1-19. https://doi.org/10.1038/s41597-022-01118-7

- De Souza, J.G. & Riris, P. 2021. Delayed demographic transition following the adoption of cultivated plants in the eastern La Plata Basin and Atlantic coast, South America. *Journal of Archaeological Science*. 125: 105293. https://doi.org/10.1016/j.jas.2020.105293

- Goldberg, A., Mychajliw, A.M. and Hadly, E.A. 2016. Post-invasion demography of prehistoric humans in South America. *Nature* 532: 232-235. https://doi.org/10.1038/nature17176 

- Iriarte, J., Elliott, S., Maezumi, S.Y., Alves, D., Gonda, R., Robinson, M., de Souza, J.G., Watling, J. and Handley, J. 2020. The origins of Amazonian landscapes: Plant cultivation, domestication and the spread of food production in tropical South America. _Quaternary Science Reviews_ 248: 106582. https://doi.org/10.1016/j.quascirev.2020.106582 

- Napolitano MF, DiNapoli RJ, Stone JH, Levin MJ, Jew NP, Lane BG, O’Connor JT, Fitzpatrick SM. 2019. Reevaluating human colonization of the Caribbean using chronometric hygiene and Bayesian modeling. _Science Advances_. 5: eaar7806. https://doi.org/10.1126/sciadv.aar7806

- Riris, P. and Arroyo-Kalin, M. 2019. Widespread population decline in South America correlates with mid-Holocene climate change. *Scientific Reports* 9: 6850. https://doi.org/10.1038/s41598-019-43086-w

- Riris, P. and De Souza, J.G. 2021. Formal tests for resistance-resilience in archaeological time series. _Frontiers in Ecology and Evolution_, 9. https://doi.org/10.3389/fevo.2021.740629

