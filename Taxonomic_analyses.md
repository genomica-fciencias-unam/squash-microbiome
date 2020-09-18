
# Bacterial diversity analyses and interaction with plant phenotype


```R
library(phyloseq)
library(ggplot2)
library(ape)
library(vegan)
library(gridExtra)
library(phylogeo)
library(ggrepel)
library(reshape2)
library(RColorBrewer)
library(ggsignif)
library(forcats)
library(tidyverse)
library(limma)
library(microbiome)
```

    Loading required package: permute
    Loading required package: lattice
    This is vegan 2.5-5
    Warning message:
    “replacing previous import ‘dplyr::combine’ by ‘gridExtra::combine’ when loading ‘phylogeo’”── [1mAttaching packages[22m ─────────────────────────────────────── tidyverse 1.2.1 ──
    [32m✔[39m [34mtibble [39m 2.1.3     [32m✔[39m [34mpurrr  [39m 0.3.3
    [32m✔[39m [34mtidyr  [39m 1.0.0     [32m✔[39m [34mdplyr  [39m 0.8.3
    [32m✔[39m [34mreadr  [39m 1.3.1     [32m✔[39m [34mstringr[39m 1.4.0
    ── [1mConflicts[22m ────────────────────────────────────────── tidyverse_conflicts() ──
    [31m✖[39m [34mdplyr[39m::[32mcombine()[39m masks [34mgridExtra[39m::combine()
    [31m✖[39m [34mdplyr[39m::[32mfilter()[39m  masks [34mstats[39m::filter()
    [31m✖[39m [34mdplyr[39m::[32mlag()[39m     masks [34mstats[39m::lag()
    
    microbiome R package (microbiome.github.com)
        
    
    
     Copyright (C) 2011-2018 Leo Lahti et al. <microbiome.github.io>
    
    
    Attaching package: ‘microbiome’
    
    The following object is masked from ‘package:base’:
    
        transform
    



```R
# Load OTU table
otu <- as.matrix(read.table("otusquash.txt", header=T, row.names=1))
colnames(otu) <- c("16S_3_Fe", "16S_3_Fr", "16S_3_Ie", "16S_3_Ir", "16S_3_Ze", "16S_Ce", "16S_2_Ze", "16S_2_Zr", "16S_1_Ze", "16S_4_Ir", "16S_4_Ze", "16S_5_Le", "16S_5_Ze", "16S_5_Zr", "16S_4_Ie", "16S_3_Le", "16S_3_Lr", "16S_3_Ls", "16S_5_Ls", "16S_3_Zr", "16S_4_Zr", "16S_Cr", "16S_2_Le", "16S_2_Ls", "16S_1_Ir", "16S_1_Ls", "16S_4_Fr", "16S_4_Le", "16S_4_Ls", "16S_5_Ir", "16S_5_Lr", "16S_1_Fr", "16S_1_Fe", "16S_1_Zr", "16S_2_Ie", "16S_4_Fe", "16S_2_Fe", "16S_2_Fr", "16S_2_Lr", "16S_1_Ie", "16S_1_Le", "16S_1_Lr", "16S_2_Ir", "16S_5_Ie", "16S_4_Lr", "16S_5_Fr", "16S_5_Fe", "16S_Cn")
OTU <- otu_table(otu, taxa_are_rows=T)
```


```R
# Load taxonomy table
taxa <- as.matrix(read.table("squash.tax", header=T, row.names=1))
TAXA <- tax_table(taxa)
```


```R
# Create phyloseq object with taxonomy and OTU tables 
calabacita <-phyloseq(OTU,TAXA)
```


```R
# Load metadata and denote variables
metadata <- read.csv("sqmetadata.csv", header=T)
sampledata <- sample_data(data.frame(sample=metadata$sample, name=metadata$name, locality=metadata$locality, latitude=metadata$Latitude, longitude=metadata$Longitude, type=metadata$type, treatment=metadata$treatment, surface=metadata$surface, alenght=metadata$alength, stlenght=metadata$stlenght, flowers=metadata$flowers, leaves=metadata$leaves, rlenght=metadata$rlenght, sdiameter=metadata$sdiameter, chlorophyll=metadata$chlorophyll, tchlorophyll=metadata$tchlorophyll, carotenoids=metadata$carotenoids, tcarotenoids=metadata$tcarotenoids, abiomass=metadata$abiomass, sla=metadata$sla, toc=metadata$TOC, tn=metadata$TN, nh4=metadata$NH4, no3=metadata$NO3, pt=metadata$PT, hpo4=metadata$HPO4, cn=metadata$C.N, cp=metadata$C.N, np=metadata$N.P, ph=metadata$pH, mat=metadata$MAT, map=metadata$MAP, ai=metadata$AI, climate=metadata$climate, soil=metadata$soil, climetrat=metadata$climetrat, samplec=metadata$samplec, samplet=metadata$samplet, local_commongarden=metadata$local_commongarden, specie=metadata$specie,row.names=sample_names(calabacita)))
sampledata
```


<table>
<thead><tr><th></th><th scope=col>sample</th><th scope=col>name</th><th scope=col>locality</th><th scope=col>latitude</th><th scope=col>longitude</th><th scope=col>type</th><th scope=col>treatment</th><th scope=col>surface</th><th scope=col>alenght</th><th scope=col>stlenght</th><th scope=col>⋯</th><th scope=col>mat</th><th scope=col>map</th><th scope=col>ai</th><th scope=col>climate</th><th scope=col>soil</th><th scope=col>climetrat</th><th scope=col>samplec</th><th scope=col>samplet</th><th scope=col>local_commongarden</th><th scope=col>specie</th></tr></thead>
<tbody>
	<tr><th scope=row>16S_3_Fe</th><td>16S_3_Fe                    </td><td>3_Fe                        </td><td>3                           </td><td>     NA                     </td><td>       NA                   </td><td>endosphere                  </td><td>fertilized                  </td><td>312.864                     </td><td>14.49                       </td><td>5.5                         </td><td>⋯                           </td><td>18.40                       </td><td> 705                        </td><td>38.30                       </td><td>arid                        </td><td>vertisol                    </td><td>arid                        </td><td>endosphere_arid             </td><td>fertilized_endosphere       </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_3_Fr</th><td>16S_3_Fr                    </td><td>3_Fr                        </td><td>3                           </td><td>     NA                     </td><td>       NA                   </td><td>rhizosphere                 </td><td>fertilized                  </td><td>312.864                     </td><td>14.49                       </td><td>5.5                         </td><td>⋯                           </td><td>18.40                       </td><td> 705                        </td><td>38.30                       </td><td>arid                        </td><td>vertisol                    </td><td>arid                        </td><td>rhizosphere_arid            </td><td>fertilized_rhizosphere      </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_3_Ie</th><td>16S_3_Ie                    </td><td>3_Ie                        </td><td>3                           </td><td>     NA                     </td><td>       NA                   </td><td>endosphere                  </td><td>inoculum                    </td><td>213.190                     </td><td>14.61                       </td><td>7.5                         </td><td>⋯                           </td><td>18.40                       </td><td> 705                        </td><td>38.30                       </td><td>arid                        </td><td>vertisol                    </td><td>arid                        </td><td>endosphere_arid             </td><td>inoculum_endosphere         </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_3_Ir</th><td>16S_3_Ir                    </td><td>3_Ir                        </td><td>3                           </td><td>     NA                     </td><td>       NA                   </td><td>rhizosphere                 </td><td>inoculum                    </td><td>213.190                     </td><td>14.61                       </td><td>7.5                         </td><td>⋯                           </td><td>18.40                       </td><td> 705                        </td><td>38.30                       </td><td>arid                        </td><td>vertisol                    </td><td>arid                        </td><td>rhizosphere_arid            </td><td>inoculum_rhizosphere        </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_3_Ze</th><td>16S_3_Ze                    </td><td>3_Ze                        </td><td>3                           </td><td>     NA                     </td><td>       NA                   </td><td>endosphere                  </td><td>sterilized                  </td><td>465.203                     </td><td>16.04                       </td><td>7.0                         </td><td>⋯                           </td><td>18.40                       </td><td> 705                        </td><td>38.30                       </td><td>arid                        </td><td>vertisol                    </td><td>arid_sterilized             </td><td>endosphere_arid_sterilized  </td><td>sterilized_endosphere       </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_Ce</th><td>16S_Ce                      </td><td>Ce                          </td><td>c                           </td><td>     NA                     </td><td>       NA                   </td><td>endosphere                  </td><td>sterilized                  </td><td> 17.469                     </td><td> 4.70                       </td><td>2.8                         </td><td>⋯                           </td><td>18.40                       </td><td>  NA                        </td><td>   NA                       </td><td>NA                          </td><td>substrate                   </td><td>none                        </td><td>endosphere_none             </td><td>sterilized_endosphere       </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_2_Ze</th><td>16S_2_Ze                    </td><td>2_Ze                        </td><td>2                           </td><td>     NA                     </td><td>       NA                   </td><td>endosphere                  </td><td>sterilized                  </td><td>100.904                     </td><td>12.75                       </td><td>3.1                         </td><td>⋯                           </td><td>15.75                       </td><td> 801                        </td><td>50.80                       </td><td>humid                       </td><td>phaeozem                    </td><td>humid_sterilized            </td><td>endosphere_humid_sterilized </td><td>sterilized_endosphere       </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_2_Zr</th><td>16S_2_Zr                    </td><td>2_Zr                        </td><td>2                           </td><td>     NA                     </td><td>       NA                   </td><td>rhizosphere                 </td><td>sterilized                  </td><td>100.904                     </td><td>12.75                       </td><td>3.1                         </td><td>⋯                           </td><td>15.75                       </td><td> 801                        </td><td>50.80                       </td><td>humid                       </td><td>phaeozem                    </td><td>humid_sterilized            </td><td>rhizosphere_humid_sterilized</td><td>sterilized_rhizosphere      </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_1_Ze</th><td>16S_1_Ze                    </td><td>1_Ze                        </td><td>1                           </td><td>     NA                     </td><td>       NA                   </td><td>endosphere                  </td><td>sterilized                  </td><td>140.477                     </td><td>12.59                       </td><td>4.3                         </td><td>⋯                           </td><td>19.70                       </td><td>1958                        </td><td>99.40                       </td><td>humid                       </td><td>vertisol                    </td><td>humid_sterilized            </td><td>endosphere_humid_sterilized </td><td>sterilized_endosphere       </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_4_Ir</th><td>16S_4_Ir                    </td><td>4_Ir                        </td><td>4                           </td><td>     NA                     </td><td>       NA                   </td><td>rhizosphere                 </td><td>inoculum                    </td><td>141.820                     </td><td>18.03                       </td><td>8.5                         </td><td>⋯                           </td><td>17.27                       </td><td> 427                        </td><td>24.70                       </td><td>arid                        </td><td>phaeozem                    </td><td>arid                        </td><td>rhizosphere_arid            </td><td>inoculum_rhizosphere        </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_4_Ze</th><td>16S_4_Ze                    </td><td>4_Ze                        </td><td>4                           </td><td>     NA                     </td><td>       NA                   </td><td>endosphere                  </td><td>sterilized                  </td><td>118.356                     </td><td>16.51                       </td><td>8.0                         </td><td>⋯                           </td><td>17.27                       </td><td> 427                        </td><td>24.70                       </td><td>arid                        </td><td>phaeozem                    </td><td>arid_sterilized             </td><td>endosphere_arid_sterilized  </td><td>sterilized_endosphere       </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_5_Le</th><td>16S_5_Le                    </td><td>5_Le                        </td><td>5                           </td><td>22.4528                     </td><td>-100.6976                   </td><td>endosphere                  </td><td>local                       </td><td>     NA                     </td><td>   NA                       </td><td> NA                         </td><td>⋯                           </td><td>18.33                       </td><td> 386                        </td><td>21.06                       </td><td>arid                        </td><td>xerosol                     </td><td>arid                        </td><td>endosphere_arid             </td><td>local_endosphere            </td><td>L                           </td><td>Cucurbita_pepo              </td></tr>
	<tr><th scope=row>16S_5_Ze</th><td>16S_5_Ze                    </td><td>5_Ze                        </td><td>5                           </td><td>     NA                     </td><td>       NA                   </td><td>endosphere                  </td><td>sterilized                  </td><td>214.452                     </td><td>17.58                       </td><td>8.8                         </td><td>⋯                           </td><td>18.33                       </td><td> 386                        </td><td>21.06                       </td><td>arid                        </td><td>xerosol                     </td><td>arid_sterilized             </td><td>endosphere_arid_sterilized  </td><td>sterilized_endosphere       </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_5_Zr</th><td>16S_5_Zr                    </td><td>5_Zr                        </td><td>5                           </td><td>     NA                     </td><td>       NA                   </td><td>rhizosphere                 </td><td>sterilized                  </td><td>214.452                     </td><td>17.58                       </td><td>8.8                         </td><td>⋯                           </td><td>18.33                       </td><td> 386                        </td><td>21.06                       </td><td>arid                        </td><td>xerosol                     </td><td>arid_sterilized             </td><td>rhizosphere_arid_sterilized </td><td>sterilized_rhizosphere      </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_4_Ie</th><td>16S_4_Ie                    </td><td>4_Ie                        </td><td>4                           </td><td>     NA                     </td><td>       NA                   </td><td>endosphere                  </td><td>inoculum                    </td><td>141.820                     </td><td>18.03                       </td><td>8.5                         </td><td>⋯                           </td><td>17.27                       </td><td> 427                        </td><td>24.70                       </td><td>arid                        </td><td>phaeozem                    </td><td>arid                        </td><td>endosphere_arid             </td><td>inoculum_endosphere         </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_3_Le</th><td>16S_3_Le                    </td><td>3_Le                        </td><td>3                           </td><td>20.0612                     </td><td>-100.8154                   </td><td>endosphere                  </td><td>local                       </td><td>     NA                     </td><td>   NA                       </td><td> NA                         </td><td>⋯                           </td><td>18.40                       </td><td> 705                        </td><td>38.30                       </td><td>arid                        </td><td>vertisol                    </td><td>arid                        </td><td>endosphere_arid             </td><td>local_endosphere            </td><td>L                           </td><td>Cucurbita_sp                </td></tr>
	<tr><th scope=row>16S_3_Lr</th><td>16S_3_Lr                    </td><td>3_Lr                        </td><td>3                           </td><td>20.0612                     </td><td>-100.8154                   </td><td>rhizosphere                 </td><td>local                       </td><td>     NA                     </td><td>   NA                       </td><td> NA                         </td><td>⋯                           </td><td>18.40                       </td><td> 705                        </td><td>38.30                       </td><td>arid                        </td><td>vertisol                    </td><td>arid                        </td><td>rhizosphere_arid            </td><td>local_rhizosphere           </td><td>L                           </td><td>Cucurbita_sp                </td></tr>
	<tr><th scope=row>16S_3_Ls</th><td>16S_3_Ls                    </td><td>3_Ls                        </td><td>3                           </td><td>20.0612                     </td><td>-100.8154                   </td><td>soil                        </td><td>local                       </td><td>     NA                     </td><td>   NA                       </td><td> NA                         </td><td>⋯                           </td><td>18.40                       </td><td> 705                        </td><td>38.30                       </td><td>arid                        </td><td>vertisol                    </td><td>arid                        </td><td>soil_arid                   </td><td>local_soil                  </td><td>L                           </td><td>No_plant                    </td></tr>
	<tr><th scope=row>16S_5_Ls</th><td>16S_5_Ls                    </td><td>5_Ls                        </td><td>5                           </td><td>22.4528                     </td><td>-100.6976                   </td><td>soil                        </td><td>local                       </td><td>     NA                     </td><td>   NA                       </td><td> NA                         </td><td>⋯                           </td><td>18.33                       </td><td> 386                        </td><td>21.06                       </td><td>arid                        </td><td>xerosol                     </td><td>arid                        </td><td>soil_arid                   </td><td>local_soil                  </td><td>L                           </td><td>No_plant                    </td></tr>
	<tr><th scope=row>16S_3_Zr</th><td>16S_3_Zr                    </td><td>3_Zr                        </td><td>3                           </td><td>     NA                     </td><td>       NA                   </td><td>rhizosphere                 </td><td>sterilized                  </td><td>465.203                     </td><td>16.04                       </td><td>7.0                         </td><td>⋯                           </td><td>18.40                       </td><td> 705                        </td><td>38.30                       </td><td>arid                        </td><td>vertisol                    </td><td>arid_sterilized             </td><td>rhizosphere_arid_sterilized </td><td>sterilized_rhizosphere      </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_4_Zr</th><td>16S_4_Zr                    </td><td>4_Zr                        </td><td>4                           </td><td>     NA                     </td><td>       NA                   </td><td>rhizosphere                 </td><td>sterilized                  </td><td>118.356                     </td><td>16.51                       </td><td>8.0                         </td><td>⋯                           </td><td>17.27                       </td><td> 427                        </td><td>24.70                       </td><td>arid                        </td><td>phaeozem                    </td><td>arid_sterilized             </td><td>rhizosphere_arid_sterilized </td><td>sterilized_rhizosphere      </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_Cr</th><td>16S_Cr                      </td><td>Cr                          </td><td>c                           </td><td>     NA                     </td><td>       NA                   </td><td>rhizosphere                 </td><td>sterilized                  </td><td> 17.469                     </td><td> 4.70                       </td><td>2.8                         </td><td>⋯                           </td><td>   NA                       </td><td>  NA                        </td><td>   NA                       </td><td>NA                          </td><td>substrate                   </td><td>none                        </td><td>rhizosphere_none            </td><td>sterilized_rhizosphere      </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_2_Le</th><td>16S_2_Le                    </td><td>2_Le                        </td><td>2                           </td><td>20.3874                     </td><td>-100.2825                   </td><td>endosphere                  </td><td>local                       </td><td>     NA                     </td><td>   NA                       </td><td> NA                         </td><td>⋯                           </td><td>15.75                       </td><td> 801                        </td><td>50.80                       </td><td>humid                       </td><td>phaeozem                    </td><td>humid                       </td><td>endosphere_humid            </td><td>local_endosphere            </td><td>L                           </td><td>Cucurbita_pepo              </td></tr>
	<tr><th scope=row>16S_2_Ls</th><td>16S_2_Ls                    </td><td>2_Ls                        </td><td>2                           </td><td>20.3874                     </td><td>-100.2825                   </td><td>soil                        </td><td>local                       </td><td>     NA                     </td><td>   NA                       </td><td> NA                         </td><td>⋯                           </td><td>15.75                       </td><td> 801                        </td><td>50.80                       </td><td>humid                       </td><td>phaeozem                    </td><td>humid                       </td><td>soil_humid                  </td><td>local_soil                  </td><td>L                           </td><td>No_plant                    </td></tr>
	<tr><th scope=row>16S_1_Ir</th><td>16S_1_Ir                    </td><td>1_Ir                        </td><td>1                           </td><td>     NA                     </td><td>       NA                   </td><td>rhizosphere                 </td><td>inoculum                    </td><td> 80.117                     </td><td> 8.67                       </td><td>4.2                         </td><td>⋯                           </td><td>19.70                       </td><td>1958                        </td><td>99.40                       </td><td>humid                       </td><td>vertisol                    </td><td>humid                       </td><td>rhizosphere_humid           </td><td>inoculum_rhizosphere        </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_1_Ls</th><td>16S_1_Ls                    </td><td>1_Ls                        </td><td>1                           </td><td>18.8060                     </td><td> -97.0843                   </td><td>soil                        </td><td>local                       </td><td>     NA                     </td><td>   NA                       </td><td> NA                         </td><td>⋯                           </td><td>19.70                       </td><td>1958                        </td><td>99.40                       </td><td>humid                       </td><td>vertisol                    </td><td>humid                       </td><td>soil_humid                  </td><td>local_soil                  </td><td>L                           </td><td>No_plant                    </td></tr>
	<tr><th scope=row>16S_4_Fr</th><td>16S_4_Fr                    </td><td>4_Fr                        </td><td>4                           </td><td>     NA                     </td><td>       NA                   </td><td>rhizosphere                 </td><td>fertilized                  </td><td>220.129                     </td><td>17.98                       </td><td>7.5                         </td><td>⋯                           </td><td>17.27                       </td><td> 427                        </td><td>24.70                       </td><td>arid                        </td><td>phaeozem                    </td><td>arid                        </td><td>rhizosphere_arid            </td><td>fertilized_rhizosphere      </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_4_Le</th><td>16S_4_Le                    </td><td>4_Le                        </td><td>4                           </td><td>21.5945                     </td><td>-101.0839                   </td><td>endosphere                  </td><td>local                       </td><td>     NA                     </td><td>   NA                       </td><td> NA                         </td><td>⋯                           </td><td>17.27                       </td><td> 427                        </td><td>24.70                       </td><td>arid                        </td><td>phaeozem                    </td><td>arid                        </td><td>endosphere_arid             </td><td>local_endosphere            </td><td>L                           </td><td>Cucurbita_pepo              </td></tr>
	<tr><th scope=row>16S_4_Ls</th><td>16S_4_Ls                    </td><td>4_Ls                        </td><td>4                           </td><td>21.5945                     </td><td>-101.0839                   </td><td>soil                        </td><td>local                       </td><td>     NA                     </td><td>   NA                       </td><td> NA                         </td><td>⋯                           </td><td>17.27                       </td><td> 427                        </td><td>24.70                       </td><td>arid                        </td><td>phaeozem                    </td><td>arid                        </td><td>soil_arid                   </td><td>local_soil                  </td><td>L                           </td><td>No_plant                    </td></tr>
	<tr><th scope=row>16S_5_Ir</th><td>16S_5_Ir                    </td><td>5_Ir                        </td><td>5                           </td><td>     NA                     </td><td>       NA                   </td><td>rhizosphere                 </td><td>inoculum                    </td><td>134.080                     </td><td>14.55                       </td><td>5.7                         </td><td>⋯                           </td><td>18.33                       </td><td> 386                        </td><td>21.06                       </td><td>arid                        </td><td>xerosol                     </td><td>arid                        </td><td>rhizosphere_arid            </td><td>inoculum_rhizosphere        </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_5_Lr</th><td>16S_5_Lr                    </td><td>5_Lr                        </td><td>5                           </td><td>22.4528                     </td><td>-100.6976                   </td><td>rhizosphere                 </td><td>local                       </td><td>     NA                     </td><td>   NA                       </td><td> NA                         </td><td>⋯                           </td><td>18.33                       </td><td> 386                        </td><td>21.06                       </td><td>arid                        </td><td>xerosol                     </td><td>arid                        </td><td>rhizosphere_arid            </td><td>local_rhizosphere           </td><td>L                           </td><td>Cucurbita_pepo              </td></tr>
	<tr><th scope=row>16S_1_Fr</th><td>16S_1_Fr                    </td><td>1_Fr                        </td><td>1                           </td><td>     NA                     </td><td>       NA                   </td><td>rhizosphere                 </td><td>fertilized                  </td><td>146.714                     </td><td>14.57                       </td><td>6.5                         </td><td>⋯                           </td><td>19.70                       </td><td>1958                        </td><td>99.40                       </td><td>humid                       </td><td>vertisol                    </td><td>humid                       </td><td>rhizosphere_humid           </td><td>fertilized_rhizosphere      </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_1_Fe</th><td>16S_1_Fe                    </td><td>1_Fe                        </td><td>1                           </td><td>     NA                     </td><td>       NA                   </td><td>endosphere                  </td><td>fertilized                  </td><td>146.714                     </td><td>14.57                       </td><td>6.5                         </td><td>⋯                           </td><td>19.70                       </td><td>1958                        </td><td>99.40                       </td><td>humid                       </td><td>vertisol                    </td><td>humid                       </td><td>endosphere_humid            </td><td>fertilized_endosphere       </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_1_Zr</th><td>16S_1_Zr                    </td><td>1_Zr                        </td><td>1                           </td><td>     NA                     </td><td>       NA                   </td><td>rhizosphere                 </td><td>sterilized                  </td><td>140.477                     </td><td>12.59                       </td><td>4.3                         </td><td>⋯                           </td><td>19.70                       </td><td>1958                        </td><td>99.40                       </td><td>humid                       </td><td>vertisol                    </td><td>humid_sterilized            </td><td>rhizosphere_humid_sterilized</td><td>sterilized_rhizosphere      </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_2_Ie</th><td>16S_2_Ie                    </td><td>2_Ie                        </td><td>2                           </td><td>     NA                     </td><td>       NA                   </td><td>endosphere                  </td><td>inoculum                    </td><td>109.261                     </td><td>15.15                       </td><td>8.0                         </td><td>⋯                           </td><td>15.75                       </td><td> 801                        </td><td>50.80                       </td><td>humid                       </td><td>phaeozem                    </td><td>humid                       </td><td>endosphere_humid            </td><td>inoculum_endosphere         </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_4_Fe</th><td>16S_4_Fe                    </td><td>4_Fe                        </td><td>4                           </td><td>     NA                     </td><td>       NA                   </td><td>endosphere                  </td><td>fertilized                  </td><td>220.129                     </td><td>17.98                       </td><td>7.5                         </td><td>⋯                           </td><td>17.27                       </td><td> 427                        </td><td>24.70                       </td><td>arid                        </td><td>phaeozem                    </td><td>arid                        </td><td>endosphere_arid             </td><td>fertilized_endosphere       </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_2_Fe</th><td>16S_2_Fe                    </td><td>2_Fe                        </td><td>2                           </td><td>     NA                     </td><td>       NA                   </td><td>endosphere                  </td><td>fertilized                  </td><td>216.119                     </td><td>15.86                       </td><td>4.5                         </td><td>⋯                           </td><td>15.75                       </td><td> 801                        </td><td>50.80                       </td><td>humid                       </td><td>phaeozem                    </td><td>humid                       </td><td>endosphere_humid            </td><td>fertilized_endosphere       </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_2_Fr</th><td>16S_2_Fr                    </td><td>2_Fr                        </td><td>2                           </td><td>     NA                     </td><td>       NA                   </td><td>rhizosphere                 </td><td>fertilized                  </td><td>216.119                     </td><td>15.86                       </td><td>4.5                         </td><td>⋯                           </td><td>15.75                       </td><td> 801                        </td><td>50.80                       </td><td>humid                       </td><td>phaeozem                    </td><td>humid                       </td><td>rhizosphere_humid           </td><td>fertilized_rhizosphere      </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_2_Lr</th><td>16S_2_Lr                    </td><td>2_Lr                        </td><td>2                           </td><td>20.3874                     </td><td>-100.2825                   </td><td>rhizosphere                 </td><td>local                       </td><td>     NA                     </td><td>   NA                       </td><td> NA                         </td><td>⋯                           </td><td>15.75                       </td><td> 801                        </td><td>50.80                       </td><td>humid                       </td><td>phaeozem                    </td><td>humid                       </td><td>rhizosphere_humid           </td><td>local_rhizosphere           </td><td>L                           </td><td>Cucurbita_pepo              </td></tr>
	<tr><th scope=row>16S_1_Ie</th><td>16S_1_Ie                    </td><td>1_Ie                        </td><td>1                           </td><td>     NA                     </td><td>       NA                   </td><td>endosphere                  </td><td>inoculum                    </td><td> 80.117                     </td><td> 8.67                       </td><td>4.2                         </td><td>⋯                           </td><td>19.70                       </td><td>1958                        </td><td>99.40                       </td><td>humid                       </td><td>vertisol                    </td><td>humid                       </td><td>endosphere_humid            </td><td>inoculum_endosphere         </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_1_Le</th><td>16S_1_Le                    </td><td>1_Le                        </td><td>1                           </td><td>18.8060                     </td><td> -97.0843                   </td><td>endosphere                  </td><td>local                       </td><td>     NA                     </td><td>   NA                       </td><td> NA                         </td><td>⋯                           </td><td>19.70                       </td><td>1958                        </td><td>99.40                       </td><td>humid                       </td><td>vertisol                    </td><td>humid                       </td><td>endosphere_humid            </td><td>local_endosphere            </td><td>L                           </td><td>Cucurbita_pepo              </td></tr>
	<tr><th scope=row>16S_1_Lr</th><td>16S_1_Lr                    </td><td>1_Lr                        </td><td>1                           </td><td>18.8060                     </td><td> -97.0843                   </td><td>rhizosphere                 </td><td>local                       </td><td>     NA                     </td><td>   NA                       </td><td> NA                         </td><td>⋯                           </td><td>19.70                       </td><td>1958                        </td><td>99.40                       </td><td>humid                       </td><td>vertisol                    </td><td>humid                       </td><td>rhizosphere_humid           </td><td>local_rhizosphere           </td><td>L                           </td><td>Cucurbita_pepo              </td></tr>
	<tr><th scope=row>16S_2_Ir</th><td>16S_2_Ir                    </td><td>2_Ir                        </td><td>2                           </td><td>     NA                     </td><td>       NA                   </td><td>rhizosphere                 </td><td>inoculum                    </td><td>109.261                     </td><td>15.15                       </td><td>8.0                         </td><td>⋯                           </td><td>15.75                       </td><td> 801                        </td><td>50.80                       </td><td>humid                       </td><td>phaeozem                    </td><td>humid                       </td><td>rhizosphere_humid           </td><td>inoculum_rhizosphere        </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_5_Ie</th><td>16S_5_Ie                    </td><td>5_Ie                        </td><td>5                           </td><td>     NA                     </td><td>       NA                   </td><td>endosphere                  </td><td>inoculum                    </td><td>134.080                     </td><td>14.55                       </td><td>5.7                         </td><td>⋯                           </td><td>18.33                       </td><td> 386                        </td><td>21.06                       </td><td>arid                        </td><td>xerosol                     </td><td>arid                        </td><td>endosphere_arid             </td><td>inoculum_endosphere         </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_4_Lr</th><td>16S_4_Lr                    </td><td>4_Lr                        </td><td>4                           </td><td>21.5945                     </td><td>-101.0839                   </td><td>rhizosphere                 </td><td>local                       </td><td>     NA                     </td><td>   NA                       </td><td> NA                         </td><td>⋯                           </td><td>17.27                       </td><td> 427                        </td><td>24.70                       </td><td>arid                        </td><td>phaeozem                    </td><td>arid                        </td><td>rhizosphere_arid            </td><td>local_rhizosphere           </td><td>L                           </td><td>Cucurbita_pepo              </td></tr>
	<tr><th scope=row>16S_5_Fr</th><td>16S_5_Fe                    </td><td>5_Fe                        </td><td>5                           </td><td>     NA                     </td><td>       NA                   </td><td>endosphere                  </td><td>fertilized                  </td><td>276.373                     </td><td>15.68                       </td><td>6.5                         </td><td>⋯                           </td><td>18.33                       </td><td> 386                        </td><td>21.06                       </td><td>arid                        </td><td>xerosol                     </td><td>arid                        </td><td>endosphere_arid             </td><td>fertilized_endosphere       </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_5_Fe</th><td>16S_5_Fr                    </td><td>5_Fr                        </td><td>5                           </td><td>     NA                     </td><td>       NA                   </td><td>rhizosphere                 </td><td>fertilized                  </td><td>276.373                     </td><td>15.68                       </td><td>6.5                         </td><td>⋯                           </td><td>18.33                       </td><td> 386                        </td><td>21.06                       </td><td>arid                        </td><td>xerosol                     </td><td>arid                        </td><td>rhizosphere_arid            </td><td>fertilized_rhizosphere      </td><td>C                           </td><td>Cucurbita_pepo_var_Zuchini  </td></tr>
	<tr><th scope=row>16S_Cn</th><td>16S_Cn                      </td><td>Cn                          </td><td>cc                          </td><td>     NA                     </td><td>       NA                   </td><td>substrate                   </td><td>substrate                   </td><td>     NA                     </td><td>   NA                       </td><td> NA                         </td><td>⋯                           </td><td>   NA                       </td><td>  NA                        </td><td>   NA                       </td><td>NA                          </td><td>substrate                   </td><td>none                        </td><td>substrate                   </td><td>sterilized_substrate        </td><td>C                           </td><td>No_plant                    </td></tr>
</tbody>
</table>




```R
# Load tree
squash_tree <- read.tree("squash.tree")
```


```R
# Create phyloseq object with metadata, taxonomy, OTU table and tree
calabacita <- phyloseq(OTU, TAXA, sampledata, squash_tree)
```

##ANALYSES OF ALPHA AND BETA DIVERSITY


```R
# Sampling locality map
mapl <- map_phyloseq(calabacita, region="mexico", color="climate") + 
        scale_colour_manual(values = c("#ff9a58ff", "#5c5cafff", "black", "black")) + 
        theme_bw()
        
ggsave("mapl.pdf", width = 25, height = 15.85, units = "cm")
"Sampling map by aridity"
mapl
```

    Warning message:
    “Removed 33 rows containing missing values (geom_point).”


'Sampling map by aridity'


    Warning message:
    “Removed 33 rows containing missing values (geom_point).”


![png](output_9_3.png)



```R
# Aridity index map
mapai <- map_phyloseq(calabacita, region="mexico", color="ai") +
        theme_bw()

print("Sampling map by aridity index")
mapai 
```

    [1] "Sampling map by aridity index"


    Warning message:
    “Removed 33 rows containing missing values (geom_point).”


![png](output_10_2.png)



```R
# Shannon diversity estimations
#Get diversity measurments
richness <- estimate_richness(calabacita)
write.csv(richness, "calabacita_richness.csv")

head(richness)
```


<table>
<thead><tr><th></th><th scope=col>Observed</th><th scope=col>Chao1</th><th scope=col>se.chao1</th><th scope=col>ACE</th><th scope=col>se.ACE</th><th scope=col>Shannon</th><th scope=col>Simpson</th><th scope=col>InvSimpson</th><th scope=col>Fisher</th></tr></thead>
<tbody>
	<tr><th scope=row>X16S_3_Fe</th><td> 4695    </td><td> 8805.709</td><td>232.1314 </td><td> 9492.881</td><td> 62.13737</td><td>4.669608 </td><td>0.9236138</td><td>13.09138 </td><td>1089.100 </td></tr>
	<tr><th scope=row>X16S_3_Fr</th><td> 6393    </td><td>13004.281</td><td>303.4851 </td><td>14771.529</td><td> 81.31639</td><td>5.126511 </td><td>0.9480213</td><td>19.23866 </td><td>1628.907 </td></tr>
	<tr><th scope=row>X16S_3_Ie</th><td>12976    </td><td>24019.184</td><td>361.2264 </td><td>27464.543</td><td>109.03461</td><td>6.614258 </td><td>0.9854015</td><td>68.49996 </td><td>3552.983 </td></tr>
	<tr><th scope=row>X16S_3_Ir</th><td>13466    </td><td>24584.621</td><td>352.6581 </td><td>28784.424</td><td>115.21287</td><td>6.733750 </td><td>0.9857913</td><td>70.37946 </td><td>3937.753 </td></tr>
	<tr><th scope=row>X16S_3_Ze</th><td> 7963    </td><td>14550.805</td><td>271.5872 </td><td>16573.097</td><td> 85.54807</td><td>5.098612 </td><td>0.9438262</td><td>17.80188 </td><td>2052.203 </td></tr>
	<tr><th scope=row>X16S_Ce</th><td> 5219    </td><td> 9099.639</td><td>202.7416 </td><td>10217.439</td><td> 65.31094</td><td>4.699306 </td><td>0.9607047</td><td>25.44832 </td><td>1111.184 </td></tr>
</tbody>
</table>




```R
#Plot shannon diversity
div_tab <- read.csv("calabacita_richness.csv")
div_tab_data <- cbind(div_tab, sampledata)
div_tab_data$name <- factor(div_tab_data$name, levels=c("1_Ls", "1_Lr", "1_Le", "1_Ir", "1_Ie", "1_Fr", "1_Fe", "1_Zr", "1_Ze", "2_Ls", "2_Lr", "2_Le", "2_Ir", "2_Ie", "2_Fr", "2_Fe", "2_Zr", "2_Ze", "3_Ls", "3_Lr", "3_Le", "3_Ir", "3_Ie", "3_Fr", "3_Fe", "3_Zr", "3_Ze", "4_Ls", "4_Lr", "4_Le", "4_Ir", "4_Ie", "4_Fr", "4_Fe", "4_Zr", "4_Ze", "5_Ls", "5_Lr", "5_Le", "5_Ir", "5_Ie", "5_Fr", "5_Fe", "5_Zr", "5_Ze", "Cr", "Ce", "Cn"))
shannon <- ggplot(div_tab_data, aes(name, Shannon, color=climate)) + geom_point(stat="identity", size=2.5) +
                  scale_colour_manual(values = c("#ff9a58ff", "#5c5cafff"), na.value = "black") + 
                  theme_bw()  + theme(axis.text.x = element_text(angle = 90))
ggsave("shannon.pdf", width = 35, height = 7.5, units = "cm")

print("Shannon diversity by sample")
shannon
```

    [1] "Shannon diversity by sample"



![png](output_12_1.png)



```R
#Estimate relative abundance
rel_phy <- transform_sample_counts(calabacita, function(x) x / sum(x))
                                   
#Most abundant phylum
#Tax glom at Phylum level
phy_calabacita <- tax_glom(rel_phy, "Phylum")
phy_calabacita
```


    phyloseq-class experiment-level object
    otu_table()   OTU Table:         [ 51 taxa and 48 samples ]
    sample_data() Sample Data:       [ 48 samples by 40 sample variables ]
    tax_table()   Taxonomy Table:    [ 51 taxa by 7 taxonomic ranks ]
    phy_tree()    Phylogenetic Tree: [ 51 tips and 50 internal nodes ]



```R
#Export the file with taxa names and taxa counts
calabacita_otu_phy <- otu_table(phy_calabacita)
calabacita_tax_phy <- tax_table(phy_calabacita)
calabacita_otu_tax <- cbind(calabacita_otu_phy, calabacita_tax_phy)
calabacita_otu_tax <- calabacita_otu_tax[, 1:50]
write.table(calabacita_otu_tax, "calabacita_otu_tax.tsv", sep = "\t")

calabacita_otu_tax <- read.table("calabacita_otu_tax.tsv", header=TRUE, row.names=1, 
                                 stringsAsFactors = FALSE, check.names=FALSE)

m_calabacita_otu_tax <- melt(calabacita_otu_tax, measure.vars = c(
    "16S_1_Ls", "16S_1_Lr", "16S_1_Le", "16S_1_Ir", "16S_1_Ie", "16S_1_Fr", "16S_1_Fe", "16S_1_Zr", "16S_1_Ze", 
    "16S_2_Ls", "16S_2_Lr", "16S_2_Le", "16S_2_Ir", "16S_2_Ie", "16S_2_Fr", "16S_2_Fe", "16S_2_Zr", "16S_2_Ze", 
    "16S_3_Ls", "16S_3_Lr", "16S_3_Le", "16S_3_Ir", "16S_3_Ie", "16S_3_Fr", "16S_3_Fe", "16S_3_Zr", "16S_3_Ze", 
    "16S_4_Ls", "16S_4_Lr", "16S_4_Le", "16S_4_Ir", "16S_4_Ie", "16S_4_Fr", "16S_4_Fe", "16S_4_Zr", "16S_4_Ze", 
    "16S_5_Ls", "16S_5_Lr", "16S_5_Le", "16S_5_Ir", "16S_5_Ie", "16S_5_Fr", "16S_5_Fe", "16S_5_Zr", "16S_5_Ze", 
    "16S_Cr", "16S_Ce", "16S_Cn"))
colnames(m_calabacita_otu_tax) <- c("Kingdom", "Phylum", "Sample", "Relative_abundance")
m_calabacita_otu_tax <- m_calabacita_otu_tax[, 2:4]

#Collapse Phylum with low abundance
m_calabacita_otu_tax$Phylum[m_calabacita_otu_tax$Relative_abundance <= 0.015] <- "Low_abundance"
cm_calabacita_otu_tax <- aggregate(m_calabacita_otu_tax$Relative_abundance 
                                          ,by=list(m_calabacita_otu_tax$Phylum, 
                                           m_calabacita_otu_tax$Sample),sum)
colnames(cm_calabacita_otu_tax) <- c("Phylum", "Sample", "Relative_abundance")

#Plot most abundant Phyla
phylum_plot <- ggplot(cm_calabacita_otu_tax, aes(x=Sample, y=Relative_abundance,
                                          fill=fct_reorder(Phylum, Relative_abundance))) + 
        scale_fill_manual(values = c("p__Acidobacteria" = "#dd87b4ff", 
                              "p__Actinobacteria" = "#d26d7aff",
                              "p__Bacteroidetes" = "#b6742aff", 
                              "p__Chloroflexi" = "#f6ef32ff",                        
                              "p__Firmicutes" = "#ffad12ff", 
                              "p__Gemmatimonadetes" = "#d96d3bff",
                              "p__Nitrospirae" = "#91569aff", 
                              "p__Planctomycetes" = "#5a9d5aff", 
                              "p__Proteobacteria" = "#419486ff", 
                              "p__Verrucomicrobia" = "#66628dff",
                               "Low_abundance" = "#999999ff")) +
       geom_bar(stat="identity", color="black", width=0.8) + scale_y_continuous(expand = c(0 ,0)) + 
       theme_light(base_size = 7) + theme(axis.text.y = element_text(size = 18),
                                          axis.text.x = element_text(size= 8, angle= 90)) 
ggsave("bact_phylum_plot.pdf", width=35, height=20, units="cm")

print("Most abundant Phylum")
phylum_plot
```

    [1] "Most abundant Phylum"



![png](output_14_1.png)



```R
#Shannon diversity by soil, rhizosphere and endosphere
div_tab_data_ns <- subset(div_tab_data, type!= "substrate")
div_tab_data_ns$type <- factor(div_tab_data_ns$type, levels=c("soil", "rhizosphere", "endosphere"))

shannon_type <- ggplot(div_tab_data_ns, aes(type, Shannon, fill=type)) + geom_boxplot() + geom_point(size = 1.5) + 
                scale_fill_manual(values = c("#ff5599ff", "#00aad4ff", "#afdde9ff")) + 
                stat_signif(test = t.test, map_signif_level = TRUE, comparisons = 
                            list(c("endosphere", "rhizosphere"), c("endosphere", "soil"), 
                                 c("rhizosphere", "soil"))) + 
                                theme_light() + theme(text = element_text(size=20), legend.position = "none") +
                                scale_x_discrete(labels=c("Soil", "Rhizosphere", "Endosphere")) 

ggsave("shannon_type.pdf", width = 20, height = 15, units = "cm")

print("Shannon diversity by sample type")
shannon_type
```

    [1] "Shannon diversity by sample type"



![png](output_15_1.png)



```R
#Plot relative abundance of Phylum by sample type

#Get relative abundances by sample type
Samples_toKeep <- c("soil", "rhizosphere", "endosphere")
phy_calabacita_ns <- subset_samples(phy_calabacita, type %in% Samples_toKeep)
calabacita_type <- merge_samples(phy_calabacita_ns, "type")
calabacita_type <- transform_sample_counts(calabacita_type, function(x) x / sum(x))

#Prepare table with information about taxonomy and relative abundances
calabacita_type_otu_phy <- t(otu_table(calabacita_type))
calabacita_type_tax_phy <- tax_table(calabacita_type)
calabacita_type_otu_tax <- cbind(calabacita_type_otu_phy, calabacita_type_tax_phy)
calabacita_type_otu_tax <- calabacita_type_otu_tax[, 1:5]
write.table(calabacita_type_otu_tax, "calabacita_type_otu_tax.tsv", sep = "\t")

calabacita_type_otu_tax <- read.table("calabacita_type_otu_tax.tsv", header=TRUE, row.names=1, 
                                 stringsAsFactors = FALSE, check.names=FALSE)
m_calabacita_type_otu_tax <- melt(calabacita_type_otu_tax, measure.vars = c("soil", "rhizosphere", "endosphere"))
colnames(m_calabacita_type_otu_tax) <- c("Kingdom", "Phylum", "Type", "Relative_abundance")

m_calabacita_type_otu_tax$Phylum[m_calabacita_type_otu_tax$Relative_abundance <= 0.01] <- "Low_abundance"

cm_calabacita_type_otu_tax <- aggregate(m_calabacita_type_otu_tax$Relative_abundance 
                                          ,by=list(m_calabacita_type_otu_tax$Phylum, 
                                           m_calabacita_type_otu_tax$Type),sum)
colnames(cm_calabacita_type_otu_tax) <- c("Phylum", "Type", "Relative_abundance")

#Create the plot in ggplot2
print("Most abundant Phylum")
phylum__type_plot <- ggplot(cm_calabacita_type_otu_tax, aes(x=Type, y=Relative_abundance,
                                          fill=fct_reorder(Phylum, Relative_abundance))) + 
               scale_fill_manual(values = c("p__Acidobacteria" = "#dd87b4ff", 
                              "p__Actinobacteria" = "#d26d7aff",
                              "p__Bacteroidetes" = "#b6742aff", 
                              "p__Chloroflexi" = "#f6ef32ff",                        
                              "p__Firmicutes" = "#ffad12ff", 
                              "p__Gemmatimonadetes" = "#d96d3bff",
                              "p__Nitrospirae" = "#91569aff", 
                              "p__Planctomycetes" = "#5a9d5aff", 
                              "p__Proteobacteria" = "#419486ff", 
                              "p__Verrucomicrobia" = "#66628dff",
                               "Low_abundance" = "#999999ff")) +
       geom_bar(stat="identity", color="black", width=0.8) + scale_y_continuous(expand = c(0 ,0)) + 
       theme_light(base_size = 7) + theme(axis.text.y = element_text(size = 20),
                                          axis.text.x = element_text(size= 17)) +
                                          scale_x_discrete(labels=c("Soil", "Rhizosphere", "Endosphere")) 
ggsave("phylum__type_plot.pdf", width=35, height=40, units="cm") 

phylum__type_plot
```

    [1] "Most abundant Phylum"



![png](output_16_1.png)



```R
#Shannon diversity by climate and sample type
div_tab_data_nz <- subset(div_tab_data, treatment!= "sterilized")
div_tab_data_nz <- subset(div_tab_data_nz, treatment!= "substrate")
div_tab_data_nz$samplec <- factor(div_tab_data_nz$samplec, levels=c("soil_humid", "soil_arid",
                            "rhizosphere_humid", "rhizosphere_arid", "endosphere_humid", "endosphere_arid"))
div_tab_data_nz$type <- factor(div_tab_data_nz$type, levels=c("soil", "rhizosphere", "endosphere"))

#t-test for significance between groups
print("t-test paired comparisons")
shannon_type_clime <- ggplot(div_tab_data_nz, aes(samplec, Shannon, fill=climate)) + geom_boxplot() + 
                geom_point(size = 1) +  scale_fill_manual(values = c("#ff9a58ff", "#5c5cafff")) + ylim(4,9) + 
                                theme_bw() + theme(text = element_text(size=20)) + 
                                theme(axis.text.x = element_text(angle = 90))                

s_sig <- shannon_type_clime + stat_signif(test = t.test, map_signif_level = TRUE, comparisons = 
                                list(c("soil_humid", "soil_arid"), c("rhizosphere_humid", "rhizosphere_arid"), 
                                 c("endosphere_humid", "endosphere_arid")))
s_sig

#Facet plot by sample type
s_facet <- shannon_type_clime + facet_wrap(~type, scale="free") + theme(axis.text.x = element_text(angle = 0))  +
              scale_x_discrete(labels=c("Humid", "Arid")) +
              theme(strip.background = element_rect(fill="white", color="black")) +
              theme(strip.text = element_text(colour="black", size=18.5))



ggsave("shannon_facet_type.pdf", width = 30, height = 10, units = "cm")

s_facet
print("Facet plot by sample type")

                                
```

    [1] "t-test paired comparisons"



![png](output_17_1.png)


    [1] "Facet plot by sample type"



![png](output_17_3.png)



```R
#Plot relative abundance of Phylum by sample type and climate

#Get relative abundances by sample type and climate
Samples_toKeep <- c("local", "inoculum", "fertilized")
phy_calabacita_nst <- subset_samples(phy_calabacita, treatment %in% Samples_toKeep)
calabacita_samplec <- merge_samples(phy_calabacita_nst, "samplec")
calabacita_samplec <- transform_sample_counts(calabacita_samplec, function(x) x / sum(x))

#Prepare table with information about taxonomy and relative abundances
calabacita_samplec_otu_phy <- t(otu_table(calabacita_samplec))
calabacita_samplec_tax_phy <- tax_table(calabacita_samplec)
calabacita_samplec_otu_tax <- cbind(calabacita_samplec_otu_phy, calabacita_samplec_tax_phy)
calabacita_samplec_otu_tax <- calabacita_samplec_otu_tax[, 1:8]

write.table(calabacita_samplec_otu_tax, "calabacita_samplec_otu_tax.tsv", sep = "\t")

calabacita_samplec_otu_tax <- read.table("calabacita_samplec_otu_tax.tsv", header=TRUE, row.names=1, 
                                 stringsAsFactors = FALSE, check.names=FALSE)
                                              
m_calabacita_samplec_otu_tax <- melt(calabacita_samplec_otu_tax, 
                                     measure.vars = c("soil_humid",  "soil_arid",                               
                                                     "rhizosphere_humid", "rhizosphere_arid",                                                  
                                                     "endosphere_humid", "endosphere_arid"))
                                              
colnames(m_calabacita_samplec_otu_tax) <- c("Kingdom", "Phylum", "Type_climate", "Relative_abundance")
m_calabacita_samplec_otu_tax$Phylum[m_calabacita_samplec_otu_tax$Relative_abundance <= 0.01] <- "Low_abundance"
                                              
cm_calabacita_samplec_otu_tax <- aggregate(m_calabacita_samplec_otu_tax$Relative_abundance 
                                          ,by=list(m_calabacita_samplec_otu_tax$Phylum, 
                                           m_calabacita_samplec_otu_tax$Type),sum)
                                              
#Add information about sample type before making facet plot in ggplot2                                             
colnames(cm_calabacita_samplec_otu_tax) <- c("Phylum", "Type_climate", "Relative_abundance")
cm_calabacita_samplec_otu_tax$Type <- factor(c(rep(c("soil"), times=21), rep(c("rhizosphere"), 
                                      times=20), rep(c("endosphere"), times=19)), 
                                            levels= c("soil", "rhizosphere", "endosphere"))

#Create the plot in ggplot2
print("Most abundant Phylum, by sample type, climate and treatment")
phylum__type_clime_plot <- ggplot(cm_calabacita_samplec_otu_tax, aes(x=Type_climate, y=Relative_abundance,
                                          fill=fct_reorder(Phylum, Relative_abundance))) + 
       scale_fill_manual(values = c("p__Acidobacteria" = "#dd87b4ff", 
                              "p__Actinobacteria" = "#d26d7aff",
                              "p__Bacteroidetes" = "#b6742aff", 
                              "p__Chloroflexi" = "#f6ef32ff",                        
                              "p__Firmicutes" = "#ffad12ff", 
                              "p__Gemmatimonadetes" = "#d96d3bff",
                              "p__Nitrospirae" = "#91569aff", 
                              "p__Planctomycetes" = "#5a9d5aff", 
                              "p__Proteobacteria" = "#419486ff", 
                              "p__Verrucomicrobia" = "#66628dff",
                               "Low_abundance" = "#999999ff")) +
       geom_bar(stat="identity", color="black", width=0.8) + scale_y_continuous(expand = c(0 ,0)) + 
       theme_light(base_size = 7) + theme(axis.text.y = element_text(size = 20),
                                          axis.text.x = element_text(size = 15)) + 
       facet_grid(~Type, scale="free") + scale_x_discrete(labels=c("Humid", "Arid"))  +
       theme(strip.background = element_rect(fill="white", color="black")) +
       theme(strip.text = element_text(colour="black", size=18.5))                      

ggsave("phylum__type_clime_plot.pdf", width=35, height=40, units="cm") 

phylum__type_clime_plot
```

    [1] "Most abundant Phylum, by sample type, climate and treatment"



![png](output_18_1.png)



```R
#Test statistically signifficant differences in relative abundance of phylum by climate

#Get relative abundances by sample type and climate
Samples_toKeep <- c("local", "inoculum", "fertilized")
phy_calabacita_nst <- subset_samples(phy_calabacita, treatment %in% Samples_toKeep)

#Join taxonomy table, otu table and metadata
#Get count table for the abundance of each genus
taxa_phy <- as.data.frame(as(tax_table(phy_calabacita_nst), "matrix"),row.names = FALSE)
otu_phy <- as.data.frame(as(t(otu_table(phy_calabacita_nst)), "matrix"),row.names = FALSE)
colnames(otu_phy) <- as.vector(taxa_phy$Phylum)

#Load climatic metadata
meta_phy_nst <- as.data.frame(meta(phy_calabacita_nst))
type <- meta_phy_nst$type
clime <- meta_phy_nst$clime
phy_nst_tab <- cbind(type, clime, otu_phy)

mphy_nst_tab <- melt(phy_nst_tab)
colnames(mphy_nst_tab) <- c("Sample_type", "Climate", "Phylum", "Relative_abundance")
head(mphy_nst_tab)
```

    Using type, clime as id variables



<table>
<thead><tr><th scope=col>Sample_type</th><th scope=col>Climate</th><th scope=col>Phylum</th><th scope=col>Relative_abundance</th></tr></thead>
<tbody>
	<tr><td>endosphere        </td><td>arid              </td><td>p__Armatimonadetes</td><td>0.0001374021      </td></tr>
	<tr><td>rhizosphere       </td><td>arid              </td><td>p__Armatimonadetes</td><td>0.0001113076      </td></tr>
	<tr><td>endosphere        </td><td>arid              </td><td>p__Armatimonadetes</td><td>0.0002922662      </td></tr>
	<tr><td>rhizosphere       </td><td>arid              </td><td>p__Armatimonadetes</td><td>0.0002233619      </td></tr>
	<tr><td>rhizosphere       </td><td>arid              </td><td>p__Armatimonadetes</td><td>0.0012776423      </td></tr>
	<tr><td>endosphere        </td><td>arid              </td><td>p__Armatimonadetes</td><td>0.0007838711      </td></tr>
</tbody>
</table>




```R
#Test statistically signifficant differences in relative abundance of phylum by climate

#Test effect of climate in relative abundance at phylum level in soil samples
print("t-test  in soil samples")
soil <- filter(mphy_nst_tab, Sample_type == "soil")

print("Firmicutes")
t.test(Relative_abundance~Climate, filter(soil, Phylum == "p__Firmicutes"))

print("Nitrospirae")
t.test(Relative_abundance~Climate, filter(soil, Phylum == "p__Nitrospirae"))

print("Planctomycetes")
t.test(Relative_abundance~Climate, filter(soil, Phylum == "p__Planctomycetes"))

print("Verrucomicrobia")
t.test(Relative_abundance~Climate, filter(soil, Phylum == "p__Verrucomicrobia"))

print("Gemmatimonadetes")
t.test(Relative_abundance~Climate, filter(soil, Phylum == "p__Gemmatimonadetes"))

print("Chloroflexi")
t.test(Relative_abundance~Climate, filter(soil, Phylum == "p__Chloroflexi"))

print("Acidobacteria")
t.test(Relative_abundance~Climate, filter(soil, Phylum == "p__Acidobacteria"))

print("Bacteroidetes")
t.test(Relative_abundance~Climate, filter(soil, Phylum == "p__Bacteroidetes"))

print("Actinobacteria")
t.test(Relative_abundance~Climate, filter(soil, Phylum == "p__Actinobacteria"))

print("Proteobacteria")
t.test(Relative_abundance~Climate, filter(soil, Phylum == "p__Proteobacteria"))
```

    [1] "t-test  in soil samples"
    [1] "Firmicutes"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = 2.3051, df = 2.1514, p-value = 0.1387
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.01449852  0.05336122
    sample estimates:
     mean in group arid mean in group humid 
            0.029119796         0.009688444 



    [1] "Nitrospirae"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = -0.67212, df = 1.0031, p-value = 0.623
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.08245622  0.07411306
    sample estimates:
     mean in group arid mean in group humid 
             0.01193509          0.01610667 



    [1] "Planctomycetes"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = 0.0045196, df = 2.1587, p-value = 0.9968
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.02542383  0.02548116
    sample estimates:
     mean in group arid mean in group humid 
             0.04440562          0.04437696 



    [1] "Verrucomicrobia"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = -3.6527, df = 1.0575, p-value = 0.1594
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.07872455  0.03996158
    sample estimates:
     mean in group arid mean in group humid 
             0.01402155          0.03340304 



    [1] "Gemmatimonadetes"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = -0.25965, df = 1.3753, p-value = 0.8286
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.07998800  0.07414023
    sample estimates:
     mean in group arid mean in group humid 
             0.05715907          0.06008295 



    [1] "Chloroflexi"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = 0.77732, df = 2.9537, p-value = 0.4944
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.03380317  0.05539802
    sample estimates:
     mean in group arid mean in group humid 
             0.06736229          0.05656486 



    [1] "Acidobacteria"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = -1.3344, df = 1.2922, p-value = 0.3699
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.2448517  0.1715292
    sample estimates:
     mean in group arid mean in group humid 
              0.1518937           0.1885549 



    [1] "Bacteroidetes"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = -0.6031, df = 1.3158, p-value = 0.6337
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.07578646  0.06429252
    sample estimates:
     mean in group arid mean in group humid 
             0.03991687          0.04566384 



    [1] "Actinobacteria"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = 5.3733, df = 2.6576, p-value = 0.0171
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     0.02818294 0.12747165
    sample estimates:
     mean in group arid mean in group humid 
              0.2358889           0.1580616 



    [1] "Proteobacteria"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = -0.70968, df = 1.927, p-value = 0.5539
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.2764938  0.2006199
    sample estimates:
     mean in group arid mean in group humid 
              0.3275289           0.3654659 




```R
#Test statistically signifficant differences in relative abundance of phylum by climate

#Test effect of climate in relative abundance at phylum level in rhizosphere samples
print("t-test  in rhizosphere samples")
rhizosphere <- filter(mphy_nst_tab, Sample_type == "rhizosphere" )

print("Firmicutes")
t.test(Relative_abundance~Climate, filter(rhizosphere, Phylum == "p__Firmicutes"))

print("Nitrospirae")
t.test(Relative_abundance~Climate, filter(rhizosphere, Phylum == "p__Nitrospirae"))

print("Planctomycetes")
t.test(Relative_abundance~Climate, filter(rhizosphere, Phylum == "p__Planctomycetes"))

print("Verrucomicrobia")
t.test(Relative_abundance~Climate, filter(rhizosphere, Phylum == "p__Verrucomicrobia"))

print("Gemmatimonadetes")
t.test(Relative_abundance~Climate, filter(rhizosphere, Phylum == "p__Gemmatimonadetes"))

print("Chloroflexi")
t.test(Relative_abundance~Climate, filter(rhizosphere, Phylum == "p__Chloroflexi"))

print("Acidobacteria")
t.test(Relative_abundance~Climate, filter(rhizosphere, Phylum == "p__Acidobacteria"))

print("Bacteroidetes")
t.test(Relative_abundance~Climate, filter(rhizosphere, Phylum == "p__Bacteroidetes"))

print("Actinobacteria")
t.test(Relative_abundance~Climate, filter(rhizosphere, Phylum == "p__Actinobacteria"))

print("Proteobacteria")
t.test(Relative_abundance~Climate, filter(rhizosphere, Phylum == "p__Proteobacteria"))
```

    [1] "t-test  in rhizosphere samples"
    [1] "Firmicutes"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = -0.081591, df = 6.3585, p-value = 0.9375
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.01577467  0.01474314
    sample estimates:
     mean in group arid mean in group humid 
             0.01221693          0.01273270 



    [1] "Nitrospirae"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = -1.3001, df = 6.1219, p-value = 0.2404
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.009865845  0.002998083
    sample estimates:
     mean in group arid mean in group humid 
            0.003595893         0.007029774 



    [1] "Planctomycetes"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = 0.18584, df = 11.274, p-value = 0.8559
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.01202086  0.01424524
    sample estimates:
     mean in group arid mean in group humid 
             0.02149373          0.02038154 



    [1] "Verrucomicrobia"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = -1.1753, df = 11.361, p-value = 0.2639
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.014218787  0.004294625
    sample estimates:
     mean in group arid mean in group humid 
             0.01797199          0.02293408 



    [1] "Gemmatimonadetes"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = -0.50636, df = 8.6875, p-value = 0.6252
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.01672952  0.01063736
    sample estimates:
     mean in group arid mean in group humid 
             0.02103093          0.02407701 



    [1] "Chloroflexi"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = -0.3554, df = 7.3132, p-value = 0.7323
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.02303400  0.01696923
    sample estimates:
     mean in group arid mean in group humid 
             0.02110669          0.02413907 



    [1] "Acidobacteria"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = -0.44686, df = 10.372, p-value = 0.6642
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.05383573  0.03577635
    sample estimates:
     mean in group arid mean in group humid 
             0.06843147          0.07746116 



    [1] "Bacteroidetes"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = 1.5439, df = 12.934, p-value = 0.1467
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.009804399  0.058823379
    sample estimates:
     mean in group arid mean in group humid 
             0.10222728          0.07771779 



    [1] "Actinobacteria"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = 1.3795, df = 11.323, p-value = 0.1944
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.01888958  0.08292969
    sample estimates:
     mean in group arid mean in group humid 
             0.11680097          0.08478092 



    [1] "Proteobacteria"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = -0.59686, df = 10.601, p-value = 0.5631
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.14952266  0.08595779
    sample estimates:
     mean in group arid mean in group humid 
              0.6067196           0.6385020 




```R
#Test statistically signifficant differences in relative abundance of phylum by climate

#Test effect of climate in relative abundance at phylum level in endosphere samples
print("t-test  in endosphere samples")
endosphere <- filter(mphy_nst_tab, Sample_type == "endosphere" )

print("Firmicutes")
t.test(Relative_abundance~Climate, filter(endosphere, Phylum == "p__Firmicutes"))

print("Nitrospirae")
t.test(Relative_abundance~Climate, filter(endosphere, Phylum == "p__Nitrospirae"))

print("Planctomycetes")
t.test(Relative_abundance~Climate, filter(endosphere, Phylum == "p__Planctomycetes"))

print("Verrucomicrobia")
t.test(Relative_abundance~Climate, filter(endosphere, Phylum == "p__Verrucomicrobia"))

print("Gemmatimonadetes")
t.test(Relative_abundance~Climate, filter(endosphere, Phylum == "p__Gemmatimonadetes"))

print("Chloroflexi")
t.test(Relative_abundance~Climate, filter(endosphere, Phylum == "p__Chloroflexi"))

print("Acidobacteria")
t.test(Relative_abundance~Climate, filter(endosphere, Phylum == "p__Acidobacteria"))

print("Bacteroidetes")
t.test(Relative_abundance~Climate, filter(endosphere, Phylum == "p__Bacteroidetes"))

print("Actinobacteria")
t.test(Relative_abundance~Climate, filter(endosphere, Phylum == "p__Actinobacteria"))

print("Proteobacteria")
t.test(Relative_abundance~Climate, filter(endosphere, Phylum == "p__Proteobacteria"))
```

    [1] "t-test  in endosphere samples"
    [1] "Firmicutes"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = 2.0175, df = 12.802, p-value = 0.06511
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.000419877  0.012004703
    sample estimates:
     mean in group arid mean in group humid 
            0.012693995         0.006901581 



    [1] "Nitrospirae"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = -0.75824, df = 6.3323, p-value = 0.4756
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.005570703  0.002909430
    sample estimates:
     mean in group arid mean in group humid 
            0.003192007         0.004522644 



    [1] "Planctomycetes"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = 2.0651, df = 12.539, p-value = 0.06023
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.0004582429  0.0187654499
    sample estimates:
     mean in group arid mean in group humid 
             0.01925186          0.01009826 



    [1] "Verrucomicrobia"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = 0.34962, df = 12.648, p-value = 0.7324
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.00997888  0.01381929
    sample estimates:
     mean in group arid mean in group humid 
             0.02287139          0.02095119 



    [1] "Gemmatimonadetes"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = 1.4053, df = 12.905, p-value = 0.1835
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.003292646  0.015523949
    sample estimates:
     mean in group arid mean in group humid 
             0.01734843          0.01123278 



    [1] "Chloroflexi"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = 1.5092, df = 12.878, p-value = 0.1554
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.002970883  0.016698900
    sample estimates:
     mean in group arid mean in group humid 
              0.0183799           0.0115159 



    [1] "Acidobacteria"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = 1.2537, df = 11.951, p-value = 0.2339
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.01898239  0.07037471
    sample estimates:
     mean in group arid mean in group humid 
             0.07199263          0.04629647 



    [1] "Bacteroidetes"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = 0.4269, df = 12.794, p-value = 0.6765
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.02460873  0.03670491
    sample estimates:
     mean in group arid mean in group humid 
             0.10052851          0.09448042 



    [1] "Actinobacteria"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = 1.9188, df = 12.003, p-value = 0.0791
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.004733615  0.074607429
    sample estimates:
     mean in group arid mean in group humid 
             0.08920273          0.05426583 



    [1] "Proteobacteria"



    
    	Welch Two Sample t-test
    
    data:  Relative_abundance by Climate
    t = -1.7546, df = 12.988, p-value = 0.1029
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.22121048  0.02293715
    sample estimates:
     mean in group arid mean in group humid 
              0.6347184           0.7338550 




```R
#Shannon diversity by climate, sample type, and local vs common garden experiment
div_tab_data_nz <- subset(div_tab_data, treatment!= "sterilized")
div_tab_data_nz <- subset(div_tab_data_nz, treatment!= "substrate")
div_tab_data_nz$lc_type_climate <- paste(div_tab_data_nz$local_commongarden,"_", div_tab_data_nz$samplec)
div_tab_data_nz$lc_type_climate <- factor(div_tab_data_nz$lc_type_climate, levels=c("L _ soil_humid", 
            "C _ soil_humid", "L _ soil_arid", "C _ soil_arid", "L _ rhizosphere_humid", "C _ rhizosphere_humid",
            "L _ rhizosphere_arid", "C _ rhizosphere_arid", "L _ endosphere_humid", "C _ endosphere_humid",
            "L _ endosphere_arid", "C _ endosphere_arid"))

#t-test for significance between groups
print("t-test paired comparisons")
shannon_type_clime_lc <- ggplot(div_tab_data_nz, aes(lc_type_climate, Shannon, fill=climate)) + geom_boxplot() + 
                geom_point(size = 1) +  scale_fill_manual(values = c("#ff9a58ff", "#5c5cafff")) + ylim(4,9) + 
                                theme_bw() + theme(text = element_text(size=20)) + 
                                theme(axis.text.x = element_text(angle = 90)) 

s_sig_lc <- shannon_type_clime_lc + stat_signif(test = t.test, map_signif_level = TRUE, comparisons = 
                                list(c("L _ rhizosphere_humid", "C _ rhizosphere_humid"), 
                                     c("L _ rhizosphere_arid", "C _ rhizosphere_arid"), 
                                     c("L _ endosphere_humid", "C _ endosphere_humid"),
                                     c("L _ endosphere_arid", "C _ endosphere_arid")))
s_sig_lc

#Facet plot by sample type
s_facet_lc <- s_sig_lc + facet_wrap(~type, scale="free") + theme(axis.text.x = element_text(angle = 0))  +
              theme(strip.background = element_rect(fill="white", color="black")) +
              theme(strip.text = element_text(colour="black", size=18.5)) +
              theme(axis.text.x = element_text(angle = 90))

ggsave("shannon_facet_type_lc.pdf", width = 30, height = 10, units = "cm")

s_facet_lc
print("Facet plot by sample type, local vs common garden")
```

    [1] "t-test paired comparisons"


    Warning message:
    “Computation failed in `stat_signif()`:
    missing value where TRUE/FALSE needed”Warning message:
    “Computation failed in `stat_signif()`:
    missing value where TRUE/FALSE needed”Warning message:
    “Computation failed in `stat_signif()`:
    missing value where TRUE/FALSE needed”Warning message:
    “Computation failed in `stat_signif()`:
    missing value where TRUE/FALSE needed”Warning message:
    “Computation failed in `stat_signif()`:
    missing value where TRUE/FALSE needed”Warning message:
    “Computation failed in `stat_signif()`:
    missing value where TRUE/FALSE needed”


![png](output_23_2.png)


    [1] "Facet plot by sample type, local vs common garden"



![png](output_23_4.png)



```R
#Estimmate diversity for soil and plant treatment: local, inoculum, fertilized inoculum, sterilized, substrate

#Subset soil samples
soil_div <- div_tab_data %>% subset(type == "soil")
soil_div <- cbind(as.data.frame(soil_div$type), as.data.frame(soil_div$Shannon))
colnames(soil_div) <- c("Treatment", "Shannon")

#Subset local, inoculum, fertilized, sterilized
trat_div <- div_tab_data %>% subset(locality != "c") %>% subset(type != "soil") %>% subset(treatment != "substrate")
trat_div <- cbind(as.data.frame(trat_div$treatment), as.data.frame(trat_div$Shannon))
colnames(trat_div) <- c("Treatment", "Shannon")

#Subset substrate samples
substrate_div <- div_tab_data %>% subset(soil == "substrate")
substrate_div <- cbind(as.data.frame(substrate_div$soil), as.data.frame(substrate_div$Shannon))
colnames(substrate_div) <- c("Treatment", "Shannon")

soil_trat_div <- rbind(soil_div, trat_div, substrate_div) 

treatment <- factor(soil_trat_div$Treatment, levels= c ("soil", "local", "inoculum", "fertilized", 
                                                   "sterilized", "substrate"))

shannon_trat <- ggplot(soil_trat_div, aes(treatment, Shannon, fill=treatment)) + geom_boxplot() + geom_point(size = 1.5) + 
                    scale_fill_manual(values = c("#aaccffff", "#aaccffff", "#ccaaffff", "#ccaaffff",
                                                "#e6e6e6ff", "#e6e6e6ff")) +
                                                 theme_light() + 
                                                 theme(text = element_text(size=15)) 

ggsave("shannon_trat.pdf", width = 20, height = 15, units = "cm")

print("Shannon diversity by treatment")
shannon_trat
```

    [1] "Shannon diversity by treatment"



![png](output_24_1.png)



```R
#Test statistically differences of groups via ANOVA
library(multcompView) 

#ANOVA test for soil, local root, inoculum, fertilized, sterilized and substrate
model <- lm(soil_trat_div$Shannon ~ soil_trat_div$Treatment)
ANOVA <- aov(model)

#Make Tukey test
TUKEY <- TukeyHSD(ANOVA, "soil_trat_div$Treatment", conf.level=0.99)

# I need to group the treatments that are not different each other together.
generate_label_df <- function(TUKEY, variable){
 
     # Extract labels and factor levels from Tukey post-hoc 
     Tukey.levels <- TUKEY[[variable]][,4]
     Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
     
     #I need to put the labels in the same order as in the boxplot :
     Tukey.labels$treatment=rownames(Tukey.labels)
     Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
     return(Tukey.labels)
     }
 
# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY , "soil_trat_div$Treatment")
LABELS
```


<table>
<thead><tr><th></th><th scope=col>Letters</th><th scope=col>treatment</th></tr></thead>
<tbody>
	<tr><th scope=row>fertilized</th><td>b         </td><td>fertilized</td></tr>
	<tr><th scope=row>inoculum</th><td>bc        </td><td>inoculum  </td></tr>
	<tr><th scope=row>local</th><td>b         </td><td>local     </td></tr>
	<tr><th scope=row>soil</th><td>c         </td><td>soil      </td></tr>
	<tr><th scope=row>sterilized</th><td>a         </td><td>sterilized</td></tr>
	<tr><th scope=row>substrate</th><td>a         </td><td>substrate </td></tr>
</tbody>
</table>




```R
#Plot relative abundance of Phylum of local, inocumul, fertilized and sterilized samples

#Get relative abundances of groups
phy_calabacita_substrate <- subset_samples(phy_calabacita, soil=="substrate")
phy_calabacita_substrate <- merge_samples(phy_calabacita_substrate, "treatment")

phy_calabacita_soil <- subset_samples(phy_calabacita, type=="soil")
phy_calabacita_soil <- merge_samples(phy_calabacita_soil, "type")

Samples_toKeep <- c("rhizosphere", "endosphere")
phy_calabacita_ns <- subset_samples(phy_calabacita, type %in% Samples_toKeep)
calabacita_samplec <- merge_samples(phy_calabacita_ns, "treatment")

phy_calabacita_complete <- merge_phyloseq(phy_calabacita_substrate, phy_calabacita_soil, calabacita_samplec)

calabacita_complete <- transform_sample_counts(phy_calabacita_complete, function(x) x / sum(x))

calabacita_complete_otu_phy <- t(otu_table(calabacita_complete))
calabacita_complete_tax_phy <- tax_table(calabacita_complete)
calabacita_complete_otu_tax <- cbind(calabacita_complete_otu_phy, calabacita_complete_tax_phy)
calabacita_complete_otu_tax <- calabacita_complete_otu_tax[, 1:8]

write.table(calabacita_complete_otu_tax, "calabacita_complete_otu_tax.tsv", sep = "\t")

calabacita_complete_otu_tax <- read.table("calabacita_complete_otu_tax.tsv", header=TRUE, row.names=1, 
                                 stringsAsFactors = FALSE, check.names=FALSE)
                                              
m_calabacita_complete_otu_tax <- melt(calabacita_complete_otu_tax, 
                                     measure.vars = c("soil", "local", "inoculum", "fertilized", 
                                                      "sterilized", "substrate"))

colnames(m_calabacita_complete_otu_tax) <- c("Kingdom", "Phylum", "Sample_group", "Relative_abundance")
m_calabacita_complete_otu_tax$Phylum[m_calabacita_complete_otu_tax$Relative_abundance <= 0.01] <- "Low_abundance"
                                               
cm_calabacita_complete_otu_tax <- aggregate(m_calabacita_complete_otu_tax$Relative_abundance 
                                          ,by=list(m_calabacita_complete_otu_tax$Phylum, 
                                           m_calabacita_complete_otu_tax$Sample_group),sum)
colnames(cm_calabacita_complete_otu_tax) <- c("Phylum", "Sample_group", "Relative_abundance")
                                              
#Create the plot
print("Most abundant Phylum of local, inoculum, fertilized and sterilized samples")
phylum_local_cg_plot <- ggplot(cm_calabacita_complete_otu_tax, aes(x=Sample_group, y=Relative_abundance,
                                          fill=fct_reorder(Phylum, Relative_abundance))) + 
        scale_fill_manual(values = c("p__Acidobacteria" = "#dd87b4ff", 
                              "p__Actinobacteria" = "#d26d7aff",
                              "p__Bacteroidetes" = "#b6742aff", 
                              "p__Chloroflexi" = "#f6ef32ff",                        
                              "p__Firmicutes" = "#ffad12ff", 
                              "p__Gemmatimonadetes" = "#d96d3bff",
                              "p__Nitrospirae" = "#91569aff", 
                              "p__Planctomycetes" = "#5a9d5aff", 
                              "p__Proteobacteria" = "#419486ff", 
                              "p__Verrucomicrobia" = "#66628dff",
                               "Low_abundance" = "#999999ff")) +
       geom_bar(stat="identity", color="black", width=0.8) + scale_y_continuous(expand = c(0 ,0)) + 
       theme_light(base_size = 7) + theme(axis.text.y = element_text(size = 18),
                                          axis.text.x = element_text(size = 15, angle = 90))

phylum_local_cg_plot 
```

    [1] "Most abundant Phylum of local, inoculum, fertilized and sterilized samples"



![png](output_26_1.png)



```R
#Estimate unifrac distance

u_dist <- distance(calabacita, "unifrac")
```

    Warning message in UniFrac(physeq, ...):
    “Randomly assigning root as -- 678871 -- in the phylogenetic tree in the data you provided.”


```R
# Whole samples MDS
tcalabacitau.mds <- ordinate(calabacita, "MDS", "u_dist")
ordum <- plot_ordination(calabacita, tcalabacitau.mds, type="sample", color="climetrat", shape="type") + 
         geom_text_repel(aes(label = name), size = 5, color = "gray") + geom_point(size=5) + theme_light() + 
         theme(text = element_text(size=15), legend.position = "bottom") + 
         scale_colour_manual(values = c("#ff9955ff", "#ffccaaff", "#8181c0ff", "#afc6e9ff", "black"))
ggsave("ordum.pdf", width = 30, height = 25, units = "cm")
ordum
```


![png](output_28_0.png)



```R
#PERMANOVA analysis by treatment
adonis(u_dist ~sampledata$treatment)
```


    
    Call:
    adonis(formula = u_dist ~ sampledata$treatment) 
    
    Permutation: free
    Number of permutations: 999
    
    Terms added sequentially (first to last)
    
                         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    sampledata$treatment  4     3.493 0.87326  3.4477 0.24283  0.001 ***
    Residuals            43    10.892 0.25329         0.75717           
    Total                47    14.384                 1.00000           
    ---
    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



```R
#PERMANOVA analysis by climate, without sterilized controls

#Filter samples sterilized samples and no soil controls
Samples_toKeep <- c("local", "inoculum", "fertilized")
ccalabacita <- subset_samples(calabacita, treatment %in% Samples_toKeep)
c_data <- metadata %>% filter(treatment == "local" | treatment == "inoculum" | treatment == "fertilized")

#Run PERMANOVA
uc_dist <- distance(ccalabacita, "unifrac")
adonis(uc_dist ~c_data$climate)
```

    Warning message in UniFrac(physeq, ...):
    “Randomly assigning root as -- 735160 -- in the phylogenetic tree in the data you provided.”


    
    Call:
    adonis(formula = uc_dist ~ c_data$climate) 
    
    Permutation: free
    Number of permutations: 999
    
    Terms added sequentially (first to last)
    
                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    c_data$climate  1    0.8108 0.81083  3.2898 0.09065  0.001 ***
    Residuals      33    8.1333 0.24646         0.90935           
    Total          34    8.9441                 1.00000           
    ---
    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



```R
#Canonical Analysis of Principal coordinate (CAP) for local samples
lcalabacita <- subset_samples(calabacita, treatment=="local")
lcalabacita.cap <- ordinate(lcalabacita, "CAP", "unifrac", 
                            ~ type + ph + ai + toc + tn + pt + no3 + nh4 + hpo4 + cp +cn + np + mat + map + 
                            soil + climate + specie) 
lcalabacita.cap

cap_plot <- plot_ordination(lcalabacita, lcalabacita.cap, color="climetrat", 
                            axes =c(1,2)) + aes(shape= type) + geom_text_repel(aes(label = name), 
                            size = 5, color = "gray") + geom_point(aes(colour = climetrat), size=5) + 
                            scale_colour_manual(values = c("#ff9955ff", "#8181c0ff")) + theme_light() + 
                            theme(text = element_text(size=15), legend.position = "bottom")   

arrowmat <- vegan::scores(lcalabacita.cap, display="bp")
arrowdf <- data.frame(labels=rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1,
                 yend = CAP2,
                 x = 0,
                 y = 0,
                 shape = NULL,
                 color = black,
                 label = labels)

label_map <- aes(x = 1.3 * CAP1,
                 y = 1.3 * CAP2,
                 shape = NULL,
                 color = NULL,
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

capl_plot <- cap_plot  +
  geom_segment(
    mapping = arrow_map,
    size = .5,
    data = arrowdf,
    color = "black",
    arrow = arrowhead
  ) +
  geom_text(
    mapping = label_map,
    size = 4, 
    data = arrowdf,
    show.legend = FALSE
  )
ggsave("capl_plot.pdf", width = 30, height = 30, units = "cm")
capl_plot

anl <- anova(lcalabacita.cap, permutations=9999)
anl
```

    Warning message in UniFrac(physeq, ...):
    “Randomly assigning root as -- 638888 -- in the phylogenetic tree in the data you provided.”


    Call: capscale(formula = distance ~ type + ph + ai + toc + tn + pt +
    no3 + nh4 + hpo4 + cp + cn + np + mat + map + soil + climate + specie,
    data = data)
    
                  Inertia Proportion Rank
    Total          3.5991     1.0000     
    Constrained    2.2562     0.6269    7
    Unconstrained  1.3429     0.3731    7
    Inertia is squared Unknown distance 
    Some constraints were aliased because they were collinear (redundant)
    
    Eigenvalues for constrained axes:
      CAP1   CAP2   CAP3   CAP4   CAP5   CAP6   CAP7 
    0.5861 0.4482 0.3323 0.3087 0.2283 0.1835 0.1690 
    
    Eigenvalues for unconstrained axes:
       MDS1    MDS2    MDS3    MDS4    MDS5    MDS6    MDS7 
    0.25013 0.23173 0.22367 0.17229 0.16282 0.15649 0.14578 



    Warning message:
    “Ignoring unknown aesthetics: label”


<table>
<thead><tr><th></th><th scope=col>Df</th><th scope=col>SumOfSqs</th><th scope=col>F</th><th scope=col>Pr(&gt;F)</th></tr></thead>
<tbody>
	<tr><th scope=row>Model</th><td>7       </td><td>2.256184</td><td>1.680059</td><td>1e-04   </td></tr>
	<tr><th scope=row>Residual</th><td>7       </td><td>1.342920</td><td>      NA</td><td>   NA   </td></tr>
</tbody>
</table>




![png](output_31_4.png)



```R
# Subsampling for common garden experiment
cgcalabacita <- merge_phyloseq(subset_samples(calabacita, treatment=="inoculum"), 
                              subset_samples(calabacita, treatment=="fertilized"), 
                              subset_samples(calabacita, treatment=="sterilized"))
cgcalabacita.cap <- ordinate(cgcalabacita, "CAP", "unifrac", ~ treatment + carotenoids + chlorophyll + 
                              + abiomass + sla)
cgcalabacita.cap

cap_plotg <- plot_ordination(cgcalabacita, cgcalabacita.cap, color="climetrat", axes =c(1,2)) + 
             aes(shape= type) + geom_text_repel(aes(label = name), size = 5, color = "gray") + 
             geom_point(aes(colour = climetrat), size=5) + theme_light() +
             theme(text = element_text(size=15), legend.position = "bottom") + 
            scale_colour_manual(values = c("#ff9955ff", "#8181c0ff", "#ffccaaff", "#afc6e9ff", "black"))
arrowmat <- vegan::scores(cgcalabacita.cap, display="bp")

arrowdf <- data.frame(labels=rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1,
                 yend = CAP2,
                 x = 0,
                 y = 0,
                 shape = NULL,
                 color = black,
                 label = labels)

label_map <- aes(x = 1.3 * CAP1,
                 y = 1.3 * CAP2,
                 shape = NULL,
                 color = NULL,
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

cap_plotg <- cap_plotg  +
  geom_segment(
    mapping = arrow_map,
    size = .5,
    data = arrowdf,
    color = "black",
    arrow = arrowhead
  ) +
  geom_text(
    mapping = label_map,
    size = 3, 
    data = arrowdf,
    show.legend = FALSE
  )
ggsave("cap_plotg.pdf", width = 30, height = 25, units = "cm")
cap_plotg

anl <- anova(cgcalabacita.cap, permutations=9999)
anl
```

    Warning message in UniFrac(physeq, ...):
    “Randomly assigning root as -- 437224 -- in the phylogenetic tree in the data you provided.”


    Call: capscale(formula = distance ~ treatment + carotenoids +
    chlorophyll + +abiomass + sla, data = data)
    
                  Inertia Proportion Rank
    Total           9.361      1.000     
    Constrained     3.548      0.379    6
    Unconstrained   5.813      0.621   25
    Inertia is squared Unknown distance 
    
    Eigenvalues for constrained axes:
      CAP1   CAP2   CAP3   CAP4   CAP5   CAP6 
    1.8245 0.5384 0.3767 0.2971 0.2615 0.2496 
    
    Eigenvalues for unconstrained axes:
      MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
    0.6287 0.4802 0.4341 0.3421 0.3182 0.2708 0.2582 0.2416 
    (Showing 8 of 25 unconstrained eigenvalues)



    Warning message:
    “Ignoring unknown aesthetics: label”


<table>
<thead><tr><th></th><th scope=col>Df</th><th scope=col>SumOfSqs</th><th scope=col>F</th><th scope=col>Pr(&gt;F)</th></tr></thead>
<tbody>
	<tr><th scope=row>Model</th><td> 6      </td><td>3.547869</td><td>2.543209</td><td>1e-04   </td></tr>
	<tr><th scope=row>Residual</th><td>25      </td><td>5.812652</td><td>      NA</td><td>   NA   </td></tr>
</tbody>
</table>




![png](output_32_4.png)



```R
#Most abundant genera

#Agglomerate OTUs by genus
gcalabacita <- tax_glom(calabacita, taxrank="Genus")
```


```R
#Get number of bacterial genus (removing unclassified: "g_")
gen_tab <- as.data.frame(as((tax_table(gcalabacita)), "matrix"),row.names = FALSE)
gen_tab <- gen_tab %>% subset(Genus != "g__")
gen_list <- gen_tab$Genus
gen_list <- unique(gen_list)
length(gen_list)
```


604



```R
#Most abundant genera

#Agglomerate OTUs by genus
gcalabacita <- tax_glom(calabacita, taxrank="Genus")
#Genus taxonomic labels

top20 <- prune_taxa(names(sort(taxa_sums(gcalabacita),TRUE)[1:20]),gcalabacita)
top20heatmap <- plot_heatmap(top20, sample.label = "name",  
                             taxa.label = "Genus", low = "#ffffffff", high = "#0a091fff",  na.value = "white",
                             taxa.order  = names(sort(taxa_sums(gcalabacita))),
                             sample.order = c("16S_1_Ls", "16S_1_Lr", "16S_1_Le", "16S_1_Ir", "16S_1_Ie", 
                                              "16S_1_Fr", "16S_1_Fe", "16S_1_Zr", "16S_1_Ze", "16S_2_Ls", 
                                              "16S_2_Lr", "16S_2_Le", "16S_2_Ir", "16S_2_Ie", "16S_2_Fr", 
                                              "16S_2_Fe", "16S_2_Zr", "16S_2_Ze", "16S_3_Ls", "16S_3_Lr", 
                                              "16S_3_Le", "16S_3_Ir", "16S_3_Ie", "16S_3_Fr", "16S_3_Fe", 
                                              "16S_3_Zr", "16S_3_Ze", "16S_4_Ls", "16S_4_Lr", "16S_4_Le", 
                                              "16S_4_Ir", "16S_4_Ie", "16S_4_Fr", "16S_4_Fe", "16S_4_Zr", 
                                              "16S_4_Ze", "16S_5_Ls", "16S_5_Lr", "16S_5_Le", "16S_5_Ir", 
                                              "16S_5_Ie", "16S_5_Fr", "16S_5_Fe", "16S_5_Zr", "16S_5_Ze", 
                                              "16S_Cr", "16S_Ce", "16S_Cn")) +
                             theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

ggsave("most_abundant_genus_order.pdf", width=60, height=30, units="cm")
top20heatmap
```

    Warning message:
    “Transformation introduced infinite values in discrete y-axis”Warning message:
    “Transformation introduced infinite values in discrete y-axis”


![png](output_35_1.png)



```R
#Most abundant genera

#Family taxonomic labels

top20 <- prune_taxa(names(sort(taxa_sums(gcalabacita),TRUE)[1:20]),gcalabacita)
top20heatmap <- plot_heatmap(top20, sample.label = "name",  
                             taxa.label = "Family", low = "#ffffffff", high = "#0a091fff",  na.value = "white",
                             taxa.order  = names(sort(taxa_sums(gcalabacita))),
                             sample.order = c("16S_1_Ls", "16S_1_Lr", "16S_1_Le", "16S_1_Ir", "16S_1_Ie", 
                                              "16S_1_Fr", "16S_1_Fe", "16S_1_Zr", "16S_1_Ze", "16S_2_Ls", 
                                              "16S_2_Lr", "16S_2_Le", "16S_2_Ir", "16S_2_Ie", "16S_2_Fr", 
                                              "16S_2_Fe", "16S_2_Zr", "16S_2_Ze", "16S_3_Ls", "16S_3_Lr", 
                                              "16S_3_Le", "16S_3_Ir", "16S_3_Ie", "16S_3_Fr", "16S_3_Fe", 
                                              "16S_3_Zr", "16S_3_Ze", "16S_4_Ls", "16S_4_Lr", "16S_4_Le", 
                                              "16S_4_Ir", "16S_4_Ie", "16S_4_Fr", "16S_4_Fe", "16S_4_Zr", 
                                              "16S_4_Ze", "16S_5_Ls", "16S_5_Lr", "16S_5_Le", "16S_5_Ir", 
                                              "16S_5_Ie", "16S_5_Fr", "16S_5_Fe", "16S_5_Zr", "16S_5_Ze", 
                                              "16S_Cr", "16S_Ce", "16S_Cn")) +
                             theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

ggsave("most_abundant_genus_order.pdf", width=60, height=30, units="cm")
top20heatmap
```

    Warning message:
    “Transformation introduced infinite values in discrete y-axis”Warning message:
    “Transformation introduced infinite values in discrete y-axis”


![png](output_36_1.png)



```R
#Most abundant genera


#Order taxonomic labels

top20 <- prune_taxa(names(sort(taxa_sums(gcalabacita),TRUE)[1:20]),gcalabacita)
top20heatmap <- plot_heatmap(top20, sample.label = "name",  
                             taxa.label = "Order", low = "#ffffffff", high = "#0a091fff",  na.value = "white",
                             taxa.order  = names(sort(taxa_sums(gcalabacita))),
                             sample.order = c("16S_1_Ls", "16S_1_Lr", "16S_1_Le", "16S_1_Ir", "16S_1_Ie", 
                                              "16S_1_Fr", "16S_1_Fe", "16S_1_Zr", "16S_1_Ze", "16S_2_Ls", 
                                              "16S_2_Lr", "16S_2_Le", "16S_2_Ir", "16S_2_Ie", "16S_2_Fr", 
                                              "16S_2_Fe", "16S_2_Zr", "16S_2_Ze", "16S_3_Ls", "16S_3_Lr", 
                                              "16S_3_Le", "16S_3_Ir", "16S_3_Ie", "16S_3_Fr", "16S_3_Fe", 
                                              "16S_3_Zr", "16S_3_Ze", "16S_4_Ls", "16S_4_Lr", "16S_4_Le", 
                                              "16S_4_Ir", "16S_4_Ie", "16S_4_Fr", "16S_4_Fe", "16S_4_Zr", 
                                              "16S_4_Ze", "16S_5_Ls", "16S_5_Lr", "16S_5_Le", "16S_5_Ir", 
                                              "16S_5_Ie", "16S_5_Fr", "16S_5_Fe", "16S_5_Zr", "16S_5_Ze", 
                                              "16S_Cr", "16S_Ce", "16S_Cn")) +
                             theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

ggsave("most_abundant_genus_order.pdf", width=60, height=30, units="cm")
top20heatmap
```

    Warning message:
    “Transformation introduced infinite values in discrete y-axis”Warning message:
    “Transformation introduced infinite values in discrete y-axis”


![png](output_37_1.png)



```R
#Comparisons between arid and humid groups

phumid <- subset_samples(gcalabacita, samplec=="rhizosphere_humid" | samplec=="endosphere_humid")
parid <- subset_samples(gcalabacita, samplec=="rhizosphere_arid" | samplec=="endosphere_arid")

m_arid_humid <- merge_samples(merge_phyloseq(phumid, parid), "climate")
m_arid_humid <- filter_taxa(m_arid_humid, function (x) {sum(x > 0) > 0}, prune=TRUE)
arid_humid_d <- t(as.data.frame(as(otu_table(m_arid_humid), "matrix"),row.names = FALSE))
colnames(arid_humid_d) <- c("Humid", "Arid")
arid_humid_t <- as.data.frame(as(tax_table(m_arid_humid), "matrix"),row.names = FALSE)
arid_humid_mi <- cbind(arid_humid_t, arid_humid_d)
arid_humid_mi <- arid_humid_mi %>% select(Genus, Humid, Arid)
arid_humid_mi[arid_humid_mi > 0] = 1
arid_humid_dat <- as.data.frame(arid_humid_mi)
row.names(arid_humid_dat) <- NULL
arid_humid_dat <- arid_humid_dat %>% subset(Genus != "g__")
arid_humid_dat <- unique(arid_humid_dat)

print("Genus per hystoric climatic conditions")
head(arid_humid_dat, 10)

#Venn diagram
print("Venn diagram")
arid_humid_dat$Genus <- NULL
vennDiagram(vennCounts(arid_humid_dat), circle.col=c("#5c5cafff", "#ff9a58ff"))

#Get genus list from arid and humid climates

#Arid genera list
farid <- filter_taxa(parid, function (x) {sum(x > 0) > 0}, prune=TRUE)
aridt <- as.data.frame(as(tax_table(farid), "matrix"),row.names = FALSE)
aridt <- aridt %>% subset(Genus != "g__")
aridt <- unique(aridt)
a_genus <- as.vector(aridt$Genus)
write.table(a_genus, "a_genus.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Humid genera list
fhumid <- filter_taxa(phumid, function (x) {sum(x > 0) > 0}, prune=TRUE)
humidt <- as.data.frame(as(tax_table(fhumid), "matrix"),row.names = FALSE, quote = FALSE, col.names = FALSE)
humidt <- humidt %>% subset(Genus != "g__")
humidt <- unique(humidt)
h_genus <- as.vector(humidt$Genus)
write.table(h_genus, "h_genus.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Make Venn diagrams with online tool:
#Alternatively, venn diagrams can be made from previous lists (a_genus.txt and h_genus.txt) with online tool:

#http://bioinformatics.psb.ugent.be/webtools/Venn/
```

    Warning message in Ops.factor(left, right):
    “‘>’ not meaningful for factors”

    [1] "Genus per hystoric climatic conditions"



<table>
<thead><tr><th></th><th scope=col>Genus</th><th scope=col>Humid</th><th scope=col>Arid</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>g__Kerstersia       </td><td>1                   </td><td>1                   </td></tr>
	<tr><th scope=row>2</th><td>g__Pigmentiphaga    </td><td>1                   </td><td>1                   </td></tr>
	<tr><th scope=row>3</th><td>g__Tetrathiobacter  </td><td>0                   </td><td>1                   </td></tr>
	<tr><th scope=row>4</th><td>g__Massilia         </td><td>1                   </td><td>1                   </td></tr>
	<tr><th scope=row>5</th><td>g__Janthinobacterium</td><td>1                   </td><td>1                   </td></tr>
	<tr><th scope=row>6</th><td>g__Burkholderia     </td><td>1                   </td><td>1                   </td></tr>
	<tr><th scope=row>10</th><td>g__Herbaspirillum   </td><td>1                   </td><td>1                   </td></tr>
	<tr><th scope=row>11</th><td>g__Ralstonia        </td><td>1                   </td><td>1                   </td></tr>
	<tr><th scope=row>12</th><td>g__Paucimonas       </td><td>1                   </td><td>0                   </td></tr>
	<tr><th scope=row>13</th><td>g__Collimonas       </td><td>1                   </td><td>1                   </td></tr>
</tbody>
</table>



    [1] "Venn diagram"



![png](output_38_4.png)



```R
#Plot most abundant genus associated to plants exclusive from each climatic condition

#Top genus in arid plant-associated samples

#Get exclusive genus list from arid climate
ea_genus <- setdiff(a_genus, h_genus)

e_parid <- subset_taxa(parid, Genus %in% ea_genus)

print("Most abundant genus in arid associated samples")
top10arid <- prune_taxa(names(sort(taxa_sums(e_parid),TRUE)[1:10]), e_parid)
top10_arid_heatmap <- plot_heatmap(top10arid, sample.label = "name", trans=NULL,  
                             taxa.order  = names(sort(taxa_sums(top10arid))),
                             taxa.label = "Genus", low = "#ffffffff", high = "#0a091fff",  na.value = "white",
                             sample.order = c("16S_3_Lr", "16S_3_Le", "16S_3_Ir", "16S_3_Ie", "16S_3_Fr", 
                                              "16S_3_Fe", "16S_4_Lr", "16S_4_Le", "16S_4_Ir", "16S_4_Ie",
                                              "16S_4_Fr", "16S_4_Fe", "16S_5_Lr", "16S_5_Le", "16S_5_Ir", 
                                              "16S_5_Ie", "16S_5_Fr", "16S_5_Fe")) +
                             theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

ggsave("top10_arid_heatmap.pdf", width=15, height=10, units="cm")
top10_arid_heatmap


#Get exclusive genus list from arid climate
eh_genus <- setdiff(h_genus, a_genus)

e_phumid <- subset_taxa(phumid, Genus %in% eh_genus)

top10humid <- prune_taxa(names(sort(taxa_sums(e_phumid),TRUE)[1:10]), e_phumid)
top10_humid_heatmap <- plot_heatmap(top10humid, sample.label = "name", trans=NULL,  
                             taxa.order  = names(sort(taxa_sums(top10humid))),
                             taxa.label = "Genus", low = "#ffffffff", high = "#0a091fff",  na.value = "white",
                             sample.order = c("16S_1_Lr", "16S_1_Le", "16S_1_Ir", "16S_1_Ie", "16S_1_Fr",  
                                              "16S_1_Fe", "16S_2_Lr", "16S_2_Le", "16S_2_Ir", "16S_2_Ie",     
                                              "16S_2_Fr", "16S_2_Fe")) +
                             theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

ggsave("top10_humid_heatmap.pdf", width=14, height=10, units="cm")
top10_humid_heatmap

print("Most abundant genus in humid associated samples")
```

    [1] "Most abundant genus in arid associated samples"



![png](output_39_1.png)


    [1] "Most abundant genus in humid associated samples"



![png](output_39_3.png)



```R
# Make object of core genera associated to plants, with treshold of 0 and prevalence in 99% of samples. 
library(microbiome)

#Subset plant associated samples
pgcalabacita <- subset_samples(gcalabacita, type=="rhizosphere" | type=="endosphere")

pccore <- core(pgcalabacita, detection = 0, prevalence = 0.99)

print("Core genus")
plot_heatmap(pccore, , sample.label = "name",  
                             taxa.label = "Genus", low = "#ffffffff", high = "#0a091fff",  na.value = "white",
                             taxa.order  = names(sort(taxa_sums(gcalabacita))), 
                             sample.order = c("16S_1_Ls", "16S_1_Lr", "16S_1_Le", "16S_1_Ir", "16S_1_Ie", 
                                              "16S_1_Fr", "16S_1_Fe", "16S_1_Zr", "16S_1_Ze", "16S_2_Ls", 
                                              "16S_2_Lr", "16S_2_Le", "16S_2_Ir", "16S_2_Ie", "16S_2_Fr", 
                                              "16S_2_Fe", "16S_2_Zr", "16S_2_Ze", "16S_3_Ls", "16S_3_Lr", 
                                              "16S_3_Le", "16S_3_Ir", "16S_3_Ie", "16S_3_Fr", "16S_3_Fe", 
                                              "16S_3_Zr", "16S_3_Ze", "16S_4_Ls", "16S_4_Lr", "16S_4_Le", 
                                              "16S_4_Ir", "16S_4_Ie", "16S_4_Fr", "16S_4_Fe", "16S_4_Zr", 
                                              "16S_4_Ze", "16S_5_Ls", "16S_5_Lr", "16S_5_Le", "16S_5_Ir", 
                                              "16S_5_Ie", "16S_5_Fr", "16S_5_Fe", "16S_5_Zr", "16S_5_Ze", 
                                              "16S_Cr", "16S_Ce", "16S_Cn")) +
                             theme(axis.text.x = element_text(size=8), axis.text.y = element_text(size=8))

#Get list of core genera
lcore <- as.data.frame(as(tax_table(pccore), "matrix"),row.names = FALSE, quote = FALSE, col.names = FALSE)
lcore <- lcore %>% subset(Genus != "g__")
lcore <- unique(lcore)
lcore <- as.vector(lcore$Genus)
write.table(lcore, "core_list.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
```

    [1] "Core genus"



![png](output_40_1.png)



```R
#Enrichement analysis of humid vs arid soils
library(DESeq2)

alpha = 0.05
merge_phyloseq(subset_samples(calabacita, samplec=="soil_arid"), 
               subset_samples(calabacita,samplec=="soil_humid"))->sar.shu.phy
sar.shu.ds <- phyloseq_to_deseq2(sar.shu.phy, ~samplec)
sar.shu.ds<-DESeq(sar.shu.ds, test="Wald", fitType="local")
sar.shu.ds.res <- results(sar.shu.ds, cooksCutoff = FALSE)
sar.shu.ds.res
write.csv(sar.shu.ds.res, "sar.shu.ds.res.csv" )
sigtab.sar.shu<-sar.shu.ds.res[which(sar.shu.ds.res$padj < alpha), ]
sigtab.sar.shu<-cbind(as(sigtab.sar.shu, "data.frame"), as(tax_table(sar.shu.phy)[rownames(sigtab.sar.shu), ], 
                                                           "matrix"))
sigtab.sar.shu.x=tapply(sigtab.sar.shu$log2FoldChange, sigtab.sar.shu$Genus, function(x) max(x))
sigtab.sar.shu.x=sort(sigtab.sar.shu.x, TRUE)
sigtab.sar.shu$Genus = factor(as.character(sigtab.sar.shu$Genus), levels=names(sigtab.sar.shu.x))

#Plot significant enriched OTUs in soils by climate
print("Enriched OTUs in humid and arid soil samples")                       
write.csv(sigtab.sar.shu, "sar.shu.csv" )
sarhudes <- ggplot(sigtab.sar.shu, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
                        scale_color_manual(values = c("p__Acidobacteria" = "#dd87b4ff", 
                              "p__Actinobacteria" = "#d26d7aff",
                              "p__Armatimonadetes" = "#ffa10dff",  
                              "p__Bacteroidetes" = "#b6742aff", 
                              "p__Chloroflexi" = "#f6ef32ff",                        
                              "p__Fibrobacteres" = "#e1c62fff",
                              "p__Firmicutes" = "#ffad12ff", 
                              "p__Gemmatimonadetes" = "#d96d3bff",
                              "p__Nitrospirae" = "#91569aff", 
                              "p__Planctomycetes" = "#5a9d5aff", 
                              "p__Proteobacteria" = "#419486ff", 
                              "p__Verrucomicrobia" = "#66628dff")) +
                     geom_point(size=3.5, alpha=1)+theme_light() + 
                        theme(text = element_text(size=10),
                              axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
                        ggtitle("Humid soil vs arid soil ; a=0.05")

ggsave("sarhudes.pdf", width = 15, height = 15, units = "cm")   
                        
sarhudes
                        
```

    Loading required package: S4Vectors
    Loading required package: stats4
    Loading required package: BiocGenerics
    Loading required package: parallel
    
    Attaching package: ‘BiocGenerics’
    
    The following objects are masked from ‘package:parallel’:
    
        clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
        clusterExport, clusterMap, parApply, parCapply, parLapply,
        parLapplyLB, parRapply, parSapply, parSapplyLB
    
    The following object is masked from ‘package:limma’:
    
        plotMA
    
    The following objects are masked from ‘package:dplyr’:
    
        combine, intersect, setdiff, union
    
    The following object is masked from ‘package:gridExtra’:
    
        combine
    
    The following objects are masked from ‘package:stats’:
    
        IQR, mad, sd, var, xtabs
    
    The following objects are masked from ‘package:base’:
    
        anyDuplicated, append, as.data.frame, basename, cbind, colMeans,
        colnames, colSums, dirname, do.call, duplicated, eval, evalq,
        Filter, Find, get, grep, grepl, intersect, is.unsorted, lapply,
        lengths, Map, mapply, match, mget, order, paste, pmax, pmax.int,
        pmin, pmin.int, Position, rank, rbind, Reduce, rowMeans, rownames,
        rowSums, sapply, setdiff, sort, table, tapply, union, unique,
        unsplit, which, which.max, which.min
    
    
    Attaching package: ‘S4Vectors’
    
    The following objects are masked from ‘package:dplyr’:
    
        first, rename
    
    The following object is masked from ‘package:tidyr’:
    
        expand
    
    The following object is masked from ‘package:base’:
    
        expand.grid
    
    Loading required package: IRanges
    
    Attaching package: ‘IRanges’
    
    The following object is masked from ‘package:microbiome’:
    
        coverage
    
    The following objects are masked from ‘package:dplyr’:
    
        collapse, desc, slice
    
    The following object is masked from ‘package:purrr’:
    
        reduce
    
    The following object is masked from ‘package:phyloseq’:
    
        distance
    
    Loading required package: GenomicRanges
    Loading required package: GenomeInfoDb
    Loading required package: SummarizedExperiment
    Loading required package: Biobase
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    
    Attaching package: ‘Biobase’
    
    The following object is masked from ‘package:phyloseq’:
    
        sampleNames
    
    Loading required package: DelayedArray
    Loading required package: matrixStats
    
    Attaching package: ‘matrixStats’
    
    The following objects are masked from ‘package:Biobase’:
    
        anyMissing, rowMedians
    
    The following object is masked from ‘package:dplyr’:
    
        count
    
    Loading required package: BiocParallel
    
    Attaching package: ‘DelayedArray’
    
    The following objects are masked from ‘package:matrixStats’:
    
        colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges
    
    The following object is masked from ‘package:purrr’:
    
        simplify
    
    The following objects are masked from ‘package:base’:
    
        aperm, apply
    
    converting counts to integer mode
    estimating size factors
    estimating dispersions
    gene-wise dispersion estimates
    mean-dispersion relationship
    final dispersion estimates
    fitting model and testing



    log2 fold change (MLE): samplec soil humid vs soil arid 
    Wald test p-value: samplec soil humid vs soil arid 
    DataFrame with 135110 rows and 6 columns
                    baseMean      log2FoldChange            lfcSE
                   <numeric>           <numeric>        <numeric>
    722089                 0                  NA               NA
    798506                 0                  NA               NA
    799936                 0                  NA               NA
    745779 0.175373106194631    1.25315857557631 4.10071413948114
    735451 0.170056600083439   -1.12216431992545 4.08943246846449
    ...                  ...                 ...              ...
    708451 0.399796271082846   -1.95050737981725 4.30481834676439
    717197                 0                  NA               NA
    738217                 0                  NA               NA
    759193                 0                  NA               NA
    679576  14.3202656150045 -0.0199870638037591 1.15608816938041
                          stat            pvalue              padj
                     <numeric>         <numeric>         <numeric>
    722089                  NA                NA                NA
    798506                  NA                NA                NA
    799936                  NA                NA                NA
    745779   0.305595204384295 0.759912872258698                NA
    735451  -0.274405881153187 0.783772723909017                NA
    ...                    ...               ...               ...
    708451  -0.453098649629038  0.65047770768766                NA
    717197                  NA                NA                NA
    738217                  NA                NA                NA
    759193                  NA                NA                NA
    679576 -0.0172885289661521 0.986206436796122 0.993247010533022


    [1] "Enriched OTUs in humid and arid soil samples"



![png](output_41_3.png)



```R
#Enrichement analysis of humid vs arid rhizospheres
alpha = 0.05
merge_phyloseq(subset_samples(calabacita, samplec=="rhizosphere_arid"), 
               subset_samples(calabacita, samplec=="rhizosphere_humid"))->rar.rhu.phy
rar.rhu.ds <- phyloseq_to_deseq2(rar.rhu.phy, ~samplec)
rar.rhu.ds<-DESeq(rar.rhu.ds, test="Wald", fitType="local")
rar.rhu.ds.res <- results(rar.rhu.ds, cooksCutoff = FALSE)
rar.rhu.ds.res
sigtab.rar.rhu<-rar.rhu.ds.res[which(rar.rhu.ds.res$padj < alpha), ]
sigtab.rar.rhu<-cbind(as(sigtab.rar.rhu, "data.frame"), as(tax_table(rar.rhu.phy)[rownames(sigtab.rar.rhu), ], 
                                                           "matrix"))
sigtab.rar.rhu.x=tapply(sigtab.rar.rhu$log2FoldChange, sigtab.rar.rhu$Genus, function(x) max(x))
sigtab.rar.rhu.x=sort(sigtab.rar.rhu.x, TRUE)
sigtab.rar.rhu$Genus = factor(as.character(sigtab.rar.rhu$Genus), levels=names(sigtab.rar.rhu.x))
write.csv(sigtab.rar.rhu, "rar.rhu.csv" )

#Plot significant enriched OTUs in arid rhizospheres by climate
print("Enriched OTUs in humid and arid rhizosphere samples")
rarhudes <- ggplot(sigtab.rar.rhu, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
                        scale_color_manual(values = c("p__Acidobacteria" = "#dd87b4ff", 
                              "p__Actinobacteria" = "#d26d7aff",
                              "p__Armatimonadetes" = "#ffa10dff",  
                              "p__Bacteroidetes" = "#b6742aff", 
                              "p__Chlorobi" = "#e9ddafd9", 
                              "p__Chloroflexi" = "#f6ef32ff",
                              "p__Cyanobacteria" = "#800000d2", 
                              "p__Gemmatimonadetes" = "#d96d3bff",
                              "p__Nitrospirae" = "#91569aff", 
                              "p__Planctomycetes" = "#5a9d5aff", 
                              "p__Proteobacteria" = "#419486ff", 
                              "p__Verrucomicrobia" = "#66628dff",
                              "p__Fibrobacteres" = "#e1c62fff",
                              "p__Firmicutes" = "#ffad12ff", 
                              "p__TM7" = "#d96d3bff", 
                              "p__WS3" = "#e41a1cff")) +
                        geom_point(size=2, alpha=1)+theme_light() + 
                        theme(text = element_text(size=6),
                              axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
                        ggtitle("Humid rizosphere vs arid rizosphere ; a=0.05")
ggsave("rarhudes.pdf", width = 35, height = 15, units = "cm")       
rarhudes
```

    converting counts to integer mode
    estimating size factors
    estimating dispersions
    gene-wise dispersion estimates
    mean-dispersion relationship
    final dispersion estimates
    fitting model and testing
    -- replacing outliers and refitting for 643 genes
    -- DESeq argument 'minReplicatesForReplace' = 7 
    -- original counts are preserved in counts(dds)
    estimating dispersions
    fitting model and testing



    log2 fold change (MLE): samplec rhizosphere humid vs rhizosphere arid 
    Wald test p-value: samplec rhizosphere humid vs rhizosphere arid 
    DataFrame with 135110 rows and 6 columns
                     baseMean     log2FoldChange             lfcSE
                    <numeric>          <numeric>         <numeric>
    722089 0.0352990963344062 -0.388979545383551  1.07532094353269
    798506  0.039180192271192 -0.388979545385015   1.0753223926491
    799936  0.148621758057253  -1.13883660021557  3.12860274434898
    745779  0.039180192271192 -0.388979545385015   1.0753223926491
    735451  0.200000817984907 0.0919132861706444  2.52453926295215
    ...                   ...                ...               ...
    708451  0.507224276511389 -0.314611384552778  1.75908403938608
    717197  0.102425304886924 0.0919182234114283  1.08607571454846
    738217                  0                 NA                NA
    759193 0.0457640796727091  0.412516624459454   1.0753293072191
    679576   15.9131916503839  -0.69312004575135 0.596055043611367
                         stat            pvalue              padj
                    <numeric>         <numeric>         <numeric>
    722089 -0.361733441279083 0.717551234444787                NA
    798506 -0.361732953804438 0.717551598760733                NA
    799936 -0.364008055120639  0.71585199431572                NA
    745779 -0.361732953804438 0.717551598760733                NA
    735451 0.0364079448157058 0.970957079334006                NA
    ...                   ...               ...               ...
    708451 -0.178849547553496 0.858055839131222                NA
    717197 0.0846333475467166 0.932552886661113                NA
    738217                 NA                NA                NA
    759193  0.383618879993386 0.701260944546836                NA
    679576  -1.16284570222221 0.244892110745522 0.376277257265578


    [1] "Enriched OTUs in humid and arid rhizosphere samples"



![png](output_42_3.png)



```R
# Humid (humid + perhumid) vs arid endospheres
alpha = 0.05
merge_phyloseq(subset_samples(calabacita, samplec=="endosphere_arid"), 
               subset_samples(calabacita, samplec=="endosphere_humid"))->ear.ehu.phy
ear.ehu.ds <- phyloseq_to_deseq2(ear.ehu.phy, ~samplec)
ear.ehu.ds<-DESeq(ear.ehu.ds, test="Wald", fitType="local")
ear.ehu.ds.res <- results(ear.ehu.ds, cooksCutoff = FALSE)
ear.ehu.ds.res
sigtab.ear.ehu<-ear.ehu.ds.res[which(ear.ehu.ds.res$padj < alpha), ]
sigtab.ear.ehu<-cbind(as(sigtab.ear.ehu, "data.frame"), as(tax_table(ear.ehu.phy)[rownames(sigtab.ear.ehu), ], 
                                                           "matrix"))
sigtab.ear.ehu.x=tapply(sigtab.ear.ehu$log2FoldChange, sigtab.ear.ehu$Genus, function(x) max(x))
sigtab.ear.ehu.x=sort(sigtab.ear.ehu.x, TRUE)
sigtab.ear.ehu$Genus = factor(as.character(sigtab.ear.ehu$Genus), levels=names(sigtab.ear.ehu.x))
write.csv(sigtab.ear.ehu, "ear.ehu.csv" )
                        
#Plot significant enriched OTUs in arid endospheres by climate
print("Enriched OTUs in humid and arid endosphere samples")
earhudes <- ggplot(sigtab.ear.ehu, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
       scale_color_manual(values = c("p__Acidobacteria" = "#dd87b4ff", 
                              "p__Actinobacteria" = "#d26d7aff",
                              "p__Armatimonadetes" = "#ffa10dff",  
                              "p__Bacteroidetes" = "#b6742aff", 
                              "p__Chlorobi" = "#e9ddafd9", 
                              "p__Chloroflexi" = "#f6ef32ff",                        
                              "p__Firmicutes" = "#ffad12ff", 
                              "p__Gemmatimonadetes" = "#d96d3bff",
                              "p__Planctomycetes" = "#5a9d5aff", 
                              "p__Proteobacteria" = "#419486ff", 
                              "p__Verrucomicrobia" = "#66628dff")) +
                        geom_point(size=2.5, alpha=0.85) +
                        theme_light() + theme(text = element_text(size=10), 
                                              axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
                        ggtitle("Humid endosphere vs arid endosphere ; a=0.05")
ggsave("earhudes.pdf", width = 22, height = 15, units = "cm")       
earhudes
```

    converting counts to integer mode
    estimating size factors
    estimating dispersions
    gene-wise dispersion estimates
    mean-dispersion relationship
    final dispersion estimates
    fitting model and testing
    -- replacing outliers and refitting for 885 genes
    -- DESeq argument 'minReplicatesForReplace' = 7 
    -- original counts are preserved in counts(dds)
    estimating dispersions
    fitting model and testing



    log2 fold change (MLE): samplec endosphere humid vs endosphere arid 
    Wald test p-value: samplec endosphere humid vs endosphere arid 
    DataFrame with 135110 rows and 6 columns
                     baseMean     log2FoldChange             lfcSE
                    <numeric>          <numeric>         <numeric>
    722089  0.112003168732085  0.497634812047458  3.13503791676162
    798506 0.0722410243145451 -0.309755406231521  3.13470218377282
    799936 0.0361205121572725 0.0167411916298426  3.13503791676162
    745779                  0                 NA                NA
    735451  0.173069060019883   1.37474706950732   3.1294169368325
    ...                   ...                ...               ...
    708451  0.357860796445828   0.17682438619249   2.2585823538409
    717197 0.0361205121572725 0.0167411916298426  3.13503791676162
    738217  0.125377442978294   1.37474706950732   3.1294169368325
    759193 0.0744071901640073  0.818229298583408  3.13503791676162
    679576   14.6059830999779   -1.0875809140328 0.917401206834354
                          stat            pvalue              padj
                     <numeric>         <numeric>         <numeric>
    722089    0.15873326743094 0.873879026841014                NA
    798506 -0.0988149393696818 0.921285206646926                NA
    799936  0.0053400284380406 0.995739294004696                NA
    745779                  NA                NA                NA
    735451   0.439298149545645  0.66044551484066                NA
    ...                    ...               ...               ...
    708451  0.0782899883600818 0.937597381150745                NA
    717197  0.0053400284380406 0.995739294004696                NA
    738217   0.439298149545645  0.66044551484066                NA
    759193   0.260995024720023 0.794096344023989                NA
    679576   -1.18550194389397 0.235819046536132 0.513105790144318


    [1] "Enriched OTUs in humid and arid endosphere samples"



![png](output_43_3.png)



```R
#Infer correlations from genus overrepresented in arid and humid samples
library(corrplot)
library(psych)

#Genus from rhizosphere positevely related to plant phenotype
rhizo <- unique(sigtab.rar.rhu$Genus)
rhizo <- rhizo %>% subset(rhizo!= "g__")

rgcalabacita <- transform_sample_counts(gcalabacita, function(x) x / sum(x))
phylo_rhizo <- subset_samples(subset_samples(rgcalabacita, type=="rhizosphere"), treatment!="local")
phylo_rhizo <- subset_taxa(phylo_rhizo, Genus %in% rhizo)

                                        #Get metadata for correlations
meta_num <- unlist(lapply(meta(phylo_rhizo), is.numeric))  
metadata_rhizo <- meta(phylo_rhizo)[ , meta_num]
metadata_rhizo <- metadata_rhizo[, colSums(is.na(metadata_rhizo)) != nrow(metadata_rhizo)]
metadata_rhizo <- metadata_rhizo[, 1:16]

#Get count table for the abundance of each genus
taxa_rhizo <- as.data.frame(as(tax_table(phylo_rhizo), "matrix"),row.names = FALSE)
otu_rhizo <- as.data.frame(as(t(otu_table(phylo_rhizo)), "matrix"),row.names = FALSE)
colnames(otu_rhizo) <- as.vector(taxa_rhizo$Genus)

#Bind both tables
print("Phenotypic variables and relative abundance")
tab_rhizo <- cbind(metadata_rhizo, otu_rhizo)
write.csv(tab_rhizo, "tab_rhizo.csv")
tab_rhizo <- read.csv("tab_rhizo.csv", row.names=1)

#Correlogram for rizosphere sample matrix
#Correlation are done with spearman correlations, with BH correction for false discovery rate.
col <- colorRampPalette(c("#8b5016ff", "white", "#0b655dff"))(20)
rcorr <- corr.test(tab_rhizo, method = "spearman", adjust = "BH", alpha = 0.05)

r <- rcorr$r[1:13, 17:140]
p <- rcorr$p[1:13, 17:140]

corrplot(r, method= "color", tl.col = "black", p.mat = p, sig.level = .05, insig = "blank", 
         col = col, tl.cex=0.45)

pdf(file = "rcplot.pdf", width = 20, height = 15)
ecplot <- corrplot(r, method= "color", tl.col = "black", p.mat = p, sig.level = .05, insig = "blank", col = col)
dev.off()   

```

    corrplot 0.84 loaded
    
    Attaching package: ‘psych’
    
    The following object is masked from ‘package:IRanges’:
    
        reflect
    
    The following objects are masked from ‘package:ggplot2’:
    
        %+%, alpha
    


    [1] "Phenotypic variables and relative abundance"



<strong>png:</strong> 2



![png](output_44_3.png)



```R
#Genus from endosphere positevely related to plant phenotype
endo <- unique(sigtab.ear.ehu$Genus)
endo <- endo %>% subset(endo!= "g__")

phylo_endo <- subset_samples(subset_samples(rgcalabacita, type=="endosphere"), treatment!="local")
phylo_endo <- subset_taxa(phylo_endo, Genus %in% endo)

#Get metadata for correlations
meta_num <- unlist(lapply(meta(phylo_endo), is.numeric))  
metadata_endo <- meta(phylo_endo)[ , meta_num]
metadata_endo <- metadata_endo[, colSums(is.na(metadata_endo)) != nrow(metadata_endo)]
metadata_endo <- metadata_endo[, 1:16]

#Get count table for the abundance of each genus
taxa_endo <- as.data.frame(as(tax_table(phylo_endo), "matrix"),row.names = FALSE)
otu_endo <- as.data.frame(as(t(otu_table(phylo_endo)), "matrix"),row.names = FALSE)
colnames(otu_endo) <- as.vector(taxa_endo$Genus)

#Bind both tables
print("Phenotypic variables and relative abundance")
tab_endo <- cbind(metadata_endo, otu_endo)
write.csv(tab_endo, "tab_endo.csv")
tab_endo <- read.csv("tab_endo.csv", row.names=1)

#Correlogram for rizosphere sample matrix
#Correlation are done with spearman correlations, with BH correction for false discovery rate.
col <- colorRampPalette(c("#8b5016ff", "white", "#0b655dff"))(20)
rcorr <- corr.test(tab_endo, method = "spearman", adjust = "BH", alpha = 0.05)

r <- rcorr$r[1:13, 17:74]
p <- rcorr$p[1:13, 17:74]

corrplot(r, method= "color", tl.col = "black", p.mat = p, sig.level = .05, insig = "blank", 
         col = col, tl.cex=0.45)


pdf(file = "ecplot.pdf", width = 20, height = 15)
ecplot <- corrplot(r, method= "color", tl.col = "black", p.mat = p, sig.level = .05, insig = "blank", col = col)
dev.off()
```

    [1] "Phenotypic variables and relative abundance"



<strong>png:</strong> 2



![png](output_45_2.png)

