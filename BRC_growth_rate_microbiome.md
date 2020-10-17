How fast does biological rock crust grow?
================
17 October, 2020

  - [Setting general parameters:](#setting-general-parameters)
  - [Description](#description)
      - [Load data](#load-data)
      - [Inspect library size and number of
        OTU](#inspect-library-size-and-number-of-otu)
      - [Explore the prevalence of different taxa in the
        database](#explore-the-prevalence-of-different-taxa-in-the-database)
          - [Plot general prevalence features of the
            phyla](#plot-general-prevalence-features-of-the-phyla)
          - [Plot general prevalence features of the top 20
            orders](#plot-general-prevalence-features-of-the-top-20-orders)
          - [Unsupervised filtering by
            prevalence](#unsupervised-filtering-by-prevalence)
      - [Exploring the dataset
        features](#exploring-the-dataset-features)
      - [Alpha diversity](#alpha-diversity)
          - [Calculate and plot alpha diversity
            metrics](#calculate-and-plot-alpha-diversity-metrics)
          - [Test the differences in alpha
            diversity.](#test-the-differences-in-alpha-diversity.)
          - [Plot all alpha diversity metrics
            together](#plot-all-alpha-diversity-metrics-together)
      - [Beta diversity](#beta-diversity)
          - [Calculate ordinations](#calculate-ordinations)
          - [Test differences between samples on the phylum
            level](#test-differences-between-samples-on-the-phylum-level)
      - [Taxonmic distribution
        analysis](#taxonmic-distribution-analysis)
      - [Differential abundance models](#differential-abundance-models)
      - [Compose figures](#compose-figures)
      - [References](#references)

### Setting general parameters:

``` r
set.seed(15102020)
bootstraps <- 1000
data_path <- "./Data/"
Proj_name <- "BRC_growth_rate"

Browns <- RColorBrewer::brewer.pal(n = 9, "YlOrBr")[9:6]
Greens <- RColorBrewer::brewer.pal(n = 9, "YlGn")[9:6]
Blues <- RColorBrewer::brewer.pal(n = 9, "YlGnBu")[9:6]
Gradient.colours <- c(Browns[1], Greens[1], Browns[2], Greens[2], Browns[3], Greens[3], Browns[4], Greens[4], Blues)
```

## Description

This script reproduces all sequence analysis steps and plots included in
the paper plus some additional exploratory analyses.

### Load data

``` r
OTUmat <- t(read.csv(paste0(data_path, "Shivta_site_otuTab2.txt"), header = TRUE, row.names = 1))
sort.order <- as.numeric(gsub("OTU([0-9]+)", "\\1", colnames( OTUmat )))
OTUmat <- OTUmat[, order(sort.order )]

Metadata <- read.csv(paste0(data_path, "Shivta_metadata.csv"), row.names = 1, header = TRUE)

read_csv(paste0(data_path, "Shivta_metadata.csv"),
                     trim_ws = TRUE) %>%
  mutate_at(
    c(
      "Rock.type",
      "Location"
    ), 
    ~(factor(.))
  ) %>% 
  column_to_rownames("Sample.code") ->
  Metadata

row.names(OTUmat) <- gsub("(.*)Nimrod[0-9]+|Osnat[0-9]+", "\\1", row.names( OTUmat))
Metadata <- Metadata[order(row.names(Metadata)), ]
OTUmat <- OTUmat[order(row.names(OTUmat)), ]
# calculate sample size
Metadata$Library.size = rowSums(OTUmat)
Metadata$Location.rock <- with(Metadata, Location:Rock.type)

# Load taxonomy data
tax.file <- "Shivta_site_silva.nrv119.taxonomy"
Taxonomy <- read.table(paste0(data_path, tax.file), stringsAsFactors = FALSE) # read taxonomy file

# count how many ';' in each cell and add up to 6
for (i in 1:nrow(Taxonomy)){
  semicolons <- length(gregexpr(";", Taxonomy$V2[i] )[[1]])
  if (semicolons < 6){
    x <- paste0( rep("Unclassified;", 6 - semicolons ), collapse = "")
    Taxonomy$V2[i] <- paste0( Taxonomy$V2[i], x, sep = "")
  }
}
# split taxonomy to columns
do.call( "rbind", strsplit( Taxonomy$V1, ";", fixed = TRUE)) %>% 
  gsub( "size=([0-9]+)", "\\1", .) %>%
  data.frame( ., do.call( "rbind", strsplit( Taxonomy$V2, ";", fixed = TRUE)), stringsAsFactors = F) %>% 
  apply(., 2, function(x) gsub( "\\(.*\\)", "", x)) %>% 
  replace(., . == "unclassified", "Unclassified") -> 
  Taxonomy
colnames( Taxonomy ) <- c( "OTU", "Frequency", "Domain", "Phylum", "Class", "Order", "Family", "Genus" )
# rownames(Taxonomy) <- colnames(Rock_weathering_OTUmat)
rownames(Taxonomy) <- Taxonomy[, 1]

Tree_IQ <- read_tree(paste0(data_path, "Shivta_site_otuReps.filtered.align.treefile"))

# generate phyloseq object
Ps_obj <- phyloseq(otu_table(OTUmat, taxa_are_rows = FALSE),
                   tax_table(Taxonomy[, -c(1, 2)]),
                   sample_data(Metadata),
                   phy_tree(Tree_IQ)
)
# Reorder factors for plotting
sample_data(Ps_obj)$Location %<>% fct_relevel("Slope", "City")
```

Remove un- and mis-classified sequences, chloroplasts and mitochondria

``` r
domains2remove <- c("", "Eukaryota", "Unclassified")
classes2remove <- c("Chloroplast")
families2remove <- c("Mitochondria")
Ps_obj_filt <- subset_taxa(Ps_obj, !is.na(Phylum) &
                        !Domain %in% domains2remove &
                      !Class %in% classes2remove &
                      !Family %in% families2remove)
```

### Inspect library size and number of OTU

``` r
Ps_obj_df <-
  as.data.frame(sample_data(Ps_obj_filt)) # Put sample_data into a ggplot-friendly data.frame
Ps_obj_df <- Ps_obj_df[order(Ps_obj_df$Library.size), ]
Ps_obj_df$Index <- seq(nrow(Ps_obj_df))
ggplot(data = Ps_obj_df, 
       aes(x = Index, y = Library.size, color = Location.rock)) + 
  geom_point(size = 4) + 
  scale_colour_manual(values = ggpomological:::pomological_palette[c(2, 1, 9, 3)], name = "Location.rock")
```

![](BRC_growth_rate_figures/Library%20Sizes-1.svg)<!-- -->

``` r
summary(sample_sums(Ps_obj_filt))
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   23972   44989   54329   53254   61352   77463

``` r
summary(taxa_sums(Ps_obj_filt))
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##       3      70     191    1814     769   76180

### Explore the prevalence of different taxa in the database

``` r
prevdf <- apply(X = otu_table(Ps_obj_filt),
                 MARGIN = ifelse(taxa_are_rows(Ps_obj_filt), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(Ps_obj_filt),
                      tax_table(Ps_obj_filt))
prevdf %>%
  group_by(Phylum) %>%
  summarise(`Mean prevalence` = mean(Prevalence),
            `Sum prevalence` = sum(Prevalence)) ->
  Prevalence_phylum_summary

Prevalence_phylum_summary %>% 
  kable(., digits = c(0, 1, 0)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Phylum

</th>

<th style="text-align:right;">

Mean prevalence

</th>

<th style="text-align:right;">

Sum prevalence

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Acidobacteria

</td>

<td style="text-align:right;">

11.7

</td>

<td style="text-align:right;">

152

</td>

</tr>

<tr>

<td style="text-align:left;">

Actinobacteria

</td>

<td style="text-align:right;">

16.3

</td>

<td style="text-align:right;">

3414

</td>

</tr>

<tr>

<td style="text-align:left;">

Aquificae

</td>

<td style="text-align:right;">

4.0

</td>

<td style="text-align:right;">

12

</td>

</tr>

<tr>

<td style="text-align:left;">

Armatimonadetes

</td>

<td style="text-align:right;">

7.6

</td>

<td style="text-align:right;">

76

</td>

</tr>

<tr>

<td style="text-align:left;">

Bacteroidetes

</td>

<td style="text-align:right;">

11.8

</td>

<td style="text-align:right;">

1027

</td>

</tr>

<tr>

<td style="text-align:left;">

Candidate\_division\_KB1

</td>

<td style="text-align:right;">

8.0

</td>

<td style="text-align:right;">

8

</td>

</tr>

<tr>

<td style="text-align:left;">

Candidate\_division\_OD1

</td>

<td style="text-align:right;">

4.0

</td>

<td style="text-align:right;">

4

</td>

</tr>

<tr>

<td style="text-align:left;">

Candidate\_division\_OP11

</td>

<td style="text-align:right;">

4.5

</td>

<td style="text-align:right;">

9

</td>

</tr>

<tr>

<td style="text-align:left;">

Candidate\_division\_TM7

</td>

<td style="text-align:right;">

7.9

</td>

<td style="text-align:right;">

158

</td>

</tr>

<tr>

<td style="text-align:left;">

Chloroflexi

</td>

<td style="text-align:right;">

11.0

</td>

<td style="text-align:right;">

1124

</td>

</tr>

<tr>

<td style="text-align:left;">

Cyanobacteria

</td>

<td style="text-align:right;">

14.9

</td>

<td style="text-align:right;">

506

</td>

</tr>

<tr>

<td style="text-align:left;">

Deinococcus-Thermus

</td>

<td style="text-align:right;">

15.8

</td>

<td style="text-align:right;">

174

</td>

</tr>

<tr>

<td style="text-align:left;">

Euryarchaeota

</td>

<td style="text-align:right;">

3.9

</td>

<td style="text-align:right;">

90

</td>

</tr>

<tr>

<td style="text-align:left;">

Firmicutes

</td>

<td style="text-align:right;">

13.4

</td>

<td style="text-align:right;">

214

</td>

</tr>

<tr>

<td style="text-align:left;">

Gemmatimonadetes

</td>

<td style="text-align:right;">

11.2

</td>

<td style="text-align:right;">

617

</td>

</tr>

<tr>

<td style="text-align:left;">

Nitrospirae

</td>

<td style="text-align:right;">

6.0

</td>

<td style="text-align:right;">

12

</td>

</tr>

<tr>

<td style="text-align:left;">

Planctomycetes

</td>

<td style="text-align:right;">

10.8

</td>

<td style="text-align:right;">

140

</td>

</tr>

<tr>

<td style="text-align:left;">

Proteobacteria

</td>

<td style="text-align:right;">

14.4

</td>

<td style="text-align:right;">

1817

</td>

</tr>

<tr>

<td style="text-align:left;">

Tenericutes

</td>

<td style="text-align:right;">

3.0

</td>

<td style="text-align:right;">

3

</td>

</tr>

<tr>

<td style="text-align:left;">

Verrucomicrobia

</td>

<td style="text-align:right;">

23.0

</td>

<td style="text-align:right;">

23

</td>

</tr>

<tr>

<td style="text-align:left;">

WCHB1-60

</td>

<td style="text-align:right;">

7.0

</td>

<td style="text-align:right;">

28

</td>

</tr>

</tbody>

</table>

``` r
prevdf %>%
  group_by(Order) %>%
  summarise(`Mean prevalence` = mean(Prevalence),
            `Sum prevalence` = sum(Prevalence)) ->
  Prevalence_Order_summary

Prevalence_Order_summary %>% 
  kable(., digits = c(0, 1, 0)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Order

</th>

<th style="text-align:right;">

Mean prevalence

</th>

<th style="text-align:right;">

Sum prevalence

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

AKIW781

</td>

<td style="text-align:right;">

9.4

</td>

<td style="text-align:right;">

198

</td>

</tr>

<tr>

<td style="text-align:left;">

AKYG1722

</td>

<td style="text-align:right;">

10.7

</td>

<td style="text-align:right;">

75

</td>

</tr>

<tr>

<td style="text-align:left;">

AT425-EubC11\_terrestrial\_group

</td>

<td style="text-align:right;">

10.2

</td>

<td style="text-align:right;">

378

</td>

</tr>

<tr>

<td style="text-align:left;">

Acidimicrobiales

</td>

<td style="text-align:right;">

14.7

</td>

<td style="text-align:right;">

368

</td>

</tr>

<tr>

<td style="text-align:left;">

Alteromonadales

</td>

<td style="text-align:right;">

12.0

</td>

<td style="text-align:right;">

12

</td>

</tr>

<tr>

<td style="text-align:left;">

Aquificales

</td>

<td style="text-align:right;">

4.0

</td>

<td style="text-align:right;">

12

</td>

</tr>

<tr>

<td style="text-align:left;">

Ardenticatenales

</td>

<td style="text-align:right;">

8.0

</td>

<td style="text-align:right;">

16

</td>

</tr>

<tr>

<td style="text-align:left;">

B103G10

</td>

<td style="text-align:right;">

7.0

</td>

<td style="text-align:right;">

7

</td>

</tr>

<tr>

<td style="text-align:left;">

BD2-11\_terrestrial\_group

</td>

<td style="text-align:right;">

5.0

</td>

<td style="text-align:right;">

5

</td>

</tr>

<tr>

<td style="text-align:left;">

Bacillales

</td>

<td style="text-align:right;">

13.7

</td>

<td style="text-align:right;">

164

</td>

</tr>

<tr>

<td style="text-align:left;">

Bacteroidales

</td>

<td style="text-align:right;">

17.0

</td>

<td style="text-align:right;">

17

</td>

</tr>

<tr>

<td style="text-align:left;">

Bdellovibrionales

</td>

<td style="text-align:right;">

8.0

</td>

<td style="text-align:right;">

40

</td>

</tr>

<tr>

<td style="text-align:left;">

Burkholderiales

</td>

<td style="text-align:right;">

18.5

</td>

<td style="text-align:right;">

204

</td>

</tr>

<tr>

<td style="text-align:left;">

C0119

</td>

<td style="text-align:right;">

5.0

</td>

<td style="text-align:right;">

5

</td>

</tr>

<tr>

<td style="text-align:left;">

Caldilineales

</td>

<td style="text-align:right;">

18.0

</td>

<td style="text-align:right;">

18

</td>

</tr>

<tr>

<td style="text-align:left;">

Campylobacterales

</td>

<td style="text-align:right;">

6.0

</td>

<td style="text-align:right;">

6

</td>

</tr>

<tr>

<td style="text-align:left;">

Caulobacterales

</td>

<td style="text-align:right;">

23.7

</td>

<td style="text-align:right;">

71

</td>

</tr>

<tr>

<td style="text-align:left;">

Chromatiales

</td>

<td style="text-align:right;">

20.0

</td>

<td style="text-align:right;">

20

</td>

</tr>

<tr>

<td style="text-align:left;">

Chthoniobacterales

</td>

<td style="text-align:right;">

23.0

</td>

<td style="text-align:right;">

23

</td>

</tr>

<tr>

<td style="text-align:left;">

Clostridiales

</td>

<td style="text-align:right;">

6.5

</td>

<td style="text-align:right;">

13

</td>

</tr>

<tr>

<td style="text-align:left;">

Corynebacteriales

</td>

<td style="text-align:right;">

8.0

</td>

<td style="text-align:right;">

40

</td>

</tr>

<tr>

<td style="text-align:left;">

Cytophagales

</td>

<td style="text-align:right;">

11.8

</td>

<td style="text-align:right;">

695

</td>

</tr>

<tr>

<td style="text-align:left;">

Deinococcales

</td>

<td style="text-align:right;">

15.1

</td>

<td style="text-align:right;">

151

</td>

</tr>

<tr>

<td style="text-align:left;">

Desulfovibrionales

</td>

<td style="text-align:right;">

5.0

</td>

<td style="text-align:right;">

5

</td>

</tr>

<tr>

<td style="text-align:left;">

E6aD10

</td>

<td style="text-align:right;">

11.0

</td>

<td style="text-align:right;">

11

</td>

</tr>

<tr>

<td style="text-align:left;">

EMP-G18

</td>

<td style="text-align:right;">

3.0

</td>

<td style="text-align:right;">

3

</td>

</tr>

<tr>

<td style="text-align:left;">

Enterobacteriales

</td>

<td style="text-align:right;">

22.0

</td>

<td style="text-align:right;">

44

</td>

</tr>

<tr>

<td style="text-align:left;">

Euzebyales

</td>

<td style="text-align:right;">

13.5

</td>

<td style="text-align:right;">

202

</td>

</tr>

<tr>

<td style="text-align:left;">

Flavobacteriales

</td>

<td style="text-align:right;">

12.7

</td>

<td style="text-align:right;">

38

</td>

</tr>

<tr>

<td style="text-align:left;">

Frankiales

</td>

<td style="text-align:right;">

20.9

</td>

<td style="text-align:right;">

356

</td>

</tr>

<tr>

<td style="text-align:left;">

Gaiellales

</td>

<td style="text-align:right;">

10.4

</td>

<td style="text-align:right;">

94

</td>

</tr>

<tr>

<td style="text-align:left;">

Gemmatimonadales

</td>

<td style="text-align:right;">

14.3

</td>

<td style="text-align:right;">

214

</td>

</tr>

<tr>

<td style="text-align:left;">

Halobacteriales

</td>

<td style="text-align:right;">

3.9

</td>

<td style="text-align:right;">

90

</td>

</tr>

<tr>

<td style="text-align:left;">

JG30-KF-CM45

</td>

<td style="text-align:right;">

13.0

</td>

<td style="text-align:right;">

468

</td>

</tr>

<tr>

<td style="text-align:left;">

Kineosporiales

</td>

<td style="text-align:right;">

20.2

</td>

<td style="text-align:right;">

121

</td>

</tr>

<tr>

<td style="text-align:left;">

Lactobacillales

</td>

<td style="text-align:right;">

18.5

</td>

<td style="text-align:right;">

37

</td>

</tr>

<tr>

<td style="text-align:left;">

Micrococcales

</td>

<td style="text-align:right;">

18.0

</td>

<td style="text-align:right;">

198

</td>

</tr>

<tr>

<td style="text-align:left;">

Micromonosporales

</td>

<td style="text-align:right;">

19.5

</td>

<td style="text-align:right;">

39

</td>

</tr>

<tr>

<td style="text-align:left;">

Myxococcales

</td>

<td style="text-align:right;">

10.2

</td>

<td style="text-align:right;">

51

</td>

</tr>

<tr>

<td style="text-align:left;">

Nitriliruptorales

</td>

<td style="text-align:right;">

9.7

</td>

<td style="text-align:right;">

58

</td>

</tr>

<tr>

<td style="text-align:left;">

Nitrosomonadales

</td>

<td style="text-align:right;">

11.5

</td>

<td style="text-align:right;">

23

</td>

</tr>

<tr>

<td style="text-align:left;">

Nitrospirales

</td>

<td style="text-align:right;">

6.0

</td>

<td style="text-align:right;">

12

</td>

</tr>

<tr>

<td style="text-align:left;">

Oceanospirillales

</td>

<td style="text-align:right;">

16.0

</td>

<td style="text-align:right;">

16

</td>

</tr>

<tr>

<td style="text-align:left;">

Orbales

</td>

<td style="text-align:right;">

7.0

</td>

<td style="text-align:right;">

7

</td>

</tr>

<tr>

<td style="text-align:left;">

Order\_II

</td>

<td style="text-align:right;">

4.8

</td>

<td style="text-align:right;">

24

</td>

</tr>

<tr>

<td style="text-align:left;">

Order\_III

</td>

<td style="text-align:right;">

8.7

</td>

<td style="text-align:right;">

26

</td>

</tr>

<tr>

<td style="text-align:left;">

PAUC43f\_marine\_benthic\_group

</td>

<td style="text-align:right;">

2.0

</td>

<td style="text-align:right;">

2

</td>

</tr>

<tr>

<td style="text-align:left;">

Planctomycetales

</td>

<td style="text-align:right;">

9.2

</td>

<td style="text-align:right;">

74

</td>

</tr>

<tr>

<td style="text-align:left;">

Propionibacteriales

</td>

<td style="text-align:right;">

14.5

</td>

<td style="text-align:right;">

318

</td>

</tr>

<tr>

<td style="text-align:left;">

Pseudomonadales

</td>

<td style="text-align:right;">

14.7

</td>

<td style="text-align:right;">

103

</td>

</tr>

<tr>

<td style="text-align:left;">

Pseudonocardiales

</td>

<td style="text-align:right;">

18.0

</td>

<td style="text-align:right;">

306

</td>

</tr>

<tr>

<td style="text-align:left;">

Rhizobiales

</td>

<td style="text-align:right;">

16.1

</td>

<td style="text-align:right;">

386

</td>

</tr>

<tr>

<td style="text-align:left;">

Rhodobacterales

</td>

<td style="text-align:right;">

14.8

</td>

<td style="text-align:right;">

163

</td>

</tr>

<tr>

<td style="text-align:left;">

Rhodospirillales

</td>

<td style="text-align:right;">

11.8

</td>

<td style="text-align:right;">

271

</td>

</tr>

<tr>

<td style="text-align:left;">

Rickettsiales

</td>

<td style="text-align:right;">

11.2

</td>

<td style="text-align:right;">

146

</td>

</tr>

<tr>

<td style="text-align:left;">

Rubrobacterales

</td>

<td style="text-align:right;">

22.3

</td>

<td style="text-align:right;">

491

</td>

</tr>

<tr>

<td style="text-align:left;">

S0134\_terrestrial\_group

</td>

<td style="text-align:right;">

18.0

</td>

<td style="text-align:right;">

18

</td>

</tr>

<tr>

<td style="text-align:left;">

SBYG-4553

</td>

<td style="text-align:right;">

10.0

</td>

<td style="text-align:right;">

10

</td>

</tr>

<tr>

<td style="text-align:left;">

Solirubrobacterales

</td>

<td style="text-align:right;">

17.2

</td>

<td style="text-align:right;">

638

</td>

</tr>

<tr>

<td style="text-align:left;">

Sphaerobacterales

</td>

<td style="text-align:right;">

3.3

</td>

<td style="text-align:right;">

10

</td>

</tr>

<tr>

<td style="text-align:left;">

Sphingobacteriales

</td>

<td style="text-align:right;">

14.2

</td>

<td style="text-align:right;">

227

</td>

</tr>

<tr>

<td style="text-align:left;">

Sphingomonadales

</td>

<td style="text-align:right;">

19.6

</td>

<td style="text-align:right;">

196

</td>

</tr>

<tr>

<td style="text-align:left;">

Streptomycetales

</td>

<td style="text-align:right;">

24.5

</td>

<td style="text-align:right;">

49

</td>

</tr>

<tr>

<td style="text-align:left;">

Streptosporangiales

</td>

<td style="text-align:right;">

3.0

</td>

<td style="text-align:right;">

3

</td>

</tr>

<tr>

<td style="text-align:left;">

Subgroup\_3

</td>

<td style="text-align:right;">

15.0

</td>

<td style="text-align:right;">

15

</td>

</tr>

<tr>

<td style="text-align:left;">

Subgroup\_4

</td>

<td style="text-align:right;">

10.6

</td>

<td style="text-align:right;">

85

</td>

</tr>

<tr>

<td style="text-align:left;">

Subgroup\_6

</td>

<td style="text-align:right;">

13.0

</td>

<td style="text-align:right;">

52

</td>

</tr>

<tr>

<td style="text-align:left;">

SubsectionI

</td>

<td style="text-align:right;">

13.0

</td>

<td style="text-align:right;">

13

</td>

</tr>

<tr>

<td style="text-align:left;">

SubsectionII

</td>

<td style="text-align:right;">

15.8

</td>

<td style="text-align:right;">

315

</td>

</tr>

<tr>

<td style="text-align:left;">

SubsectionIII

</td>

<td style="text-align:right;">

14.8

</td>

<td style="text-align:right;">

118

</td>

</tr>

<tr>

<td style="text-align:left;">

SubsectionIV

</td>

<td style="text-align:right;">

13.5

</td>

<td style="text-align:right;">

54

</td>

</tr>

<tr>

<td style="text-align:left;">

TRA3-20

</td>

<td style="text-align:right;">

14.0

</td>

<td style="text-align:right;">

14

</td>

</tr>

<tr>

<td style="text-align:left;">

Thermales

</td>

<td style="text-align:right;">

23.0

</td>

<td style="text-align:right;">

23

</td>

</tr>

<tr>

<td style="text-align:left;">

Thermophilales

</td>

<td style="text-align:right;">

15.0

</td>

<td style="text-align:right;">

30

</td>

</tr>

<tr>

<td style="text-align:left;">

Unclassified

</td>

<td style="text-align:right;">

9.2

</td>

<td style="text-align:right;">

642

</td>

</tr>

<tr>

<td style="text-align:left;">

Unknown\_Order

</td>

<td style="text-align:right;">

9.8

</td>

<td style="text-align:right;">

98

</td>

</tr>

<tr>

<td style="text-align:left;">

Vampirovibrionales

</td>

<td style="text-align:right;">

6.0

</td>

<td style="text-align:right;">

6

</td>

</tr>

<tr>

<td style="text-align:left;">

WD2101\_soil\_group

</td>

<td style="text-align:right;">

13.0

</td>

<td style="text-align:right;">

39

</td>

</tr>

<tr>

<td style="text-align:left;">

Xanthomonadales

</td>

<td style="text-align:right;">

18.0

</td>

<td style="text-align:right;">

18

</td>

</tr>

</tbody>

</table>

Based on that we will remove all phyla with a prevalence of under 8

``` r
Prevalence_phylum_summary %>% 
  filter(`Sum prevalence` < 8) %>% 
  dplyr::select(Phylum) %>% 
  map(as.character) %>% 
  unlist() ->
  filterPhyla
Ps_obj_filt %<>% subset_taxa(!Phylum %in% filterPhyla)
sample_data(Ps_obj_filt)$Library.size <- rowSums(otu_table(Ps_obj_filt))
print(Ps_obj)
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 767 taxa and 25 samples ]
    ## sample_data() Sample Data:       [ 25 samples by 6 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 767 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 767 tips and 765 internal nodes ]
    ## taxa are columns

``` r
print(Ps_obj_filt)
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 732 taxa and 25 samples ]
    ## sample_data() Sample Data:       [ 25 samples by 6 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 732 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 732 tips and 730 internal nodes ]
    ## taxa are columns

#### Plot general prevalence features of the phyla

``` r
# Subset to the remaining phyla
prevdf_phylum_filt <- subset(prevdf, Phylum %in% get_taxa_unique(Ps_obj_filt, "Phylum"))
ggplot(prevdf_phylum_filt,
       aes(TotalAbundance, Prevalence / nsamples(Ps_obj_filt), color = Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05,
             alpha = 0.5,
             linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap( ~ Phylum) + theme(legend.position = "none")
```

![](BRC_growth_rate_figures/prevalence%20phylum-1.svg)<!-- -->

#### Plot general prevalence features of the top 20 orders

``` r
# Subset to the remaining phyla
prevdf_order_filt <- subset(prevdf, Order %in% get_taxa_unique(Ps_obj_filt, "Order"))
# grab the top 30 most abundant orders
prevdf_order_filt %>% 
  group_by(Order) %>%
  summarise(Combined.abundance = sum(TotalAbundance)) %>% 
  arrange(desc(Combined.abundance)) %>% 
  .[1:30, "Order"]  ->
  Orders2plot
prevdf_order_filt2 <- subset(prevdf, Order %in% Orders2plot$Order)
ggplot(prevdf_order_filt2,
       aes(TotalAbundance, Prevalence / nsamples(Ps_obj_filt), color = Order)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05,
             alpha = 0.5,
             linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap( ~ Order) + theme(legend.position = "none")
```

![](BRC_growth_rate_figures/prevalence%20order-1.svg)<!-- -->

#### Unsupervised filtering by prevalence

We will remove all sequences which appear in less than 10% of the
samples

``` r
# Define prevalence threshold as 10% of total samples
prevalenceThreshold <- 0.1 * nsamples(Ps_obj_filt)
prevalenceThreshold
```

    ## [1] 2.5

``` r
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa <-
  row.names(prevdf_phylum_filt)[(prevdf_phylum_filt$Prevalence >= prevalenceThreshold)]
Ps_obj_filt  %<>%  prune_taxa(keepTaxa, .)
sample_data(Ps_obj_filt)$Library.size <- rowSums(otu_table(Ps_obj_filt))
print(Ps_obj)
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 767 taxa and 25 samples ]
    ## sample_data() Sample Data:       [ 25 samples by 6 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 767 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 767 tips and 765 internal nodes ]
    ## taxa are columns

``` r
print(Ps_obj_filt)
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 711 taxa and 25 samples ]
    ## sample_data() Sample Data:       [ 25 samples by 6 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 711 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 711 tips and 709 internal nodes ]
    ## taxa are columns

This removed 56 or 7% of the sequences.

### Exploring the dataset features

First let’s look at the count data distribution after filtering:

``` r
PlotLibDist(Ps_obj_filt)
```

![](BRC_growth_rate_figures/plot%20abundance-1.svg)<!-- -->

``` r
sample_data(Ps_obj_filt) %>% 
  as_tibble() %>% 
  dplyr::select(Sample.name, Library.size) %>% 
  as(., "data.frame") %>% 
  kable(.) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Sample.name

</th>

<th style="text-align:right;">

Library.size

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

City Chalk 2

</td>

<td style="text-align:right;">

49453

</td>

</tr>

<tr>

<td style="text-align:left;">

City Chalk 3

</td>

<td style="text-align:right;">

54157

</td>

</tr>

<tr>

<td style="text-align:left;">

City Chalk 4

</td>

<td style="text-align:right;">

44943

</td>

</tr>

<tr>

<td style="text-align:left;">

City Chalk 5

</td>

<td style="text-align:right;">

58835

</td>

</tr>

<tr>

<td style="text-align:left;">

City Chalk 6

</td>

<td style="text-align:right;">

77463

</td>

</tr>

<tr>

<td style="text-align:left;">

City Limestone 1

</td>

<td style="text-align:right;">

36122

</td>

</tr>

<tr>

<td style="text-align:left;">

City Limestone 2

</td>

<td style="text-align:right;">

52816

</td>

</tr>

<tr>

<td style="text-align:left;">

City Limestone 3

</td>

<td style="text-align:right;">

56370

</td>

</tr>

<tr>

<td style="text-align:left;">

City Limestone 4

</td>

<td style="text-align:right;">

53092

</td>

</tr>

<tr>

<td style="text-align:left;">

City Limestone 5

</td>

<td style="text-align:right;">

61350

</td>

</tr>

<tr>

<td style="text-align:left;">

City Limestone 6

</td>

<td style="text-align:right;">

54508

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope Limestone 1

</td>

<td style="text-align:right;">

75580

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope Chalk 1

</td>

<td style="text-align:right;">

51652

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope Chalk 2

</td>

<td style="text-align:right;">

45201

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope Chalk 3

</td>

<td style="text-align:right;">

33632

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope Chalk 4

</td>

<td style="text-align:right;">

70890

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope Chalk 5

</td>

<td style="text-align:right;">

71460

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope Chalk 6

</td>

<td style="text-align:right;">

58961

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope Chalk 7

</td>

<td style="text-align:right;">

36759

</td>

</tr>

<tr>

<td style="text-align:left;">

City Chalk 1

</td>

<td style="text-align:right;">

23968

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope Limestone 2

</td>

<td style="text-align:right;">

62460

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope Limestone 3

</td>

<td style="text-align:right;">

36635

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope Limestone 4

</td>

<td style="text-align:right;">

57291

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope Limestone 5

</td>

<td style="text-align:right;">

64511

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope Limestone 6

</td>

<td style="text-align:right;">

41565

</td>

</tr>

</tbody>

</table>

The figure and table indicate only a small deviation in the number of
reads per samples.

``` r
(mod1 <- adonis(
  otu_table(Ps_obj_filt) ~ Library.size,
  data = as(sample_data(Ps_obj_filt), "data.frame"), 
  method = "horn",
  permutations = 9999
))
```

    ## 
    ## Call:
    ## adonis(formula = otu_table(Ps_obj_filt) ~ Library.size, data = as(sample_data(Ps_obj_filt),      "data.frame"), permutations = 9999, method = "horn") 
    ## 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##              Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
    ## Library.size  1    0.4358 0.43579  1.3712 0.05626 0.2005
    ## Residuals    23    7.3095 0.31781         0.94374       
    ## Total        24    7.7453                 1.00000

``` r
PlotReadHist(as(otu_table(Ps_obj_filt), "matrix"))
```

![](BRC_growth_rate_figures/mod%20abundance-1.svg)<!-- -->

``` r
notAllZero <- (rowSums(t(otu_table(Ps_obj_filt))) > 0)
vsn::meanSdPlot(as.matrix(log2(t(otu_table(Ps_obj_filt))[notAllZero, ] + 1)))
```

![](BRC_growth_rate_figures/mod%20abundance-2.svg)<!-- -->

The difference in library sizes is low and its effect on the community
composition is minimal. We’ll use the GMPR method for library size
normalisation (Chen and Chen [2017](#ref-chen_gmpr:_2017))

``` r
Ps_obj_filt_GMPR <- Ps_obj_filt
Ps_obj_filt %>%
  otu_table(.) %>%
  t() %>%
  as(., "matrix") %>%
  GMPR() ->
  GMPR_factors
```

    ## Begin GMPR size factor calculation ...
    ## Completed!
    ## Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers!

``` r
Ps_obj_filt %>%
  otu_table(.) %>%
  t() %*% diag(1 / GMPR_factors$gmpr) %>%
  t() %>%
  as.data.frame(., row.names = sample_names(Ps_obj_filt)) %>%
  otu_table(., taxa_are_rows = FALSE) ->
  otu_table(Ps_obj_filt_GMPR)
sample_data(Ps_obj_filt_GMPR)$Library.size <- sample_sums(Ps_obj_filt)
adonis(
  otu_table(Ps_obj_filt_GMPR) ~ Library.size,
  data = as(sample_data(Ps_obj_filt_GMPR), "data.frame"),
  method = "horn",
  permutations = 9999
)
```

    ## 
    ## Call:
    ## adonis(formula = otu_table(Ps_obj_filt_GMPR) ~ Library.size,      data = as(sample_data(Ps_obj_filt_GMPR), "data.frame"), permutations = 9999,      method = "horn") 
    ## 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##              Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
    ## Library.size  1    0.4358 0.43579  1.3712 0.05626 0.2068
    ## Residuals    23    7.3095 0.31781         0.94374       
    ## Total        24    7.7453                 1.00000

``` r
PlotLibDist(Ps_obj_filt_GMPR)
```

![](BRC_growth_rate_figures/GMPR-1.svg)<!-- -->

Did it improve anything?

``` r
PlotReadHist(as(otu_table(Ps_obj_filt_GMPR), "matrix"))
```

![](BRC_growth_rate_figures/GMPR%20diag%20plots-1.svg)<!-- -->

``` r
notAllZero <- (rowSums(t(otu_table(Ps_obj_filt_GMPR))) > 0)
vsn::meanSdPlot(as.matrix(log2(t(otu_table(Ps_obj_filt_GMPR))[notAllZero, ] + 1)))
```

![](BRC_growth_rate_figures/GMPR%20diag%20plots-2.svg)<!-- -->

### Alpha diversity

#### Calculate and plot alpha diversity metrics

We do that by simulating 1000 rarefaction events and calculating the
metrics each time. Then, the result is averaged.

``` r
rarefaction.mat <- matrix(0, nrow = nsamples(Ps_obj_filt), ncol = bootstraps)
rownames(rarefaction.mat) <- sample_names(Ps_obj_filt)
rich.ests <- list(S.obs = rarefaction.mat, S.chao1 = rarefaction.mat, se.chao1 = rarefaction.mat,
                   S.ACE = rarefaction.mat, se.ACE = rarefaction.mat)
for (i in seq(bootstraps)) {
  sub.OTUmat <- rrarefy(otu_table(Ps_obj_filt), min(rowSums(otu_table(Ps_obj_filt))))
  for (j in seq(length(rich.ests))) {
    rich.ests[[j]][, i] <- t(estimateR(sub.OTUmat))[, j]
  }
}
Richness <- data.frame(row.names = row.names(rich.ests[[1]]))
for (i in c(1, seq(2, length(rich.ests), 2))) {
  S <- apply(rich.ests[[i]], 1, mean)
  if (i == 1) { 
    se <- apply(rich.ests[[i]], 1, function(x) (mean(x)/sqrt(length(x))))
    } else se <- apply(rich.ests[[i + 1]], 1, mean)
  Richness <- cbind(Richness, S, se)
}
colnames(Richness) <- c("S.obs", "S.obs.se", "S.chao1", "S.chao1.se", "S.ACE", "S.ACE.se")
saveRDS(Richness, file = paste0("./Results/", Proj_name, "_richness.RDS"))
write.csv(Richness, file = paste0("./Results/", Proj_name, "_richness.csv"))
ses <- grep("\\.se", colnames(Richness))
Richness[, ses] %>% 
  gather(key = "est.se") -> se.dat
Richness[, -unique(ses)] %>% 
  gather(key = "est") -> mean.dat
n <- length(unique(mean.dat$est))
# diversity indices
diversity.inds <- list(Shannon = rarefaction.mat, inv.simpson = rarefaction.mat, BP = rarefaction.mat)
for (i in seq(bootstraps)) {
  sub.OTUmat <- rrarefy(otu_table(Ps_obj_filt), min(rowSums(otu_table(Ps_obj_filt))))
  diversity.inds$Shannon[, i] <- diversityresult(sub.OTUmat, index = 'Shannon', method = 'each site', digits = 3)[, 1]
  diversity.inds$inv.simpson[, i] <- diversityresult(sub.OTUmat, index = 'inverseSimpson', method = 'each site', digits = 3)[, 1]
  diversity.inds$BP[, i] <- diversityresult(sub.OTUmat, index = 'Berger', method = 'each site', digits = 3)[, 1]
}
Diversity <- data.frame(row.names = row.names(diversity.inds[[1]]))
for (i in seq(length(diversity.inds))) {
  S <- apply(diversity.inds[[i]], 1, mean)
  se <- apply(diversity.inds[[i]], 1, function(x) (mean(x)/sqrt(length(x))))
  Diversity <- cbind(Diversity, S, se)
}
colnames(Diversity) <- c("Shannon", "Shannon.se", "Inv.simpson", "Inv.simpson.se", "BP", "BP.se")
ses <- grep("\\.se", colnames(Diversity))
Diversity[, ses] %>% gather(key = "est.se") -> se.dat
Diversity[, -unique(ses)] %>% gather(key = "est") -> mean.dat
saveRDS(Diversity, file = paste0("./Results/", Proj_name, "_diversity.RDS"))
write.csv(Diversity, file = paste0("./Results/", Proj_name, "_diversity.csv"))
```

#### Test the differences in alpha diversity.

``` r
# make combined richness diversity
Richness_Diversity <- cbind(Richness, Diversity)
ses <- grep("\\.se", colnames(Richness_Diversity))
Richness_Diversity[, ses] %>% 
  gather(key = "est.se") -> 
  se.dat
Richness_Diversity[, -unique(ses)] %>% 
  gather(key = "Metric", 
         value = "Estimate") -> 
  mean.dat
Richness_Diversity_long <-
  cbind(
    Sample = rep(rownames(Richness_Diversity), times = length(unique(mean.dat$Metric))),
    mean.dat,
    lerr = mean.dat$Estimate - se.dat$value,
    herr = mean.dat$Estimate + se.dat$value
  )
Richness_Diversity_long$Metric <-
  factor(
    Richness_Diversity_long$Metric,
    levels = c("S.obs", "S.chao1", "S.ACE", "Shannon", "Inv.simpson", "BP"),
    labels = c("S obs.", "Chao1", "ACE", "Shannon", "Inv. Simpson" , "Berger Parker")
  )
Richness_Diversity_long %<>%
  cbind(., 
        sample_data(Ps_obj_filt))

# S Obs
data2test <- Richness_Diversity_long[Richness_Diversity_long$Metric == "S obs.", ] 

(mod_obsS <- TestAlphaV3(filter(Richness_Diversity_long, Metric == "S obs.")))
```

    ## Call:
    ##    aov(formula = as.formula(paste(response_name, paste(factor_names[1], 
    ##     factor_names[2], sep = " * "), sep = " ~ ")), data = data2test)
    ## 
    ## Terms:
    ##                 Location Rock.type Location:Rock.type Residuals
    ## Sum of Squares        65       705               9787     62082
    ## Deg. of Freedom        1         1                  1        21
    ## 
    ## Residual standard error: 54.4
    ## Estimated effects may be unbalanced

![](BRC_growth_rate_figures/test%20alpha-1.svg)<!-- -->

    ## [1] "Unequal group sizes - showing SS type III"
    ## Anova Table (Type III tests)
    ## 
    ## Response: Estimate
    ##                    Sum Sq Df F value  Pr(>F)    
    ## (Intercept)        594952  1  201.25 3.1e-12 ***
    ## Location             3860  1    1.31   0.266    
    ## Rock.type            7702  1    2.61   0.121    
    ## Location:Rock.type   9787  1    3.31   0.083 .  
    ## Residuals           62082 21                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Tables of means
    ## Grand mean
    ##     
    ## 313 
    ## 
    ##  Location 
    ##     Slope City
    ##       314  311
    ## rep    13   12
    ## 
    ##  Rock.type 
    ##     Chalk Limestone
    ##       307       318
    ## rep    13        12
    ## 
    ##  Location:Rock.type 
    ##         Rock.type
    ## Location Chalk Limestone
    ##    Slope 292   340      
    ##    rep     7     6      
    ##    City  326   296      
    ##    rep     6     6

    ## Call:
    ##    aov(formula = as.formula(paste(response_name, paste(factor_names[1], 
    ##     factor_names[2], sep = " * "), sep = " ~ ")), data = data2test)
    ## 
    ## Terms:
    ##                 Location Rock.type Location:Rock.type Residuals
    ## Sum of Squares        65       705               9787     62082
    ## Deg. of Freedom        1         1                  1        21
    ## 
    ## Residual standard error: 54.4
    ## Estimated effects may be unbalanced

``` r
# Post-hoc test
marginal <- emmeans(mod_obsS,
                   ~ Location : Rock.type)
summary(marginal)
```

<div class="kable-table">

<table>

<thead>

<tr>

<th style="text-align:left;">

Location

</th>

<th style="text-align:left;">

Rock.type

</th>

<th style="text-align:right;">

emmean

</th>

<th style="text-align:right;">

SE

</th>

<th style="text-align:right;">

df

</th>

<th style="text-align:right;">

lower.CL

</th>

<th style="text-align:right;">

upper.CL

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

292

</td>

<td style="text-align:right;">

20.6

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

249

</td>

<td style="text-align:right;">

334

</td>

</tr>

<tr>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

326

</td>

<td style="text-align:right;">

22.2

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

280

</td>

<td style="text-align:right;">

372

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

340

</td>

<td style="text-align:right;">

22.2

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

294

</td>

<td style="text-align:right;">

387

</td>

</tr>

<tr>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

296

</td>

<td style="text-align:right;">

22.2

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

249

</td>

<td style="text-align:right;">

342

</td>

</tr>

</tbody>

</table>

</div>

``` r
contrast(marginal, 
         method = "pairwise", 
         adjust = "tukey")
```

    ##  contrast                         estimate   SE df t.ratio p.value
    ##  Slope Chalk - City Chalk            -34.6 30.2 21 -1.143  0.6680 
    ##  Slope Chalk - Slope Limestone       -48.8 30.2 21 -1.614  0.3930 
    ##  Slope Chalk - City Limestone         -4.1 30.2 21 -0.135  0.9990 
    ##  City Chalk - Slope Limestone        -14.3 31.4 21 -0.454  0.9680 
    ##  City Chalk - City Limestone          30.5 31.4 21  0.971  0.7670 
    ##  Slope Limestone - City Limestone     44.8 31.4 21  1.426  0.4980 
    ## 
    ## P value adjustment: tukey method for comparing a family of 4 estimates

``` r
(obsS_pairwise <- cld(marginal,
                      alpha = 0.05,
                      Letters = letters,
                      adjust = "tukey")) # works with lm but not with two-factor ART
```

<div class="kable-table">

<table>

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:left;">

Location

</th>

<th style="text-align:left;">

Rock.type

</th>

<th style="text-align:right;">

emmean

</th>

<th style="text-align:right;">

SE

</th>

<th style="text-align:right;">

df

</th>

<th style="text-align:right;">

lower.CL

</th>

<th style="text-align:right;">

upper.CL

</th>

<th style="text-align:left;">

.group

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

1

</td>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

292

</td>

<td style="text-align:right;">

20.6

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

236

</td>

<td style="text-align:right;">

347

</td>

<td style="text-align:left;">

a

</td>

</tr>

<tr>

<td style="text-align:left;">

4

</td>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

296

</td>

<td style="text-align:right;">

22.2

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

235

</td>

<td style="text-align:right;">

356

</td>

<td style="text-align:left;">

a

</td>

</tr>

<tr>

<td style="text-align:left;">

2

</td>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

326

</td>

<td style="text-align:right;">

22.2

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

266

</td>

<td style="text-align:right;">

387

</td>

<td style="text-align:left;">

a

</td>

</tr>

<tr>

<td style="text-align:left;">

3

</td>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

340

</td>

<td style="text-align:right;">

22.2

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

280

</td>

<td style="text-align:right;">

401

</td>

<td style="text-align:left;">

a

</td>

</tr>

</tbody>

</table>

</div>

``` r
(mod_obsS %>% 
  anova() %>% 
  mutate(`Part Eta Sq`=`Sum Sq`/sum(`Sum Sq`) ) ->
  mod_obsS_ANOVA)
```

<div class="kable-table">

<table>

<thead>

<tr>

<th style="text-align:right;">

Df

</th>

<th style="text-align:right;">

Sum Sq

</th>

<th style="text-align:right;">

Mean Sq

</th>

<th style="text-align:right;">

F value

</th>

<th style="text-align:right;">

Pr(\>F)

</th>

<th style="text-align:right;">

Part Eta Sq

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

64.6

</td>

<td style="text-align:right;">

64.6

</td>

<td style="text-align:right;">

0.022

</td>

<td style="text-align:right;">

0.884

</td>

<td style="text-align:right;">

0.001

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

704.8

</td>

<td style="text-align:right;">

704.8

</td>

<td style="text-align:right;">

0.238

</td>

<td style="text-align:right;">

0.630

</td>

<td style="text-align:right;">

0.010

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

9786.6

</td>

<td style="text-align:right;">

9786.6

</td>

<td style="text-align:right;">

3.310

</td>

<td style="text-align:right;">

0.083

</td>

<td style="text-align:right;">

0.135

</td>

</tr>

<tr>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

62081.9

</td>

<td style="text-align:right;">

2956.3

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

0.855

</td>

</tr>

</tbody>

</table>

</div>

``` r
# pwpp(marginal) # Pairwise P-value plot. Fails for unbalanced design
emmip(mod_obsS, Location ~ Rock.type)
```

![](BRC_growth_rate_figures/test%20alpha-2.svg)<!-- -->

``` r
# summary(as.glht(pairs(marginal))) # fails because of unbalanced design

# Shannon
(mod_Shannon <- TestAlphaV3(filter(Richness_Diversity_long, Metric == "Shannon")))
```

    ## Call:
    ##    aov(formula = as.formula(paste(response_name, paste(factor_names[1], 
    ##     factor_names[2], sep = " * "), sep = " ~ ")), data = data2test)
    ## 
    ## Terms:
    ##                 Location Rock.type Location:Rock.type Residuals
    ## Sum of Squares      1.07      0.18               0.16      4.41
    ## Deg. of Freedom        1         1                  1        21
    ## 
    ## Residual standard error: 0.458
    ## Estimated effects may be unbalanced

![](BRC_growth_rate_figures/test%20alpha-3.svg)<!-- -->

    ## [1] "Unequal group sizes - showing SS type III"
    ## Anova Table (Type III tests)
    ## 
    ## Response: Estimate
    ##                    Sum Sq Df F value Pr(>F)    
    ## (Intercept)         313.2  1 1491.16 <2e-16 ***
    ## Location              1.1  1    5.38   0.03 *  
    ## Rock.type             0.2  1    0.78   0.39    
    ## Location:Rock.type    0.2  1    0.74   0.40    
    ## Residuals             4.4 21                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Tables of means
    ## Grand mean
    ##      
    ## 3.55 
    ## 
    ##  Location 
    ##     Slope  City
    ##      3.75  3.33
    ## rep 13.00 12.00
    ## 
    ##  Rock.type 
    ##     Chalk Limestone
    ##      3.47      3.64
    ## rep 13.00     12.00
    ## 
    ##  Location:Rock.type 
    ##         Rock.type
    ## Location Chalk Limestone
    ##    Slope 3.60  3.92     
    ##    rep   7.00  6.00     
    ##    City  3.33  3.34     
    ##    rep   6.00  6.00

    ## Call:
    ##    aov(formula = as.formula(paste(response_name, paste(factor_names[1], 
    ##     factor_names[2], sep = " * "), sep = " ~ ")), data = data2test)
    ## 
    ## Terms:
    ##                 Location Rock.type Location:Rock.type Residuals
    ## Sum of Squares      1.07      0.18               0.16      4.41
    ## Deg. of Freedom        1         1                  1        21
    ## 
    ## Residual standard error: 0.458
    ## Estimated effects may be unbalanced

``` r
# Post-hoc test
marginal <- emmeans(mod_Shannon,
                   ~ Location : Rock.type)
summary(marginal)
```

<div class="kable-table">

<table>

<thead>

<tr>

<th style="text-align:left;">

Location

</th>

<th style="text-align:left;">

Rock.type

</th>

<th style="text-align:right;">

emmean

</th>

<th style="text-align:right;">

SE

</th>

<th style="text-align:right;">

df

</th>

<th style="text-align:right;">

lower.CL

</th>

<th style="text-align:right;">

upper.CL

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

3.60

</td>

<td style="text-align:right;">

0.173

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

3.24

</td>

<td style="text-align:right;">

3.96

</td>

</tr>

<tr>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

3.33

</td>

<td style="text-align:right;">

0.187

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

2.94

</td>

<td style="text-align:right;">

3.72

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

3.92

</td>

<td style="text-align:right;">

0.187

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

3.53

</td>

<td style="text-align:right;">

4.31

</td>

</tr>

<tr>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

3.34

</td>

<td style="text-align:right;">

0.187

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

2.95

</td>

<td style="text-align:right;">

3.73

</td>

</tr>

</tbody>

</table>

</div>

``` r
contrast(marginal, 
         method = "pairwise", 
         adjust = "tukey")
```

    ##  contrast                         estimate    SE df t.ratio p.value
    ##  Slope Chalk - City Chalk            0.268 0.255 21  1.050  0.7220 
    ##  Slope Chalk - Slope Limestone      -0.321 0.255 21 -1.258  0.5980 
    ##  Slope Chalk - City Limestone        0.264 0.255 21  1.035  0.7310 
    ##  City Chalk - Slope Limestone       -0.589 0.265 21 -2.225  0.1490 
    ##  City Chalk - City Limestone        -0.004 0.265 21 -0.015  1.0000 
    ##  Slope Limestone - City Limestone    0.585 0.265 21  2.210  0.1530 
    ## 
    ## P value adjustment: tukey method for comparing a family of 4 estimates

``` r
(Shannon_pairwise <- cld(marginal,
                      alpha = 0.05,
                      Letters = letters,
                      adjust = "tukey")) # works with lm but not with two-factor ART
```

<div class="kable-table">

<table>

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:left;">

Location

</th>

<th style="text-align:left;">

Rock.type

</th>

<th style="text-align:right;">

emmean

</th>

<th style="text-align:right;">

SE

</th>

<th style="text-align:right;">

df

</th>

<th style="text-align:right;">

lower.CL

</th>

<th style="text-align:right;">

upper.CL

</th>

<th style="text-align:left;">

.group

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

2

</td>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

3.33

</td>

<td style="text-align:right;">

0.187

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

2.82

</td>

<td style="text-align:right;">

3.84

</td>

<td style="text-align:left;">

a

</td>

</tr>

<tr>

<td style="text-align:left;">

4

</td>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

3.34

</td>

<td style="text-align:right;">

0.187

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

2.83

</td>

<td style="text-align:right;">

3.85

</td>

<td style="text-align:left;">

a

</td>

</tr>

<tr>

<td style="text-align:left;">

1

</td>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

3.60

</td>

<td style="text-align:right;">

0.173

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

3.13

</td>

<td style="text-align:right;">

4.07

</td>

<td style="text-align:left;">

a

</td>

</tr>

<tr>

<td style="text-align:left;">

3

</td>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

3.92

</td>

<td style="text-align:right;">

0.187

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

3.41

</td>

<td style="text-align:right;">

4.43

</td>

<td style="text-align:left;">

a

</td>

</tr>

</tbody>

</table>

</div>

``` r
(mod_Shannon %>% 
  anova() %>% 
  mutate(`Part Eta Sq`=`Sum Sq`/sum(`Sum Sq`) ) ->
  mod_Shannon_ANOVA)
```

<div class="kable-table">

<table>

<thead>

<tr>

<th style="text-align:right;">

Df

</th>

<th style="text-align:right;">

Sum Sq

</th>

<th style="text-align:right;">

Mean Sq

</th>

<th style="text-align:right;">

F value

</th>

<th style="text-align:right;">

Pr(\>F)

</th>

<th style="text-align:right;">

Part Eta Sq

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1.069

</td>

<td style="text-align:right;">

1.069

</td>

<td style="text-align:right;">

5.090

</td>

<td style="text-align:right;">

0.035

</td>

<td style="text-align:right;">

0.184

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.176

</td>

<td style="text-align:right;">

0.176

</td>

<td style="text-align:right;">

0.839

</td>

<td style="text-align:right;">

0.370

</td>

<td style="text-align:right;">

0.030

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.156

</td>

<td style="text-align:right;">

0.156

</td>

<td style="text-align:right;">

0.744

</td>

<td style="text-align:right;">

0.398

</td>

<td style="text-align:right;">

0.027

</td>

</tr>

<tr>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

4.411

</td>

<td style="text-align:right;">

0.210

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

0.759

</td>

</tr>

</tbody>

</table>

</div>

``` r
# pwpp(marginal) # Pairwise P-value plot. Fails for unbalanced design
emmip(mod_Shannon, Location ~ Rock.type)
```

![](BRC_growth_rate_figures/test%20alpha-4.svg)<!-- -->

``` r
# ACE
(mod_ACE <- TestAlphaV3(filter(Richness_Diversity_long, Metric == "ACE")))
```

    ## Call:
    ##    aov(formula = as.formula(paste(response_name, paste(factor_names[1], 
    ##     factor_names[2], sep = " * "), sep = " ~ ")), data = data2test)
    ## 
    ## Terms:
    ##                 Location Rock.type Location:Rock.type Residuals
    ## Sum of Squares         3      1924               8494     70754
    ## Deg. of Freedom        1         1                  1        21
    ## 
    ## Residual standard error: 58
    ## Estimated effects may be unbalanced

![](BRC_growth_rate_figures/test%20alpha-5.svg)<!-- -->

    ## [1] "Unequal group sizes - showing SS type III"
    ## Anova Table (Type III tests)
    ## 
    ## Response: Estimate
    ##                     Sum Sq Df F value Pr(>F)    
    ## (Intercept)        4356003  1 1292.88 <2e-16 ***
    ## Location                12  1    0.00   0.95    
    ## Rock.type             1634  1    0.48   0.49    
    ## Location:Rock.type    8494  1    2.52   0.13    
    ## Residuals            70754 21                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Tables of means
    ## Grand mean
    ##     
    ## 417 
    ## 
    ##  Location 
    ##     Slope City
    ##       417  418
    ## rep    13   12
    ## 
    ##  Rock.type 
    ##     Chalk Limestone
    ##       409       426
    ## rep    13        12
    ## 
    ##  Location:Rock.type 
    ##         Rock.type
    ## Location Chalk Limestone
    ##    Slope 392   446      
    ##    rep     7     6      
    ##    City  428   407      
    ##    rep     6     6

    ## Call:
    ##    aov(formula = as.formula(paste(response_name, paste(factor_names[1], 
    ##     factor_names[2], sep = " * "), sep = " ~ ")), data = data2test)
    ## 
    ## Terms:
    ##                 Location Rock.type Location:Rock.type Residuals
    ## Sum of Squares         3      1924               8494     70754
    ## Deg. of Freedom        1         1                  1        21
    ## 
    ## Residual standard error: 58
    ## Estimated effects may be unbalanced

``` r
# Post-hoc test
marginal <- emmeans(mod_ACE,
                   ~ Location : Rock.type)
summary(marginal)
```

<div class="kable-table">

<table>

<thead>

<tr>

<th style="text-align:left;">

Location

</th>

<th style="text-align:left;">

Rock.type

</th>

<th style="text-align:right;">

emmean

</th>

<th style="text-align:right;">

SE

</th>

<th style="text-align:right;">

df

</th>

<th style="text-align:right;">

lower.CL

</th>

<th style="text-align:right;">

upper.CL

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

392

</td>

<td style="text-align:right;">

21.9

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

347

</td>

<td style="text-align:right;">

438

</td>

</tr>

<tr>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

428

</td>

<td style="text-align:right;">

23.7

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

379

</td>

<td style="text-align:right;">

477

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

446

</td>

<td style="text-align:right;">

23.7

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

396

</td>

<td style="text-align:right;">

495

</td>

</tr>

<tr>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

407

</td>

<td style="text-align:right;">

23.7

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

358

</td>

<td style="text-align:right;">

457

</td>

</tr>

</tbody>

</table>

</div>

``` r
contrast(marginal, 
         method = "pairwise", 
         adjust = "tukey")
```

    ##  contrast                         estimate   SE df t.ratio p.value
    ##  Slope Chalk - City Chalk            -35.5 32.3 21 -1.101  0.6930 
    ##  Slope Chalk - Slope Limestone       -53.2 32.3 21 -1.646  0.3760 
    ##  Slope Chalk - City Limestone        -14.8 32.3 21 -0.458  0.9670 
    ##  City Chalk - Slope Limestone        -17.6 33.5 21 -0.525  0.9520 
    ##  City Chalk - City Limestone          20.7 33.5 21  0.619  0.9250 
    ##  Slope Limestone - City Limestone     38.4 33.5 21  1.144  0.6670 
    ## 
    ## P value adjustment: tukey method for comparing a family of 4 estimates

``` r
(ACE_pairwise <- cld(marginal,
                      alpha = 0.05,
                      Letters = letters,
                      adjust = "tukey")) # works with lm but not with two-factor ART
```

<div class="kable-table">

<table>

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:left;">

Location

</th>

<th style="text-align:left;">

Rock.type

</th>

<th style="text-align:right;">

emmean

</th>

<th style="text-align:right;">

SE

</th>

<th style="text-align:right;">

df

</th>

<th style="text-align:right;">

lower.CL

</th>

<th style="text-align:right;">

upper.CL

</th>

<th style="text-align:left;">

.group

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

1

</td>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

392

</td>

<td style="text-align:right;">

21.9

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

333

</td>

<td style="text-align:right;">

452

</td>

<td style="text-align:left;">

a

</td>

</tr>

<tr>

<td style="text-align:left;">

4

</td>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

407

</td>

<td style="text-align:right;">

23.7

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

343

</td>

<td style="text-align:right;">

472

</td>

<td style="text-align:left;">

a

</td>

</tr>

<tr>

<td style="text-align:left;">

2

</td>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

428

</td>

<td style="text-align:right;">

23.7

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

363

</td>

<td style="text-align:right;">

493

</td>

<td style="text-align:left;">

a

</td>

</tr>

<tr>

<td style="text-align:left;">

3

</td>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

446

</td>

<td style="text-align:right;">

23.7

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

381

</td>

<td style="text-align:right;">

510

</td>

<td style="text-align:left;">

a

</td>

</tr>

</tbody>

</table>

</div>

``` r
(mod_ACE %>% 
  anova() %>% 
  mutate(`Part Eta Sq`=`Sum Sq`/sum(`Sum Sq`) ) ->
  mod_ACE_ANOVA)
```

<div class="kable-table">

<table>

<thead>

<tr>

<th style="text-align:right;">

Df

</th>

<th style="text-align:right;">

Sum Sq

</th>

<th style="text-align:right;">

Mean Sq

</th>

<th style="text-align:right;">

F value

</th>

<th style="text-align:right;">

Pr(\>F)

</th>

<th style="text-align:right;">

Part Eta Sq

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

2.54

</td>

<td style="text-align:right;">

2.54

</td>

<td style="text-align:right;">

0.001

</td>

<td style="text-align:right;">

0.978

</td>

<td style="text-align:right;">

0.000

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1923.90

</td>

<td style="text-align:right;">

1923.90

</td>

<td style="text-align:right;">

0.571

</td>

<td style="text-align:right;">

0.458

</td>

<td style="text-align:right;">

0.024

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

8494.42

</td>

<td style="text-align:right;">

8494.42

</td>

<td style="text-align:right;">

2.521

</td>

<td style="text-align:right;">

0.127

</td>

<td style="text-align:right;">

0.105

</td>

</tr>

<tr>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

70753.94

</td>

<td style="text-align:right;">

3369.24

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

0.872

</td>

</tr>

</tbody>

</table>

</div>

``` r
# pwpp(marginal) # Pairwise P-value plot. Fails for unbalanced design
emmip(mod_ACE, Location ~ Rock.type)
```

![](BRC_growth_rate_figures/test%20alpha-6.svg)<!-- -->

``` r
# summary(as.glht(pairs(marginal))) # fails because of unbalanced design

#Inv. Simpson
(mod_InvSim <- TestAlphaV3(filter(Richness_Diversity_long, Metric == "Inv. Simpson")))
```

    ## Call:
    ##    aov(formula = as.formula(paste(response_name, paste(factor_names[1], 
    ##     factor_names[2], sep = " * "), sep = " ~ ")), data = data2test)
    ## 
    ## Terms:
    ##                 Location Rock.type Location:Rock.type Residuals
    ## Sum of Squares       470       125                 99      1025
    ## Deg. of Freedom        1         1                  1        21
    ## 
    ## Residual standard error: 6.99
    ## Estimated effects may be unbalanced

![](BRC_growth_rate_figures/test%20alpha-7.svg)<!-- -->

    ## [1] "Unequal group sizes - showing SS type III"
    ## Anova Table (Type III tests)
    ## 
    ## Response: Estimate
    ##                    Sum Sq Df F value Pr(>F)    
    ## (Intercept)          7498  1  153.60  4e-11 ***
    ## Location              504  1   10.33 0.0042 ** 
    ## Rock.type             117  1    2.40 0.1364    
    ## Location:Rock.type     99  1    2.03 0.1686    
    ## Residuals            1025 21                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Tables of means
    ## Grand mean
    ##      
    ## 17.4 
    ## 
    ##  Location 
    ##     Slope City
    ##      21.5 12.9
    ## rep  13.0 12.0
    ## 
    ##  Rock.type 
    ##     Chalk Limestone
    ##      15.2      19.7
    ## rep  13.0      12.0
    ## 
    ##  Location:Rock.type 
    ##         Rock.type
    ## Location Chalk Limestone
    ##    Slope 17.7  26.0     
    ##    rep    7.0   6.0     
    ##    City  12.7  13.0     
    ##    rep    6.0   6.0

    ## Call:
    ##    aov(formula = as.formula(paste(response_name, paste(factor_names[1], 
    ##     factor_names[2], sep = " * "), sep = " ~ ")), data = data2test)
    ## 
    ## Terms:
    ##                 Location Rock.type Location:Rock.type Residuals
    ## Sum of Squares       470       125                 99      1025
    ## Deg. of Freedom        1         1                  1        21
    ## 
    ## Residual standard error: 6.99
    ## Estimated effects may be unbalanced

``` r
# Post-hoc test
marginal <- emmeans(mod_InvSim,
                   ~ Location : Rock.type)
summary(marginal)
```

<div class="kable-table">

<table>

<thead>

<tr>

<th style="text-align:left;">

Location

</th>

<th style="text-align:left;">

Rock.type

</th>

<th style="text-align:right;">

emmean

</th>

<th style="text-align:right;">

SE

</th>

<th style="text-align:right;">

df

</th>

<th style="text-align:right;">

lower.CL

</th>

<th style="text-align:right;">

upper.CL

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

17.7

</td>

<td style="text-align:right;">

2.64

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

12.20

</td>

<td style="text-align:right;">

23.2

</td>

</tr>

<tr>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

12.7

</td>

<td style="text-align:right;">

2.85

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

6.75

</td>

<td style="text-align:right;">

18.6

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

26.0

</td>

<td style="text-align:right;">

2.85

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

20.09

</td>

<td style="text-align:right;">

32.0

</td>

</tr>

<tr>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

13.0

</td>

<td style="text-align:right;">

2.85

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

7.09

</td>

<td style="text-align:right;">

19.0

</td>

</tr>

</tbody>

</table>

</div>

``` r
contrast(marginal, 
         method = "pairwise", 
         adjust = "tukey")
```

    ##  contrast                         estimate   SE df t.ratio p.value
    ##  Slope Chalk - City Chalk             5.01 3.89 21  1.290  0.5800 
    ##  Slope Chalk - Slope Limestone       -8.33 3.89 21 -2.140  0.1720 
    ##  Slope Chalk - City Limestone         4.67 3.89 21  1.200  0.6330 
    ##  City Chalk - Slope Limestone       -13.34 4.03 21 -3.310  0.0160 
    ##  City Chalk - City Limestone         -0.34 4.03 21 -0.090  1.0000 
    ##  Slope Limestone - City Limestone    13.00 4.03 21  3.220  0.0200 
    ## 
    ## P value adjustment: tukey method for comparing a family of 4 estimates

``` r
(InvSim_pairwise <- cld(marginal,
                      alpha = 0.05,
                      Letters = letters,
                      adjust = "tukey")) # works with lm but not with two-factor ART
```

<div class="kable-table">

<table>

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:left;">

Location

</th>

<th style="text-align:left;">

Rock.type

</th>

<th style="text-align:right;">

emmean

</th>

<th style="text-align:right;">

SE

</th>

<th style="text-align:right;">

df

</th>

<th style="text-align:right;">

lower.CL

</th>

<th style="text-align:right;">

upper.CL

</th>

<th style="text-align:left;">

.group

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

2

</td>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

12.7

</td>

<td style="text-align:right;">

2.85

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

4.92

</td>

<td style="text-align:right;">

20.4

</td>

<td style="text-align:left;">

a

</td>

</tr>

<tr>

<td style="text-align:left;">

4

</td>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

13.0

</td>

<td style="text-align:right;">

2.85

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

5.26

</td>

<td style="text-align:right;">

20.8

</td>

<td style="text-align:left;">

a

</td>

</tr>

<tr>

<td style="text-align:left;">

1

</td>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

17.7

</td>

<td style="text-align:right;">

2.64

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

10.50

</td>

<td style="text-align:right;">

24.9

</td>

<td style="text-align:left;">

ab

</td>

</tr>

<tr>

<td style="text-align:left;">

3

</td>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

26.0

</td>

<td style="text-align:right;">

2.85

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

18.26

</td>

<td style="text-align:right;">

33.8

</td>

<td style="text-align:left;">

b

</td>

</tr>

</tbody>

</table>

</div>

``` r
(mod_InvSim %>% 
  anova() %>% 
  mutate(`Part Eta Sq`=`Sum Sq`/sum(`Sum Sq`) ) ->
  mod_InvSim_ANOVA)
```

<div class="kable-table">

<table>

<thead>

<tr>

<th style="text-align:right;">

Df

</th>

<th style="text-align:right;">

Sum Sq

</th>

<th style="text-align:right;">

Mean Sq

</th>

<th style="text-align:right;">

F value

</th>

<th style="text-align:right;">

Pr(\>F)

</th>

<th style="text-align:right;">

Part Eta Sq

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

470.4

</td>

<td style="text-align:right;">

470.4

</td>

<td style="text-align:right;">

9.64

</td>

<td style="text-align:right;">

0.005

</td>

<td style="text-align:right;">

0.273

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

125.3

</td>

<td style="text-align:right;">

125.3

</td>

<td style="text-align:right;">

2.57

</td>

<td style="text-align:right;">

0.124

</td>

<td style="text-align:right;">

0.073

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

99.3

</td>

<td style="text-align:right;">

99.3

</td>

<td style="text-align:right;">

2.03

</td>

<td style="text-align:right;">

0.169

</td>

<td style="text-align:right;">

0.058

</td>

</tr>

<tr>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

1025.0

</td>

<td style="text-align:right;">

48.8

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

0.596

</td>

</tr>

</tbody>

</table>

</div>

``` r
# pwpp(marginal) # Pairwise P-value plot. Fails for unbalanced design
emmip(mod_InvSim, Location ~ Rock.type)
```

![](BRC_growth_rate_figures/test%20alpha-8.svg)<!-- -->

``` r
# summary(as.glht(pairs(marginal))) # fails because of unbalanced design


#Berger Parker
(mod_BP <- TestAlphaV3(filter(Richness_Diversity_long, Metric == "Berger Parker")))
```

    ## Call:
    ##    aov(formula = as.formula(paste(response_name, paste(factor_names[1], 
    ##     factor_names[2], sep = " * "), sep = " ~ ")), data = data2test)
    ## 
    ## Terms:
    ##                 Location Rock.type Location:Rock.type Residuals
    ## Sum of Squares    0.0643    0.0133             0.0000    0.1080
    ## Deg. of Freedom        1         1                  1        21
    ## 
    ## Residual standard error: 0.0717
    ## Estimated effects may be unbalanced

![](BRC_growth_rate_figures/test%20alpha-9.svg)<!-- -->

    ## [1] "Unequal group sizes - showing SS type III"
    ## Anova Table (Type III tests)
    ## 
    ## Response: Estimate
    ##                    Sum Sq Df F value  Pr(>F)    
    ## (Intercept)         0.924  1  179.69 9.2e-12 ***
    ## Location            0.066  1   12.92  0.0017 ** 
    ## Rock.type           0.013  1    2.59  0.1228    
    ## Location:Rock.type  0.000  1    0.00  0.9939    
    ## Residuals           0.108 21                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Tables of means
    ## Grand mean
    ##       
    ## 0.192 
    ## 
    ##  Location 
    ##      Slope   City
    ##      0.143  0.244
    ## rep 13.000 12.000
    ## 
    ##  Rock.type 
    ##      Chalk Limestone
    ##      0.214     0.168
    ## rep 13.000    12.000
    ## 
    ##  Location:Rock.type 
    ##         Rock.type
    ## Location Chalk Limestone
    ##    Slope 0.16  0.12     
    ##    rep   7.00  6.00     
    ##    City  0.27  0.22     
    ##    rep   6.00  6.00

    ## Call:
    ##    aov(formula = as.formula(paste(response_name, paste(factor_names[1], 
    ##     factor_names[2], sep = " * "), sep = " ~ ")), data = data2test)
    ## 
    ## Terms:
    ##                 Location Rock.type Location:Rock.type Residuals
    ## Sum of Squares    0.0643    0.0133             0.0000    0.1080
    ## Deg. of Freedom        1         1                  1        21
    ## 
    ## Residual standard error: 0.0717
    ## Estimated effects may be unbalanced

``` r
# Post-hoc test
marginal <- emmeans(mod_BP,
                   ~ Location : Rock.type)
summary(marginal)
```

<div class="kable-table">

<table>

<thead>

<tr>

<th style="text-align:left;">

Location

</th>

<th style="text-align:left;">

Rock.type

</th>

<th style="text-align:right;">

emmean

</th>

<th style="text-align:right;">

SE

</th>

<th style="text-align:right;">

df

</th>

<th style="text-align:right;">

lower.CL

</th>

<th style="text-align:right;">

upper.CL

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

0.164

</td>

<td style="text-align:right;">

0.027

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

0.108

</td>

<td style="text-align:right;">

0.220

</td>

</tr>

<tr>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

0.268

</td>

<td style="text-align:right;">

0.029

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

0.207

</td>

<td style="text-align:right;">

0.328

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

0.118

</td>

<td style="text-align:right;">

0.029

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

0.057

</td>

<td style="text-align:right;">

0.179

</td>

</tr>

<tr>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

0.221

</td>

<td style="text-align:right;">

0.029

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

0.160

</td>

<td style="text-align:right;">

0.282

</td>

</tr>

</tbody>

</table>

</div>

``` r
contrast(marginal, 
         method = "pairwise", 
         adjust = "tukey")
```

    ##  contrast                         estimate     SE df t.ratio p.value
    ##  Slope Chalk - City Chalk          -0.1035 0.0399 21 -2.600  0.0740 
    ##  Slope Chalk - Slope Limestone      0.0460 0.0399 21  1.150  0.6620 
    ##  Slope Chalk - City Limestone      -0.0571 0.0399 21 -1.430  0.4950 
    ##  City Chalk - Slope Limestone       0.1495 0.0414 21  3.610  0.0080 
    ##  City Chalk - City Limestone        0.0464 0.0414 21  1.120  0.6800 
    ##  Slope Limestone - City Limestone  -0.1031 0.0414 21 -2.490  0.0910 
    ## 
    ## P value adjustment: tukey method for comparing a family of 4 estimates

``` r
(BP_pairwise <- cld(marginal,
                      alpha = 0.05,
                      Letters = letters,
                      adjust = "tukey")) # works with lm but not with two-factor ART
```

<div class="kable-table">

<table>

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:left;">

Location

</th>

<th style="text-align:left;">

Rock.type

</th>

<th style="text-align:right;">

emmean

</th>

<th style="text-align:right;">

SE

</th>

<th style="text-align:right;">

df

</th>

<th style="text-align:right;">

lower.CL

</th>

<th style="text-align:right;">

upper.CL

</th>

<th style="text-align:left;">

.group

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

3

</td>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

0.118

</td>

<td style="text-align:right;">

0.029

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

0.038

</td>

<td style="text-align:right;">

0.198

</td>

<td style="text-align:left;">

a

</td>

</tr>

<tr>

<td style="text-align:left;">

1

</td>

<td style="text-align:left;">

Slope

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

0.164

</td>

<td style="text-align:right;">

0.027

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

0.090

</td>

<td style="text-align:right;">

0.238

</td>

<td style="text-align:left;">

ab

</td>

</tr>

<tr>

<td style="text-align:left;">

4

</td>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Limestone

</td>

<td style="text-align:right;">

0.221

</td>

<td style="text-align:right;">

0.029

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

0.141

</td>

<td style="text-align:right;">

0.301

</td>

<td style="text-align:left;">

ab

</td>

</tr>

<tr>

<td style="text-align:left;">

2

</td>

<td style="text-align:left;">

City

</td>

<td style="text-align:left;">

Chalk

</td>

<td style="text-align:right;">

0.268

</td>

<td style="text-align:right;">

0.029

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

0.188

</td>

<td style="text-align:right;">

0.347

</td>

<td style="text-align:left;">

b

</td>

</tr>

</tbody>

</table>

</div>

``` r
(mod_BP %>% 
  anova() %>% 
  mutate(`Part Eta Sq`=`Sum Sq`/sum(`Sum Sq`) ) ->
  mod_BP_ANOVA)
```

<div class="kable-table">

<table>

<thead>

<tr>

<th style="text-align:right;">

Df

</th>

<th style="text-align:right;">

Sum Sq

</th>

<th style="text-align:right;">

Mean Sq

</th>

<th style="text-align:right;">

F value

</th>

<th style="text-align:right;">

Pr(\>F)

</th>

<th style="text-align:right;">

Part Eta Sq

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.064

</td>

<td style="text-align:right;">

0.064

</td>

<td style="text-align:right;">

12.52

</td>

<td style="text-align:right;">

0.002

</td>

<td style="text-align:right;">

0.347

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.013

</td>

<td style="text-align:right;">

0.013

</td>

<td style="text-align:right;">

2.59

</td>

<td style="text-align:right;">

0.123

</td>

<td style="text-align:right;">

0.072

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:right;">

0.00

</td>

<td style="text-align:right;">

0.994

</td>

<td style="text-align:right;">

0.000

</td>

</tr>

<tr>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

0.108

</td>

<td style="text-align:right;">

0.005

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

0.582

</td>

</tr>

</tbody>

</table>

</div>

``` r
# pwpp(marginal) # Pairwise P-value plot. Fails for unbalanced design
emmip(mod_BP, Location ~ Rock.type)
```

![](BRC_growth_rate_figures/test%20alpha-10.svg)<!-- -->

``` r
# summary(as.glht(pairs(marginal))) # fails because of unbalanced design
```

#### Plot all alpha diversity metrics together

``` r
Richness_Diversity_long %>% 
  dplyr::filter(!Metric %in% c("Chao1", "ACE")) %>% 
    mutate_at(., "Metric", ~fct_recode(., "Observed S" = "S obs.", "Inv. Simpson" = "Inv. Simpson", "Berger Parker" = "Berger Parker")) %>% 
  mutate_at(., "Metric", ~fct_relevel(., "Observed S", "Inv. Simpson", "Shannon", "Berger Parker")) %>% 
  droplevels() ->
  Richness_Diversity_long2plot

p_alpha <- ggplot() +
  geom_violin(data = Richness_Diversity_long2plot,
             aes(
               x = Location,
               y = Estimate,
               ymin = lerr,
               ymax = herr
             ), colour = "grey",
              fill = "grey",
              alpha = 1 / 3) +
  geom_jitter2(data = Richness_Diversity_long2plot,
               aes(x = Location,
               y = Estimate,
               ymin = lerr,
               ymax = herr,
               colour = Location), size = 3, width = 0.2, alpha = 2/3) +
  scale_colour_manual(values = Gradient.colours[c(5, 6, 11)], name = "") +
  # geom_errorbar(alpha = 1 / 2, width = 0.3) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.9,
    hjust = 1
  ),
  legend.position="none") +
  facet_grid(Metric ~ Rock.type, scale = "free") +
  background_grid(major = "y",
                  minor = "none",
                  size.major = 0.8) 

dat_text <- data.frame(
  label = str_remove_all(c(obsS_pairwise$.group[1:4], 
                           Shannon_pairwise$.group[1:4], 
                           InvSim_pairwise$.group[1:4],
                           BP_pairwise$.group[1:4]), 
                         pattern = " "),
  Metric = fct_inorder(rep(levels(Richness_Diversity_long2plot$Metric), each = 4)),
  Rock.type = fct_c(obsS_pairwise$Rock.type[1:4], 
                           Shannon_pairwise$Rock.type[1:4], 
                           InvSim_pairwise$Rock.type[1:4],
                           BP_pairwise$Rock.type[1:4]), 
  x = fct_c(obsS_pairwise$Location[1:4], 
                           Shannon_pairwise$Location[1:4], 
                           InvSim_pairwise$Location[1:4],
                           BP_pairwise$Location[1:4]),
  # x     = as.factor(levels(Richness_Diversity_long2plot$Climate.Source)),
  y = rep(c(520, 45, 5.2, 0.52), each = 4)
  # y = rep(c(40, 140, 0.5), each = 6)
)
p_alpha <- p_alpha + geom_text(
  data = dat_text,
  mapping = aes(x = x, y = y, label = label),
  nudge_x = 0,
  nudge_y = 0
)
print(p_alpha)
```

![](BRC_growth_rate_figures/plot%20alpha-1.svg)<!-- -->

### Beta diversity

Calculate and plot beta diversity metrics.

``` r
(mod1 <-  adonis(
  otu_table(Ps_obj_filt_GMPR) ~ Location * Rock.type,
  data = as(sample_data(Ps_obj_filt_GMPR), "data.frame"),
  method = "horn",
  permutations = 9999
))
```

    ## 
    ## Call:
    ## adonis(formula = otu_table(Ps_obj_filt_GMPR) ~ Location * Rock.type,      data = as(sample_data(Ps_obj_filt_GMPR), "data.frame"), permutations = 9999,      method = "horn") 
    ## 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                    Df SumsOfSqs MeanSqs F.Model    R2 Pr(>F)  
    ## Location            1      0.63   0.628    2.11 0.081  0.023 *
    ## Rock.type           1      0.39   0.386    1.30 0.050  0.245  
    ## Location:Rock.type  1      0.49   0.485    1.63 0.063  0.098 .
    ## Residuals          21      6.25   0.297         0.806         
    ## Total              24      7.75                 1.000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(mod2 <- adonis(
  otu_table(Ps_obj_filt_GMPR) ~ Location,
  data = as(sample_data(Ps_obj_filt_GMPR), "data.frame"),
  method = "horn",
  permutations = 9999
))
```

    ## 
    ## Call:
    ## adonis(formula = otu_table(Ps_obj_filt_GMPR) ~ Location, data = as(sample_data(Ps_obj_filt_GMPR),      "data.frame"), permutations = 9999, method = "horn") 
    ## 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs MeanSqs F.Model    R2 Pr(>F)  
    ## Location   1      0.63   0.628    2.03 0.081  0.029 *
    ## Residuals 23      7.12   0.309         0.919         
    ## Total     24      7.75                 1.000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mod1_pairwise <- PairwiseAdonis(
  otu_table(Ps_obj_filt_GMPR),
  sample_data(Ps_obj_filt_GMPR)$Location,
  sim.function = "vegdist",
  sim.method = "horn",
  p.adjust.m = "BH"
)
print(mod1_pairwise)
```

    ##          pairs total.DF F.Model     R2 p.value p.adjusted sig
    ## 1 City - Slope       24    2.03 0.0811   0.034      0.034   .

``` r
(sig_pairs1 <- list(as.character(mod1_pairwise$pairs[mod1_pairwise$p.adjusted < 0.05])))
```

    ## [[1]]
    ## [1] "City - Slope"

``` r
mod2_pairwise <- PairwiseAdonis(
  otu_table(Ps_obj_filt_GMPR),
  sample_data(Ps_obj_filt_GMPR)$Rock.type,
  sim.function = "vegdist",
  sim.method = "horn",
  p.adjust.m = "BH"
)
print(mod2_pairwise)
```

    ##               pairs total.DF F.Model     R2 p.value p.adjusted sig
    ## 1 Chalk - Limestone       24    1.21 0.0501    0.32       0.32

``` r
(sig_pairs2 <- list(as.character(mod2_pairwise$pairs[mod2_pairwise$p.adjusted < 0.05])))
```

    ## [[1]]
    ## character(0)

``` r
mod3_pairwise <- PairwiseAdonis(
  otu_table(Ps_obj_filt_GMPR),
  sample_data(Ps_obj_filt_GMPR)$Location.rock,
  sim.function = "vegdist",
  sim.method = "horn",
  p.adjust.m = "BH"
)
print(mod3_pairwise)
```

    ##                              pairs total.DF F.Model     R2 p.value p.adjusted sig
    ## 1      City:Chalk - City:Limestone       11    1.20 0.1069   0.344      0.413    
    ## 2     City:Chalk - Slope:Limestone       11    1.43 0.1251   0.180      0.270    
    ## 3         City:Chalk - Slope:Chalk       12    2.92 0.2099   0.017      0.102    
    ## 4 City:Limestone - Slope:Limestone       11    1.07 0.0964   0.429      0.429    
    ## 5     City:Limestone - Slope:Chalk       12    1.98 0.1528   0.091      0.182    
    ## 6    Slope:Limestone - Slope:Chalk       12    1.90 0.1471   0.080      0.182

``` r
(sig_pairs3 <- list(as.character(mod3_pairwise$pairs[mod3_pairwise$p.adjusted < 0.1])))
```

    ## [[1]]
    ## character(0)

``` r
Unifrac_mat <- UniFrac(Ps_obj_filt, 
                       weighted = TRUE, 
                       normalized = TRUE, 
                       parallel = TRUE, 
                       fast = TRUE)

(mod1 <-  adonis(
  Unifrac_mat ~ Location * Rock.type,
  data = as(sample_data(Ps_obj_filt_GMPR), "data.frame"),
  method = "horn",
  permutations = 9999
))
```

    ## 
    ## Call:
    ## adonis(formula = Unifrac_mat ~ Location * Rock.type, data = as(sample_data(Ps_obj_filt_GMPR),      "data.frame"), permutations = 9999, method = "horn") 
    ## 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                    Df SumsOfSqs MeanSqs F.Model    R2 Pr(>F)  
    ## Location            1     0.047  0.0469    1.37 0.053  0.214  
    ## Rock.type           1     0.053  0.0529    1.54 0.060  0.147  
    ## Location:Rock.type  1     0.063  0.0632    1.84 0.072  0.068 .
    ## Residuals          21     0.719  0.0342         0.815         
    ## Total              24     0.882                 1.000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(mod2 <- adonis(
  Unifrac_mat ~ Location,
  data = as(sample_data(Ps_obj_filt_GMPR), "data.frame"),
  method = "horn",
  permutations = 9999
))
```

    ## 
    ## Call:
    ## adonis(formula = Unifrac_mat ~ Location, data = as(sample_data(Ps_obj_filt_GMPR),      "data.frame"), permutations = 9999, method = "horn") 
    ## 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs MeanSqs F.Model    R2 Pr(>F)
    ## Location   1     0.047  0.0469    1.29 0.053   0.26
    ## Residuals 23     0.835  0.0363         0.947       
    ## Total     24     0.882                 1.000

The difference between city and slope is significant based on the
Morisita-Horn distances between OTUs but not based on UniFrac.

##### Calculate ordinations

``` r
Ps_obj_ord1 <- ordinate(Ps_obj_filt_GMPR, "CAP", "horn", formula = Ps_obj_filt_GMPR ~  Location * Rock.type)
Ps_obj_ord2 <- ordinate(Ps_obj_filt_GMPR, "CAP", "horn", formula = Ps_obj_filt_GMPR ~  Location)

explained <- eigenvals(Ps_obj_ord2)/sum( eigenvals(Ps_obj_ord2)) * 100
explained <- as.numeric(format(round(explained, 1), nsmall = 1))

Ps_obj_filt_GMPR %>% 
  plot_ordination(., Ps_obj_ord2, type = "samples", shape = "Rock.type", color = "Location", justDF = TRUE) -> 
  ord_df

p_ord <- ggplot(ord_df,
             aes(
               x = CAP1,
               y = MDS1,
               shape = Rock.type,
               color = Location
             )) +
  stat_ellipse(
    aes(x = CAP1, 
        y = MDS1, 
        fill = Location
        ),
    geom = "polygon",
    alpha = 1/4,
    type = "t",
    level = 0.9,
    # linetype = 2,
    inherit.aes = FALSE
  ) +
  geom_point(size = 4, alpha = 2 / 3) +
  guides(colour = guide_legend(title = "Location"), shape = guide_legend(title = "Rock.type")) +
  scale_colour_manual(values = Gradient.colours) +
  scale_fill_manual(values = Gradient.colours, guide = "none") +
  labs(x = sprintf("CAP1 (%s%%)", explained[1]), 
       y = sprintf("CAP2 (%s%%)", explained[2])) +
  coord_fixed(ratio = sqrt(explained[2] / explained[1])) +
   theme(legend.justification = "top")
  # facet_wrap(. ~ Rock.type)
print(p_ord)
```

![](BRC_growth_rate_figures/ordinate%20all-1.svg)<!-- -->

``` r
Ps_obj_ord1 <- ordinate(Ps_obj_filt, "PCoA", "Unifrac", formula = Ps_obj_filt ~ Location * Rock.type)
Ps_obj_ord2 <- ordinate(Ps_obj_filt, "PCoA", "Unifrac", formula = Ps_obj_filt ~ Location)

explained <- Ps_obj_ord2$values$Relative_eig/sum(Ps_obj_ord2$values$Relative_eig) * 100
explained <- as.numeric(format(round(explained, 1), nsmall = 1))

Ps_obj_filt %>% 
  plot_ordination(., Ps_obj_ord2, type = "samples", shape = "Rock.type", color = "Location", justDF = TRUE) -> 
  ord_df

p_ord_phylo <- ggplot(ord_df,
             aes(
               x = Axis.1,
               y = Axis.2,
               shape = Rock.type,
               color = Location
             )) +
  stat_ellipse(
    aes(x = Axis.1, 
        y = Axis.2, 
        fill = Location
        ),
    geom = "polygon",
    alpha = 1/4,
    type = "t",
    level = 0.9,
    # linetype = 2,
    inherit.aes = FALSE
  ) +
  geom_point(size = 4, alpha = 2 / 3) +
  theme_bw(base_size = 14) +
  guides(colour = guide_legend(title = "Location"), shape = guide_legend(title = "Rock.type")) +
  scale_colour_manual(values = Gradient.colours) +
  scale_fill_manual(values = Gradient.colours, guide = "none") +
  labs(x = sprintf("CAP1 (%s%%)", explained[1]), 
       y = sprintf("CAP2 (%s%%)", explained[2])) +
  coord_fixed(ratio = sqrt(explained[2] / explained[1])) #+ 
  # facet_wrap(. ~ Rock.type)
print(p_ord_phylo)
```

![](BRC_growth_rate_figures/ordinate%20phylo-1.svg)<!-- -->

#### Test differences between samples on the phylum level

STAMPR analysis of the differences of each phylum between locations
using Aligned Rank Transformed ANOVA test and a post-hoc estimated
marginal means.

``` r
Taxa_tests_phylum1 <- STAMPR2(Ps_obj_filt, vars2test = "Location", threshold = 0.001, outputfile = paste0(Proj_name, "_Location"))

pSTAMPR1 <- plotSTAMPR(Taxa_tests_phylum1, pair = "Slope - City")
print(pSTAMPR1)
```

![](BRC_growth_rate_figures/STAMPR-1.svg)<!-- -->

``` r
Taxa_tests_phylum2 <- STAMPR2(Ps_obj_filt, vars2test = c("Location", "Rock.type"), threshold = 0.001, outputfile = paste0(Proj_name, "_Location_Rock"))

pSTAMPR2 <- plotSTAMPR(Taxa_tests_phylum2, pair = "Slope:Chalk - City:Chalk")
print(pSTAMPR2)
```

![](BRC_growth_rate_figures/STAMPR-2.svg)<!-- -->

### Taxonmic distribution analysis

Agglomerate data and tag rare taxa

``` r
Ps_obj_filt_GMPR_ra <- transform_sample_counts(Ps_obj_filt_GMPR, function(x){x / sum(x)} * 100)

Ps_obj_filt_GMPR_glom <- tax_glom(Ps_obj_filt_GMPR_ra, 
                             "Phylum", 
                             NArm = TRUE)
Ps_obj_filt_GMPR_glom_DF <- speedyseq::psmelt(Ps_obj_filt_GMPR_glom)
Ps_obj_filt_GMPR_glom_DF$Phylum %<>% as.character()
# Ps_obj_filt3_glom_DF %<>% mutate(Species = fct_relevel(Species, "NA", after = Inf))

# group dataframe by Phylum, calculate median rel. abundance
Ps_obj_filt_GMPR_glom_DF %>%
  group_by(Phylum) %>%
  summarise(median = median(Abundance)) ->
  medians

# find Phyla whose rel. abund. is less than 1%
Rare_phyla0.01 <- medians[medians$median <= 0.01, ]$Phylum

# change their name to "Rare"
Ps_obj_filt_GMPR_glom_DF[Ps_obj_filt_GMPR_glom_DF$Phylum %in% Rare_phyla0.01, ]$Phylum <- 'Rare'
# re-group
Ps_obj_filt_GMPR_glom_DF %>%
  group_by(Sample, Phylum, Location, Rock.type, Location.rock) %>%
  summarise(Abundance = sum(Abundance)) ->
  Ps_obj_filt_GMPR_glom_DF_2plot

# ab.taxonomy$Freq <- sqrt(ab.taxonomy$Freq)
# Ps_obj_filt3_glom_rel_DF$Phylum %<>% sub("unclassified", "Unclassified", .)
# Ps_obj_filt3_glom_rel_DF$Phylum %<>% sub("uncultured", "Unclassified", .)

Ps_obj_filt_GMPR_glom_DF_2plot %>% 
  group_by(Sample) %>% 
  filter(Phylum == "Rare") %>% 
  summarise(`Rares (%)` = sum(Abundance)) -> 
  Rares
```

Summarise taxonomy

``` r
# Percentage of reads classified as rare 
Rares %>%
  kable(., digits = c(2), caption = "Percentage of reads per sample classified as rare:") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Percentage of reads per sample classified as rare:

</caption>

<thead>

<tr>

<th style="text-align:left;">

Sample

</th>

<th style="text-align:right;">

Rares (%)

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

ShivSClk1S12

</td>

<td style="text-align:right;">

0.01

</td>

</tr>

<tr>

<td style="text-align:left;">

ShivSClk2S15

</td>

<td style="text-align:right;">

0.59

</td>

</tr>

<tr>

<td style="text-align:left;">

ShivSClk3S16

</td>

<td style="text-align:right;">

0.02

</td>

</tr>

<tr>

<td style="text-align:left;">

ShivSClk4S17

</td>

<td style="text-align:right;">

0.01

</td>

</tr>

<tr>

<td style="text-align:left;">

ShivSClk5S18

</td>

<td style="text-align:right;">

0.01

</td>

</tr>

<tr>

<td style="text-align:left;">

ShivSLimek1S13

</td>

<td style="text-align:right;">

0.02

</td>

</tr>

<tr>

<td style="text-align:left;">

ShivSLimek2S14

</td>

<td style="text-align:right;">

0.04

</td>

</tr>

<tr>

<td style="text-align:left;">

ShivSLimek3S19

</td>

<td style="text-align:right;">

14.27

</td>

</tr>

<tr>

<td style="text-align:left;">

ShivSLimek4S20

</td>

<td style="text-align:right;">

11.48

</td>

</tr>

<tr>

<td style="text-align:left;">

ShivSLimek5S21

</td>

<td style="text-align:right;">

0.01

</td>

</tr>

<tr>

<td style="text-align:left;">

ShivSLimek6S25

</td>

<td style="text-align:right;">

0.01

</td>

</tr>

<tr>

<td style="text-align:left;">

SlpClk1S3

</td>

<td style="text-align:right;">

0.24

</td>

</tr>

<tr>

<td style="text-align:left;">

SlpClk2S4

</td>

<td style="text-align:right;">

0.04

</td>

</tr>

<tr>

<td style="text-align:left;">

SlpClk3S5

</td>

<td style="text-align:right;">

0.04

</td>

</tr>

<tr>

<td style="text-align:left;">

SlpClk4S6

</td>

<td style="text-align:right;">

0.03

</td>

</tr>

<tr>

<td style="text-align:left;">

SlpClk5S7

</td>

<td style="text-align:right;">

0.01

</td>

</tr>

<tr>

<td style="text-align:left;">

SlpClk6S8

</td>

<td style="text-align:right;">

0.03

</td>

</tr>

<tr>

<td style="text-align:left;">

SlpClk7S9

</td>

<td style="text-align:right;">

0.01

</td>

</tr>

<tr>

<td style="text-align:left;">

SlpClk8S10

</td>

<td style="text-align:right;">

0.03

</td>

</tr>

<tr>

<td style="text-align:left;">

SlpClk9S11

</td>

<td style="text-align:right;">

0.05

</td>

</tr>

<tr>

<td style="text-align:left;">

SlpLime1S1

</td>

<td style="text-align:right;">

0.01

</td>

</tr>

<tr>

<td style="text-align:left;">

SlpLime2S2

</td>

<td style="text-align:right;">

0.04

</td>

</tr>

<tr>

<td style="text-align:left;">

SlpLime3S22

</td>

<td style="text-align:right;">

0.09

</td>

</tr>

<tr>

<td style="text-align:left;">

SlpLime4S23

</td>

<td style="text-align:right;">

0.01

</td>

</tr>

<tr>

<td style="text-align:left;">

SlpLime5S24

</td>

<td style="text-align:right;">

0.02

</td>

</tr>

</tbody>

</table>

``` r
sample_order <- match(Rares$Sample, row.names(sample_data(Ps_obj_filt_GMPR_glom)))
Rares %<>% arrange(., sample_order)

Rares %>% 
  cbind(., sample_data(Ps_obj_filt_GMPR_glom)) %>% 
  group_by(Location.rock) %>% 
  setNames(make.names(names(.), unique = TRUE)) %>% # fails for some reason without it
  summarise(`Rares (%)` = mean(`Rares....`)) -> 
  Rares_merged

# Percentage of reads classified as rare 
Rares_merged %>%
  kable(., digits = c(2), caption = "Percentage of reads per sample classified as rare:") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Percentage of reads per sample classified as rare:

</caption>

<thead>

<tr>

<th style="text-align:left;">

Location.rock

</th>

<th style="text-align:right;">

Rares (%)

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

City:Chalk

</td>

<td style="text-align:right;">

0.11

</td>

</tr>

<tr>

<td style="text-align:left;">

City:Limestone

</td>

<td style="text-align:right;">

4.31

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope:Chalk

</td>

<td style="text-align:right;">

0.03

</td>

</tr>

<tr>

<td style="text-align:left;">

Slope:Limestone

</td>

<td style="text-align:right;">

0.07

</td>

</tr>

</tbody>

</table>

Plot taxonomy box-plot

``` r
Ps_obj_filt_GMPR_glom_DF_2plot %>% 
  group_by(Phylum) %>% 
  summarise(sum.Taxa = sum(Abundance)) %>% 
  arrange(desc(sum.Taxa)) -> Taxa_rank
Ps_obj_filt_GMPR_glom_DF_2plot$Phylum %<>% 
  factor(., levels = Taxa_rank$Phylum) %>% 
  fct_relevel(., "Rare", after = Inf)
  
p_taxa_box <-
  ggplot(Ps_obj_filt_GMPR_glom_DF_2plot, aes(x = Phylum, y = (Abundance))) +
  geom_boxplot(aes(group = interaction(Phylum, Location)), position = position_dodge(width = 0.9), fatten = 1) +
  geom_point(
    aes(colour = Rock.type),
    position = position_jitterdodge(dodge.width = 1),
    alpha = 1 / 2,
    stroke = 0,
    size = 2
  ) +
  scale_colour_manual(values = Gradient.colours, name = "") +
  theme_bw()+
  # theme_cowplot(font_size = 11, font_family = f_name) +
  labs(x = NULL, y = "Relative abundance (%)") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  facet_grid(Location ~ .) +
  background_grid(major = "xy",
                  minor = "none") +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.9,
    hjust = 0.9
  ),
  legend.position = c(.99, .99),
  legend.justification = c("right", "top"),
  legend.box.just = "top",
  legend.margin = margin(0, 3, 3, 3))
print(p_taxa_box)
```

![](BRC_growth_rate_figures/taoxonomy_boxplot-1.svg)<!-- -->

### Differential abundance models

Tag rare phyla (for plotting purposes only)

``` r
Ps_obj_filt_GMPR_glom <- tax_glom(Ps_obj_filt_GMPR, 
                             "Phylum", 
                             NArm = TRUE) # glomerate to the phylum level
Ps_obj_filt_GMPR_glom_rel <- transform_sample_counts(Ps_obj_filt_GMPR_glom, function(x) x / sum(x)) # transform to rel. ab.
Ps_obj_filt_GMPR_glom_rel_DF <- speedyseq::psmelt(Ps_obj_filt_GMPR_glom_rel) # generate a df
Ps_obj_filt_GMPR_glom_rel_DF$Phylum %<>% as.character() # factor to char

# group dataframe by Phylum, calculate median rel. abundance
Ps_obj_filt_GMPR_glom_rel_DF %>%
  group_by(Phylum) %>%
  summarise(median = median(Abundance)) ->
  medians

# find Phyla whose median rel. abund. is less than 0.5%
Rare_phyla0.005 <- medians[medians$median <= 0.005, ]$Phylum

# change their name to "Rare"
Ps_obj_filt_GMPR_glom_rel_DF[Ps_obj_filt_GMPR_glom_rel_DF$Phylum %in% Rare_phyla0.005, ]$Phylum <- 'Rare'

# re-group
Ps_obj_filt_GMPR_glom_rel_DF %>%
  group_by(Phylum) %>%
  summarise(Abundance = sum(Abundance)) %>% 
  arrange(desc(Abundance)) -> Taxa_rank
```

Detect differentially abundant OTUs using corncob (Martin, Witten, and
Willis [2020](#ref-martin_modeling_2020))

``` r
comparison_string <- c("City", "Slope")

Ps_obj_filt %>%
  subset_samples(Location %in% c(comparison_string[1], comparison_string[2])) %>%
  tax_glom("Order") ->
  Ps_obj_filt_pairwise_glom

# Test differential abundance for location
da_Loc <- differentialTest(formula = ~ Location,
                           phi.formula = ~ Location,
                           formula_null = ~ 1,
                           phi.formula_null = ~ Location, 
                           test = "Wald", boot = FALSE,
                           data = Ps_obj_filt,
                           fdr_cutoff = 0.05,
                           full_output = TRUE)
da_Loc_intervals <- plot(da_Loc, level = "Class", data_only = T)
which(is.na(da_Loc$p)) %>% names
```

    ## [1] "OTU398" "OTU164" "OTU476" "OTU641" "OTU239" "OTU481" "OTU727" "OTU451" "OTU245"

``` r
Ps_obj_filt %>%
  transform_sample_counts(., function(x) x / sum(x) * 100) %>% 
  taxa_sums(.) %>% 
  map_dbl(~(.x / nsamples(Ps_obj_filt))) %>% 
  enframe(name = "OTU", value = "Mean abundance (%)") -> 
  baseMean

map(da_Loc$all_models,15) %>% 
  map(.,2) %>% 
  unlist %>%  # grab all mu.LocationSlope Estimates (differences in estimated population relative abundance)
  bind_cols(OTU = taxa_names(Ps_obj_filt), 
            tax_table(Ps_obj_filt), 
            `Differential abundance` = .,
            Significance = fct_recode(as_factor(taxa_names(Ps_obj_filt) %in% da_Loc$significant_taxa), Pass = "TRUE", Fail = "FALSE"),
            ymin = as.numeric(NA),
            ymax = as.numeric(NA)
            ) %>%
  left_join(., baseMean, by = "OTU") ->
  da_Loc_df

da_Loc_df %<>% rows_update(., tibble(ymin = da_Loc_intervals$xmin, OTU = da_Loc$significant_taxa), by = "OTU")
da_Loc_df %<>% rows_update(., tibble(ymax = da_Loc_intervals$xmax, OTU = da_Loc$significant_taxa), by = "OTU")
da_Loc_df[da_Loc_df$Phylum %in% Rare_phyla0.005, "Phylum"] <- 'Rare' # rare_phyla is

p_corncob_loc <- GGPlotCorncob(da_Loc_df, OTU_labels = FALSE, Taxa = "Phylum", Y_val = "Differential abundance", sig_level = 0.05, Rank = Taxa_rank)

corncob_summary <- tibble(Label = c(paste0("⬆", sum(da_Loc_df$`Differential abundance` > 0 &  da_Loc_df$Significance == "Pass"), " ⬇", sum(da_Loc_df$`Differential abundance` < 0 &  da_Loc_df$Significance == "Pass"), " (", nrow(da_Loc_df), ")")))

p_corncob_loc <- p_corncob_loc +
  geom_text(
    data    = corncob_summary,
    mapping = aes(x = Inf, y = Inf, label = Label),
    hjust   = 1.1,
    vjust   = 1.6
  ) +
  scale_size_continuous(name = "Mean abundance (%)",
                        limits = c(min(da_Loc_df$`Mean abundance (%)`), max(da_Loc_df$`Mean abundance (%)`))
  ) + 
  labs(title = paste(comparison_string, collapse = " - ")) +
  coord_cartesian(ylim = c(-15, 15))
print(p_corncob_loc)
```

![](BRC_growth_rate_figures/corncob%20-%20location-1.svg)<!-- -->

Modelling differential abundance and variance between locations
discovered `length(da_Loc$significant_taxa)`

``` r
comparison_string <- c("Limestone", "Chalk")
# Test differential abundance and variance for rock type
da_Rock <- differentialTest(formula = ~ Rock.type,
                                 phi.formula = ~ Rock.type,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ Rock.type, 
                                 test = "Wald", boot = FALSE,
                                 data = Ps_obj_filt,
                                 fdr_cutoff = 0.05,
                                full_output = TRUE)
da_Rock_intervals <- plot(da_Rock, level = "Class", data_only = TRUE)
which(is.na(da_Rock$p)) %>% names
```

    ##  [1] "OTU780" "OTU164" "OTU580" "OTU679" "OTU696" "OTU712" "OTU519" "OTU503" "OTU542"
    ## [10] "OTU451" "OTU114" "OTU309" "OTU502" "OTU245" "OTU229"

``` r
map(da_Rock$all_models,15) %>% 
  map(.,2) %>% 
  unlist %>%  # grab all mu.LocationSlope Estimates (differences in estimated population relative abundance)
  bind_cols(OTU = taxa_names(Ps_obj_filt), 
            tax_table(Ps_obj_filt), 
            `Differential abundance` = .,
            Significance = fct_recode(as_factor(taxa_names(Ps_obj_filt) %in% da_Rock$significant_taxa), Pass = "TRUE", Fail = "FALSE"),
            ymin = as.numeric(NA),
            ymax = as.numeric(NA)
            ) %>%
  left_join(., baseMean, by = "OTU") ->
  da_Rock_df

da_Rock_df %<>% rows_update(., tibble(ymin = da_Rock_intervals$xmin, OTU = da_Rock$significant_taxa), by = "OTU")
da_Rock_df %<>% rows_update(., tibble(ymax = da_Rock_intervals$xmax, OTU = da_Rock$significant_taxa), by = "OTU")
da_Rock_df[da_Rock_df$Phylum %in% Rare_phyla0.005, "Phylum"] <- 'Rare' # rare_phyla is

p_corncob_rock <- GGPlotCorncob(da_Rock_df, OTU_labels = FALSE, Taxa = "Phylum", Y_val = "Differential abundance", sig_level = 0.05, Rank = Taxa_rank)

corncob_summary <- tibble(Label = c(paste0("⬆", sum(da_Rock_df$`Differential abundance` > 0 &  da_Rock_df$Significance == "Pass"), " ⬇", sum(da_Rock_df$`Differential abundance` < 0 &  da_Rock_df$Significance == "Pass"), " (", nrow(da_Rock_df), ")")))

p_corncob_rock <- p_corncob_rock +
  geom_text(
    data    = corncob_summary,
    mapping = aes(x = Inf, y = Inf, label = Label),
    hjust   = 1.1,
    vjust   = 1.6
  ) +
  scale_size_continuous(name = "Mean abundance (%)",
                        limits = c(min(da_Rock_df$`Mean abundance (%)`), max(da_Rock_df$`Mean abundance (%)`))
  ) + 
  labs(title = paste(comparison_string, collapse = " - ")) +
  coord_cartesian(ylim = c(-15, 15))
print(p_corncob_rock)
```

![](BRC_growth_rate_figures/corncob%20-%20rock-1.svg)<!-- -->

Modelling differential abundance and variance between rock types
discovered `length(da_Rock$significant_taxa)`

``` r
# Test differential abundance for location, control for Rock.type for both cases
comparison_string <- c("City", "Slope")
da_Loc_exRock <- differentialTest(formula = ~ Location + Rock.type,
                                 phi.formula = ~ Location + Rock.type,
                                 formula_null = ~ Rock.type,
                                 phi.formula_null = ~ Location + Rock.type, 
                                 test = "Wald", boot = FALSE,
                                 data = Ps_obj_filt,
                                 fdr_cutoff = 0.05,
                                full_output = TRUE)
da_Loc_exRock_intervals <- plot(da_Loc_exRock, level = "Class", data_only = TRUE)

which(is.na(da_Loc_exRock$p)) %>% names
```

    ##  [1] "OTU780" "OTU398" "OTU164" "OTU476" "OTU580" "OTU679" "OTU696" "OTU712" "OTU641"
    ## [10] "OTU519" "OTU239" "OTU481" "OTU503" "OTU542" "OTU727" "OTU451" "OTU114" "OTU309"
    ## [19] "OTU502" "OTU245" "OTU229"

``` r
map(da_Loc_exRock$all_models,15) %>% 
  map(.,2) %>% 
  unlist %>%  # grab all mu.LocationSlope Estimates (differences in estimated population relative abundance)
  bind_cols(OTU = taxa_names(Ps_obj_filt), 
            tax_table(Ps_obj_filt), 
            `Differential abundance` = .,
            Significance = fct_recode(as_factor(taxa_names(Ps_obj_filt) %in% da_Loc_exRock$significant_taxa), Pass = "TRUE", Fail = "FALSE"),
            ymin = as.numeric(NA),
            ymax = as.numeric(NA)
            ) %>%
  left_join(., baseMean, by = "OTU") ->
  da_Loc_exRock_df

da_Loc_exRock_df %<>% rows_update(., tibble(ymin = da_Loc_exRock_intervals$xmin, OTU = da_Loc_exRock$significant_taxa), by = "OTU")
da_Loc_exRock_df %<>% rows_update(., tibble(ymax = da_Loc_exRock_intervals$xmax, OTU = da_Loc_exRock$significant_taxa), by = "OTU")
da_Loc_exRock_df[da_Loc_exRock_df$Phylum %in% Rare_phyla0.005, "Phylum"] <- 'Rare' # rare_phyla is

p_corncob_locExroc <- GGPlotCorncob(da_Loc_exRock_df, OTU_labels = FALSE, Taxa = "Phylum", Y_val = "Differential abundance", sig_level = 0.05, Rank = Taxa_rank)

corncob_summary <- tibble(Label = c(paste0("⬆", sum(da_Loc_exRock_df$`Differential abundance` > 0 &  da_Loc_exRock_df$Significance == "Pass"), " ⬇", sum(da_Loc_exRock_df$`Differential abundance` < 0 &  da_Loc_exRock_df$Significance == "Pass"), " (", nrow(da_Loc_exRock_df), ")")))

p_corncob_locExroc <- p_corncob_locExroc +
  geom_text(
    data    = corncob_summary,
    mapping = aes(x = Inf, y = Inf, label = Label),
    hjust   = 1.1,
    vjust   = 1.6
  ) +
  scale_size_continuous(name = "Mean abundance (%)",
                        limits = c(min(da_Loc_exRock_df$`Mean abundance (%)`), max(da_Loc_exRock_df$`Mean abundance (%)`))
  ) + 
  labs(title = paste(comparison_string, collapse = " - ")) +
  coord_cartesian(ylim = c(-15, 15))
print(p_corncob_locExroc)
```

![](BRC_growth_rate_figures/corncob%20-%20loc%20exl.%20rock-1.svg)<!-- -->

Modelling differential abundance between locations, while controlling
for rock type discovered `length(da_Loc_exRock$significant_taxa)`

``` r
mod260 <- bbdml(formula = OTU260 ~ 1,
             phi.formula = ~ 1,
             data = Ps_obj_filt)
mod260_Loc <- bbdml(formula = OTU260 ~ Location,
             phi.formula = ~ Location,
             data = Ps_obj_filt)
mod260_Loc_rock <- bbdml(formula = OTU97 ~ Location*Rock.type,
             phi.formula = ~ Location*Rock.type,
             data = Ps_obj_filt)
lrtest(mod_null = mod260, mod = mod260_Loc)
```

    ## [1] 0.0029

``` r
# lrtest(mod_null = mod260_Loc, mod = mod260_Loc_rock)
summary(mod260_Loc)
```

    ## 
    ## Call:
    ## bbdml(formula = OTU260 ~ Location, phi.formula = ~Location, data = Ps_obj_filt)
    ## 
    ## 
    ## Coefficients associated with abundance:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   -8.562      0.239  -35.78   <2e-16 ***
    ## Location1      0.676      0.239    2.82     0.01 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## 
    ## Coefficients associated with dispersion:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   -8.535      0.409  -20.89  1.6e-15 ***
    ## Location1      1.427      0.409    3.49   0.0022 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## 
    ## Log-likelihood: -80.58

``` r
plot(mod260_Loc, color = "Location", shape = "Rock.type") # add total = TRUE for total counts (i.e. not relative abundance)
```

![](BRC_growth_rate_figures/plot%20OTU%20260-1.svg)<!-- -->

### Compose figures

``` r
composite_plot <- ((p_alpha + p_taxa_box +  plot_layout(widths = c(1, 2))) /(p_ord + pSTAMPR1) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = f_size)))

plot_file <- "./Results/Microbiome_1"
svglite(paste0(plot_file, ".svg"), 
        width = 14, 
        height = 10.5)
print(composite_plot)
dev.off()
```

    ## pdf 
    ##   2

``` r
ggsave(
  paste0(plot_file, ".png"),
  composite_plot,
  device = agg_png,
  width = 14, 
  height = 10.5, 
  units = "cm", 
  res = 900,
  scaling = 0.38
)
gz(paste0(plot_file, ".svg"), paste0(plot_file, ".svgz"))
knitr::include_graphics(paste0(plot_file, ".png"))
```

<img src="./Results/Microbiome_1.png" width="1587" />

``` r
plot_file <- "./Results/Microbiome_2"
svglite(paste0(plot_file, ".svg"), 
        width = 12, 
        height = 10)
print(p_corncob_locExroc)
dev.off()
```

    ## pdf 
    ##   2

``` r
ggsave(
  paste0(plot_file, ".png"),
  p_corncob_locExroc,
  device = agg_png,
  width = 12, 
  height = 10, 
  units = "cm", 
  res = 900,
  scaling = 0.38
)
gz(paste0(plot_file, ".svg"), paste0(plot_file, ".svgz"))
```

### References

<div id="refs" class="references">

<div id="ref-chen_gmpr:_2017">

Chen, Jun, and Li Chen. 2017. “GMPR: A novel normalization method for
microbiome sequencing data.” *bioRxiv*, February, 112565.
<https://doi.org/10.1101/112565>.

</div>

<div id="ref-martin_modeling_2020">

Martin, Bryan D., Daniela Witten, and Amy D. Willis. 2020. “Modeling
Microbial Abundances and Dysbiosis with Beta-Binomial Regression.”
*Annals of Applied Statistics* 14 (1): 94–115.
<https://doi.org/10.1214/19-AOAS1283>.

</div>

</div>
