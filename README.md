# 2020-FatBodyQueens-Bimpatiens

## Contains all analyses performed on datasets for transcriptome analysis in the fat body of pre-overwintering bumble bee queens (Costa et al. 2020)

### Synopsis

We used RNA sequencing to examine how age and nectar diet (specifically, the concentration of sucrose in nectar) impact gene expression in the pre-overwintering bumble bee queen fat body, the “liver-like” organ in insects with broad functions related to nutrient storage and metabolism.

### Sample Information

| ID_Seq	| ID_Queen	| Colony	| Age	| Diet (% sucrose solution) |
| --------| ----------|---------| ----|---------------------------|
| Q648_R47_0d |	Q648	| R47	| 0	| # |
| Q675_R59_0d	|	Q675	|	R59	|	0	| # |
| Q680_R53_0d	|	Q680	|	R53	|	0	| # |
| Q903_R47_0d	|	Q903	|	R47	|	0	| # |
| Q736_R53_3d0	|	Q736	|	R53	|	3	| 0 |
| Q755_R53_3d0	|	Q755	|	R53	|	3	| 0 |
| Q828_R47_3d0	|	Q828	|	R47	|	3	| 0 |
| Q853_R47_3d0	|	Q853	|	R47	|	3	| 0 |
| Q666_R45_3d25	|	Q666	|	R45	|	3	| 25 |
| Q829_R47_3d25	|	Q829	|	R47	|	3	| 25 |
| Q831_R59_3d25	|	Q831	|	R59	|	3	| 25 |
| Q852_R49_3d25	|	Q852	|	R49	|	3	| 25 |
| Q625_R59_3d50	|	Q625	|	R59	|	3	| 50 |
| Q669_R53_3d50	|	Q669	|	R53	|	3	| 50 |
| Q878_R47_3d50	|	Q878	|	R47	|	3	| 50 |
| Q934_R53_3d50	|	Q934	|	R53	|	3	| 50 |
| Q825_R59_3d75	|	Q825	|	R59	|	3	| 75 |
| Q839_R47_3d75	|	Q839	|	R47	|	3	| 75 |
| Q899_R45_3d75	|	Q899	|	R45	|	3	| 75 |
| Q931_R47_3d75	|	Q931	|	R47	|	3	| 75 |
| Q765_R53_9d0	|	Q765	|	R53	|	9	| 0 |
| Q820_R45_9d0	|	Q820	|	R45	|	9	| 0 |
| Q867_R47_9d0	|	Q867	|	R47	|	9	| 0 |
| Q879_R59_9d0	|	Q879	|	R59	|	9	| 0 |
| Q767_R59_9d25	|	Q767	|	R59	|	9	| 25 |
| Q785_R51_9d25	|	Q785	|	R51	|	9	| 25 |
| Q815_R45_9d25	|	Q815	|	R45	|	9	| 25 |
| Q921_R53_9d25	|	Q921	|	R53	|	9	| 25 |
| Q588_R56_9d50	|	Q588	|	R56	|	9	| 50 |
| Q787_R45_9d50	|	Q787	|	R45	|	9	| 50 |
| Q869_R59_9d50	|	Q869	|	R59	|	9	| 50 |
| Q888_R47_9d50	|	Q888	|	R47	|	9	| 50 |
| Q633_R56_9d75	|	Q633	|	R56	|	9	| 75 |
| Q781_R59_9d75	|	Q781	|	R59	|	9	| 75 |
| Q788_R45_9d75	|	Q788	|	R45	|	9	| 75 |
| Q922_R53_9d75	|	Q922	|	R53	|	9	| 75 |
| Q623_R53_12d0	|	Q623	|	R53	|	12	| 0 |
| Q706_R45_12d0	|	Q706	|	R45	|	12	| 0 |
| Q793_R59_12d0	|	Q793	|	R59	|	12	| 0 |
| Q836_R47_12d0	|	Q836	|	R47	|	12	| 0 |
| Q436_R59_12d25	|	Q436	|	R59	|	12	| 25 |
| Q438_R45_12d25	|	Q438	|	R45	|	12	| 25 |
| Q806_R59_12d25	|	Q806	|	R59	|	12	| 25 |
| Q837_R47_12d25	|	Q837	|	R47	|	12	| 25 |
| Q441_R59_12d50	|	Q441	|	R59	|	12	| 50 |
| Q595_R53_12d50	|	Q595	|	R53	|	12	| 50 |
| Q707_R59_12d50	|	Q707	|	R59	|	12	| 50 |
| Q902_R53_12d50	|	Q902	|	R53	|	12	| 50 |
| Q437_R59_12d75	|	Q437	|	R59	|	12	| 75 |
| Q464_R47_12d75	|	Q464	|	R47	|	12	| 75 |
| Q749_R53_12d75	|	Q749	|	R53	|	12	| 75 |
| Q866_R59_12d75	|	Q866	|	R59	|	12	| 75 |


_NOTES: Newly eclosed callow queens, which are less than 24 hr of age, were removed from their natal colonies and either collected immediately into liquid nitrogen (hereafter referred to as “0 day” group) or subjected to one of four diet treatments, where they were given 0%, 25%, 50% or 75% (w/v) sucrose solution; these are hereafter referred to as “sucrose solution” of various concentrations._


### Alignment details

In brief, the workflow:

````
1. Sequence data quality was evaluated using FASTQC;
2. Sequencing adapters and low quality bases were trimmed/filtered using Trimmomatic;
3. Reads were mapped to the B. impatiens v2.0 reference genome (Sadd et al. 2015) using HiSat2;
	 - Reference genome index was created using OGS annotation (splice sites & exons)
4. Exonic expression (read counts) was extracted and summed for genes using featureCount.
````

_see [readmapping.sh](https://github.com/claudinpcosta/2020-FatBodyQueens-Bimpatiens/blob/master/readmapping.sh) for additional information_

### Identification of differentially expressed genes

Differential gene expression was analysed using the DESeq2 package. The 0 day group was excluded from this factorial analysis because these queens were collected on the day of eclosion and did not receive a diet treatment. Expression data from the 0 day group were analysed for a subset of genes. Read counts were transformed using the variance stabilizing transformation with diet, age, colony, or age × diet as factors in the design for sample clustering and visualization of gene expression levels.

_see [DESEq2code.R](https://github.com/claudinpcosta/2020-FatBodyQueens-Bimpatiens/blob/master/DESEq2code.R) for the codes, and we can find the input in the folder [Data](https://github.com/claudinpcosta/2020-FatBodyQueens-Bimpatiens/tree/master/Data)_











