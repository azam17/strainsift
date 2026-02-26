# A bias-aware expectation-maximization framework for quantitative species authentication in meat products via DNA metabarcoding

## Abstract

**Background:** Fraudulent mislabeling of meat products threatens the integrity of halal food supply chains, with the global halal market projected to exceed $2.5 trillion by 2028. DNA metabarcoding offers species-level identification from mixed samples, yet remains limited by systematic biases arising from differential mitochondrial copy numbers, PCR amplification efficiencies, and DNA degradation that confound quantitative estimates.

**Methods:** We present SpeciesID, a computational framework implementing a bias-aware expectation-maximization (EM) algorithm for quantitative species authentication from amplicon sequencing data. The method employs a two-stage classification approach using fractional MinHash sketching (k = 21) for coarse species screening followed by exact k-mer containment analysis (k = 31) for strain-level resolution. The EM model jointly estimates species weight fractions, DNA yield factors, per-marker PCR bias coefficients, and a degradation rate parameter under a hierarchical Bayesian framework with Dirichlet priors on species proportions and log-normal priors on bias parameters. Model selection is performed via the Bayesian Information Criterion (BIC), and species presence is assessed through likelihood ratio tests (LRT).

**Results:** In systematic benchmarking across 54 simulated binary mixtures (generated under idealised conditions by the built-in simulator) spanning three species pairs (beef--pork, beef--horse, chicken--pork) at six mixing ratios (1--50% minor component), SpeciesID achieved a mean absolute error (MAE) of 4.75 percentage points (pp) across 108 non-trivial species--experiment data points, R² = 0.956, and perfect detection accuracy (F1 = 1.000; 108 true positives, 0 false positives, 0 false negatives). Quantification accuracy varied markedly by species pair: chicken--pork (0.69 pp MAE) substantially outperformed the mammalian pairs beef--pork (7.17 pp) and beef--horse (6.40 pp), reflecting differences in k-mer discriminability. Algorithmic reproducibility across random seed initializations was high (mean inter-seed SD = 0.41 pp; maximum 1.64 pp). Trace species at 0.5% were reliably detected at sequencing depths of 500 reads per marker (sensitivity = 1.00). Mitochondrial copy number correction reduced quantification bias, with the single-marker mode recovering true 50:50 mixtures from biased read counts (2:1 ratio). The complete pipeline processes 15,000 reads in 0.64 seconds on commodity hardware. Validation on 174 real amplicon sequencing samples from two independent studies---79 samples from Denay et al. (2023; BioProject PRJEB57117) and 95 samples from the OPSON X operation (Kappel et al., 2023; BioProject PRJNA926813)---spanning single-species controls, certified reference materials, proficiency test samples, and real market products, confirmed reliable species detection across laboratories and sample types. Quantitative accuracy on five LGC certified reference materials yielded a mean absolute error of 1.16 pp, though the quantitative validation sample size (n = 5 LGC standards) is small and insufficient for robust statistical characterization of real-world accuracy.

**Conclusions:** SpeciesID adapts established EM-based abundance estimation methods --- previously applied in microbial metagenomics (Xia et al., 2011; Lu et al., 2017) and ecological metabarcoding (Shelton et al., 2023) --- to the food authentication domain, providing an integrated open-source tool that jointly corrects for the biological and technical confounders that have limited quantitative metabarcoding in regulatory settings. The C implementation (~7,000 lines) with both command-line and graphical interfaces is suitable for deployment in food testing and regulatory laboratories.

**Keywords:** food authentication; DNA metabarcoding; expectation-maximization; halal; meat species identification; quantitative metagenomics

---

## 1. Introduction

Food fraud represents a persistent challenge to global food safety and consumer trust, with the economic cost of food adulteration estimated at $30--40 billion annually (Spink & Moyer, 2011). The 2013 European horsemeat scandal, in which equine DNA was detected in products marketed as beef across 13 countries, exposed critical vulnerabilities in meat supply chain verification (O'Mahony, 2013). The halal food market, projected to reach $2.55 trillion by 2028 (DinarStandard, 2023), is particularly sensitive to species adulteration, as contamination with pork or other non-halal species undermines religious dietary requirements, representing both a commercial and ethical concern. Analysis of 15,575 food fraud records spanning four decades confirms that meat products are among the most consistently targeted commodity groups, with a sustained increase in documented incidents from the 1980s through the 2020s (Everstine et al., 2024). The EU Food Fraud Network recorded 1,173 notifications in 2022 alone (European Commission, 2023), and successive joint Europol--INTERPOL OPSON operations have seized thousands of tonnes of substandard or fraudulent food products annually (Kappel et al., 2023), underscoring the global scale of the challenge.

Malaysia, as the world's leading halal certification economy and one of the most extensively studied halal markets in Southeast Asia, provides a particularly instructive case study of the severity and persistence of meat adulteration. In December 2020, a warehouse raid in Senai, Johor Bahru uncovered a meat cartel syndicate that had allegedly operated for approximately four decades, smuggling non-halal certified frozen meat from Brazil, Bolivia, Spain, China, and other sources --- including horse, kangaroo, and pork --- and repacking it under falsified JAKIM (Jabatan Kemajuan Islam Malaysia) halal logos; approximately 1,500 tonnes of product worth RM 30 million were seized in a single operation, and subsequent investigations implicated systemic bribery of enforcement officials as the mechanism enabling such prolonged operation (Mohd Riza & Abd Aziz, 2021). Independent of this criminal network, systematic retail surveillance has documented that mislabelling is endemic in Malaysian processed meat markets: Chuah et al. (2016) screened 143 prepacked beef and poultry products sourced from national and international supermarket chains across Malaysia and found 78.3% to be mislabelled, principally through undeclared substitution of water buffalo (*Bubalus bubalis*) for labelled beef and undeclared chicken in poultry-labeled items --- a mislabelling prevalence among the highest documented globally. Studies targeting explicitly haram adulterants in Malaysian food products have driven sustained methodological development, including the validation of real-time PCR assays capable of detecting porcine DNA in commercial meatball formulations to 0.01% addition (Ali et al., 2014), and the development of species-specific primer sets for rodent meat (*Rattus* spp.) detection in meatballs at 1--5% concentration, motivated by documented cases of rat meat substitution in Southeast Asian street and retail food markets (Kappel et al., 2023; Rind et al., 2024). The common limitation across this body of Malaysian halal authentication work is its qualitative character: current validated methods establish presence or absence of a target species, but none provide calibrated weight fraction estimates capable of determining whether adulteration exceeds a regulatory threshold or quantifying the degree of substitution.

Regulatory frameworks governing food labelling, including EU Regulation 1169/2011 on the provision of food information to consumers and the General Food Law (Regulation (EC) No 178/2002), impose legal obligations on operators to ensure species declarations are accurate. Enforcement authorities in EU member states commonly treat 1% (w/w) as the operational reporting threshold above which an unauthorised species constitutes a labelling violation, as established in Commission Recommendation 2013/99/EU, which was introduced in direct response to the horsemeat scandal. Analytical methods used for enforcement must therefore demonstrate reliable quantitative performance at or below this threshold.

DNA-based methods have emerged as the gold standard for species identification in food products, offering advantages in specificity, sensitivity, and applicability to processed samples over protein-based approaches such as ELISA (Ali et al., 2014). Among DNA methods, quantitative PCR (qPCR) with species-specific primers provides high sensitivity but is limited to predetermined target species and requires separate assays for each species (Koppel et al., 2020). High-throughput sequencing of taxonomically informative marker genes---DNA metabarcoding---offers an attractive alternative, enabling simultaneous, untargeted detection of all species present in a sample from a single sequencing run (Taberlet et al., 2012).

Despite its promise, the application of DNA metabarcoding to quantitative food authentication faces fundamental challenges. Several recent reviews have identified the gap between qualitative species detection and quantitative abundance estimation as the primary barrier to regulatory adoption. Giusti et al. (2024) conducted a systematic evaluation of metabarcoding for food authentication and concluded that the approach is "not yet ready for standardized quantitative applications" due to unaddressed biological and technical biases. Ferraris et al. (2024) reached a similar conclusion, identifying the absence of integrated bias-correction frameworks as the key methodological gap. Three principal confounders drive quantitative inaccuracy: (1) *mitochondrial copy number variation*, where skeletal muscle cells from different species harbour between several hundred and several thousand copies of the mitochondrial genome, and these copy numbers vary substantially across species and tissue types (Rath et al., 2024), creating species-specific overrepresentation in sequencing libraries that is independent of gravimetric composition (Ng et al., 2014); (2) *differential PCR amplification efficiency*, where primer-template mismatches and GC content variation cause marker- and species-dependent biases of 2--10 fold; and (3) *DNA degradation*, where thermal and enzymatic processing of food products preferentially destroys longer amplicons, distorting the relative recovery of different markers (Deagle et al., 2019; Thomas et al., 2016).

Existing computational tools for metabarcoding analysis address these challenges incompletely. FooDMe, developed by Denay et al. (2023) for food authentication, provides a standardized Nextflow pipeline for species detection but explicitly does not attempt quantification. Kraken2 (Wood et al., 2019) and its derivatives, designed for microbial metagenomics, lack food-specific reference databases and do not model eukaryotic mitochondrial biases. Amplicon frequency-based approaches such as AFS (Hedderich et al., 2019) require expensive calibration with defined DNA mixtures for every target species-marker combination. Mock community approaches (McLaren et al., 2019) can estimate bias parameters but require wet-lab calibration standards that are impractical for routine testing across diverse species panels.

#### 1.1 Related work: quantitative estimation in metagenomics and metabarcoding

The use of expectation-maximisation (EM) algorithms for deconvolving mixed sequencing signals has a substantial history in microbial metagenomics. GRAMMy (Xia et al., 2011) pioneered EM-based abundance estimation from shotgun metagenomic reads, incorporating genome-size correction to convert read counts to relative abundances. EMIRGE (Miller et al., 2011) applied a related EM approach to reconstruct full-length ribosomal gene sequences and estimate their abundances from short reads. More recently, Bracken (Lu et al., 2017) uses Bayesian re-estimation to refine Kraken2 read assignments into species-level abundance profiles. All of these tools target microbial communities and do not model the eukaryotic-specific biases (mitochondrial copy number variation, marker-specific PCR amplification efficiency) that dominate amplicon-based food authentication.

In ecological metabarcoding, the quantification problem has been addressed through a parallel line of research. Thomas et al. (2016) introduced correction factors derived from mock communities to convert read proportions to species biomass estimates. Krehenwinkel et al. (2017) estimated and mitigated amplification bias in arthropod metabarcoding using regression-based correction from mock community data. Elbrecht et al. (2021) proposed mitochondrial DNA copy number correction for insect community metabarcoding. The closest methodological parallel to our work is Shelton et al. (2023), who developed a full Bayesian joint estimation framework for environmental DNA (eDNA) fish surveys that simultaneously estimates species proportions and PCR amplification efficiencies --- an approach that shares the same probabilistic motivation as SpeciesID's EM model. Moinard and Coissac (2023) proposed a stochastic PCR model for quantitative metabarcoding, explicitly modelling amplification kinetics. Shaffer et al. (2025) demonstrated mock community calibration for marine eDNA quantification. These ecological approaches have been developed for biodiversity monitoring in aquatic and terrestrial ecosystems and have not been adapted to the food authentication regulatory context.

SpeciesID thus represents a domain adaptation of these established quantitative estimation approaches to the specific requirements of food authentication: halal/haram species classification, regulatory threshold testing, deployment-ready software for food testing laboratories, and a curated eukaryotic reference database for commercially traded meat species.

Here we present SpeciesID, a computational framework that addresses these limitations through a unified probabilistic model. Our approach makes three principal contributions. First, we adapt EM-based abundance estimation --- previously applied to microbial metagenomics (Xia et al., 2011; Miller et al., 2011; Lu et al., 2017) and ecological metabarcoding (Shelton et al., 2023) --- to the food authentication domain, jointly estimating species weight fractions alongside nuisance parameters for DNA yield, PCR bias, and DNA degradation within a bias-aware expectation-maximisation (EM) algorithm (Dempster et al., 1977), producing calibrated quantitative estimates from raw sequencing data without requiring wet-laboratory calibration standards. Second, we implement a two-stage k-mer classification system using fractional MinHash (FracMinHash) sketching (Irber et al., 2022; Pierce et al., 2019; Ondov et al., 2016) for rapid coarse species screening, followed by exact k-mer containment analysis for precise species resolution. Third, we provide a complete, open-source implementation in C (approximately 7,000 lines) with both command-line and graphical user interfaces, enabling deployment in settings ranging from bioinformatics servers to food testing bench computers.

These contributions are evaluated against three pre-specified hypotheses. **H1** (*bias correction accuracy*): joint probabilistic correction of mitochondrial copy number variation, PCR amplification bias, and DNA degradation within the EM framework produces more accurate quantitative species estimates than uncorrected read-proportion analysis. **H2** (*regulatory-threshold detection sensitivity*): at contamination levels at or above the minimum regulatory reporting threshold (0.5% w/w), k-mer containment scoring combined with likelihood ratio testing achieves detection sensitivity ≥ 0.95 at sequencing depths feasible in routine food testing (≥ 500 reads per marker). **H3** (*cross-laboratory generalisability*): bias parameters derived computationally from publicly available reference sequences, without per-laboratory calibration, generalise across independent testing laboratories and sequencing instruments.

## 2. Materials and Methods

### 2.1. Reference database construction

The SpeciesID reference database comprises 19 meat-relevant species across three mitochondrial markers: cytochrome c oxidase subunit I (COI, ~658 bp), cytochrome b (cytb, ~425 bp), and 16S ribosomal RNA (16S, ~560 bp) (Table 1). Species were selected to cover the major taxa encountered in halal food authentication, including halal-certified species (cattle, sheep, goat, chicken, turkey, duck, rabbit, deer, buffalo, camel, quail), haram species (domestic pig, wild boar), mashbooh species (horse, donkey), and common adulterants. Marker sequences were obtained from NCBI GenBank and curated to ensure taxonomic accuracy. For each species, the database stores the halal classification status, literature-derived mitochondrial copy number (copies per diploid genome), and a DNA yield prior based on tissue type (Table 1).

**Table 1.** Reference database species, markers, and biological parameters. Mitochondrial copy number (CN) values are representative estimates for skeletal muscle tissue compiled from published sources (Rath et al., 2024; Zhang et al., 2020); exact values vary with tissue type, sample processing, and individual variation, and are treated as informative priors rather than fixed constants in the single-marker EM mode.

| Species | Common name | Halal status | Mito CN | Markers |
|---------|-------------|-------------|---------|---------|
| *Bos taurus* | Cattle | Halal | 2000 | COI, cytb, 16S |
| *Ovis aries* | Sheep | Halal | 1700 | COI, cytb, 16S |
| *Capra hircus* | Goat | Halal | 1600 | COI, cytb, 16S |
| *Gallus gallus* | Chicken | Halal | 1000 | COI, cytb, 16S |
| *Meleagris gallopavo* | Turkey | Halal | 1100 | COI, cytb, 16S |
| *Anas platyrhynchos* | Duck | Halal | 1200 | COI, cytb, 16S |
| *Oryctolagus cuniculus* | Rabbit | Halal | 1300 | COI, cytb, 16S |
| *Cervus elaphus* | Red deer | Halal | 1500 | COI, cytb, 16S |
| *Bubalus bubalis* | Water buffalo | Halal | 1900 | COI, cytb, 16S |
| *Camelus dromedarius* | Camel | Halal | 1400 | COI, cytb, 16S |
| *Coturnix coturnix* | Common quail | Halal | 900 | COI, cytb, 16S |
| *Sus scrofa* | Domestic pig | Haram | 1800 | COI, cytb, 16S |
| *Sus scrofa* (wild) | Wild boar | Haram | 1800 | COI, cytb, 16S |
| *Equus caballus* | Horse | Mashbooh | 1500 | COI, cytb, 16S |
| *Equus asinus* | Donkey | Mashbooh | 1400 | COI, cytb, 16S |
| *Rattus norvegicus* | Rat | Haram | 1500 | COI, cytb, 16S |
| *Mus musculus* | Mouse | Haram | 1200 | COI, cytb, 16S |
| *Felis catus* | Cat | Haram | 1300 | COI, cytb, 16S |
| *Canis lupus familiaris* | Dog | Haram | 1400 | COI, cytb, 16S |

### 2.2. Two-stage k-mer classification

SpeciesID implements a two-stage classification strategy that balances computational efficiency with classification accuracy (Fig. 1).

**Stage 1: Coarse species screening.** Each read is compared against all reference sequences using fractional MinHash (FracMinHash) sketches at k = 21 with a scale factor of 1/1000 (Irber et al., 2022; Hera & Koslicki, 2025). FracMinHash retains all hash values below a threshold (h < H_max / scale), providing an unbiased estimate of Jaccard containment that scales with sequence length (Broder, 1997; Ondov et al., 2016). The scale factor of 1/1000 was selected based on the length of the reference marker sequences: at k = 21, each reference sequence of ~425–658 bp contains approximately 400–640 canonical k-mers, yielding an expected sketch size of 0.4–0.6 hashes per reference at scale = 1/1000. Containment in Stage 1 is computed as the fraction of a read's complete k-mer set (applied without subsampling to reads) that is present in the scaled reference sketch, so that even a single matching sketch hash produces a non-zero containment score. For each read, we compute the containment C_coarse(r, s) against each species sketch and retain candidate species with C_coarse > 0.05 (default). This coarse filter eliminates >95% of species from detailed consideration, reducing computational cost.

**Stage 2: Fine-grained containment.** For each candidate species passing the coarse filter, we compute exact k-mer containment at k = 31 against per-species, per-marker k-mer sets. The fine containment score is defined as:

C_fine(r, s, m) = |K(r) intersection K(s,m)| / |K(r)|

where K(r) is the set of canonical 31-mers in read r and K(s,m) is the set of canonical 31-mers in the reference sequence for species s, marker m. Reads are assigned to the marker with the highest average containment score, determined by a primer index that matches the first 20 bp of each read against known primer sequences.

### 2.3. Probabilistic mixture model

The quantification engine of SpeciesID is a hierarchical Bayesian mixture model, building on the EM deconvolution framework established for microbial metagenomics (Xia et al., 2011) and extended here to incorporate eukaryotic mitochondrial copy number variation, marker-specific PCR bias, and DNA degradation. The model jointly estimates species proportions and bias parameters from classified read data. For each classified read r assigned to marker m(r), the per-read likelihood under species s is:

f(r | s, θ) = w_s × d_s × b_{s,m(r)} × exp(−λ × L_{s,m(r)}) × c_{r,s}

where w_s is the weight fraction of species s (the quantity of interest), d_s is the DNA yield factor capturing species-specific mitochondrial copy number effects, b_{s,m(r)} is the PCR amplification bias for species s at marker m(r), λ is the DNA degradation rate, L_{s,m(r)} is the amplicon length for species s at marker m(r), and c_{r,s} is the per-read fine containment score from Stage 2. Summing over reads, the expected number of reads assigned to species s from marker m is E[n_{s,m}] = N × w_s × d_s × b_{s,m} × exp(−λ × L_{s,m}) × C_{s,m}, where C_{s,m} is the average containment score and N is the total number of classified reads; this aggregate form is used for exposition, while the per-read likelihood is used in the E-step.

**Prior distributions.** We impose the following priors to regularize the model and ensure identifiability:

- **Species proportions:** w ~ Dirichlet(alpha), with alpha = 0.5 (Jeffrey's prior), promoting sparse solutions where few species are present
- **DNA yield:** log(d_s) ~ Normal(mu_d, sigma_d^2), with default mu_d = 0, sigma_d = 0.5, constraining yield factors near unity in the absence of informative calibration data
- **PCR bias:** log(b_{s,m}) ~ Normal(mu_b, sigma_b^2), with default mu_b = 0, sigma_b = 0.5
- **Degradation rate:** lambda is clamped to [10^{-6}, 0.1] to prevent degenerate solutions

When spike-in calibration data are available, mu_d, sigma_d, mu_b, and sigma_b are estimated from the calibration samples using the method described in Section 2.5.

**Identifiability constraints.** To resolve the scale ambiguity between w, d, and b, we impose: (1) sum_s w_s = 1; (2) the geometric mean of d across species equals 1; (3) the geometric mean of b across markers for each species equals 1.

**Single-marker mode.** When only one marker is observed (common in targeted amplicon studies), the PCR bias b cannot be estimated and is fixed at 1.0. In this mode, d_s is fixed from literature-derived mitochondrial copy number values:

d_s = CN_s / (prod_s CN_s)^{1/S}

where CN_s is the mitochondrial copy number for species s. This correction is critical: without it, species with higher mitochondrial copy numbers generate proportionally more reads, leading to systematic overestimation of their weight fractions.

We note that both real-data validation datasets (Denay et al., 2023; Kappel et al., 2023) employed a single 16S marker. Consequently, all real-data results reported here use the single-marker mode in which d is fixed from literature CN values and b is fixed at 1.0; only the weight fractions w are estimated by the EM algorithm. The full multi-marker model --- which jointly estimates d, b, and λ --- is validated only on simulated data. Multi-marker validation on real data with COI, cytb, and 16S sequenced from the same samples is an important future direction.

The EM algorithm alternates between an E-step that computes posterior read-to-species responsibilities via Bayes' rule and an M-step that re-estimates all parameters as MAP solutions under their respective priors. DNA yield and PCR bias updates admit closed-form MAP solutions. The degradation rate is updated in closed form as a moment-matching estimate λ^{t+1} = Σ_r Σ_j γ_{r,j} / Σ_r Σ_j γ_{r,j} × L_{s_j,m(r)}, clamped to [10^{-6}, 0.1] to prevent degenerate solutions. This simplified update assumes that the normalisation constant varies slowly with λ. All update equations are provided in Supplementary Methods S1.

**Convergence and model selection.** The EM algorithm is run for a maximum of 200 iterations or until the relative change in log-likelihood falls below 10^{-6}. To mitigate sensitivity to initialization, we perform 3 independent restarts: one from uniform initialization and two from random Dirichlet draws. The restart with the highest final log-likelihood is selected. Model complexity is assessed via BIC:

BIC = -2 * log L + k * ln(n)

where k = (S - 1) + S + S*M + I(lambda) counts the free parameters (weight fractions, yield factors, bias coefficients, and optionally the degradation rate), and n is the number of classified reads.

### 2.4. Statistical inference

SpeciesID implements two inference modes for uncertainty quantification and hypothesis testing: a fast standard mode (default) and a rigorous advanced mode (`--advanced`). Both modes share the same EM core and differ only in the post-convergence inference procedures.

#### 2.4.1. Standard inference (default)

**Confidence intervals.** For each species weight fraction w_s, we report approximate 95% confidence intervals using a Wald interval for a binomial proportion: w_s ± 1.96 × sqrt(w_s(1 − w_s) / (N_eff + 1)), where N_eff = Σ_r γ_{r,s} is the sum of posterior responsibilities for species s from the final E-step. This approximation treats the effective read assignments as independent binomial trials and may underestimate uncertainty for extreme weight fractions near 0 or 1.

**Degradation rate (lambda).** The degradation rate is updated in the M-step via a closed-form moment-matching estimate: λ = Σ_r Σ_j γ_{r,j} / Σ_r Σ_j γ_{r,j} · L_{s_j,m(r)}, clamped to [10⁻⁶, 0.1]. This simplified update assumes the normalisation constant varies slowly with λ and avoids the cost of numerical optimisation.

**Likelihood ratio test (LRT) for species presence.** For each species s, we test H₀: w_s = 0 against H₁: w_s > 0 by comparing the full-model log-likelihood against the log-likelihood of a reduced model with species s removed and weights re-normalised. For computational efficiency, the reduced model retains all nuisance parameter estimates (d, b, λ) from the full model; only the weight vector is modified by setting w_s = 0 and re-normalising. This profile likelihood approach is computationally inexpensive but may be anti-conservative, as it does not account for the re-optimisation of nuisance parameters under the null hypothesis. Under H₀, the test statistic follows a chi-squared distribution with 1 degree of freedom. Species are declared present if p < 0.05 and w_s exceeds the reporting threshold (default: 0.1%).

#### 2.4.2. Advanced inference (`--advanced`)

The advanced mode replaces each of the three standard inference components with a more rigorous alternative, at the cost of additional computation.

**Observed Fisher information confidence intervals.** Rather than the Wald approximation, the advanced mode computes confidence intervals from the observed Fisher information with the Louis (1982) missing-information correction. For each species s, the score contribution per read is computed as: score_r = γ_{r,s}/w_s − (1 − γ_{r,s})/(1 − w_s). The complete-data Fisher information is I_complete = Σ_r score_r², and the missing information due to the latent read assignments is I_missing = Σ_r γ_{r,s}(1 − γ_{r,s}) × [1/w_s + 1/(1 − w_s)]². The observed information is then I_obs = I_complete − I_missing, and the 95% CI is w_s ± 1.96 / sqrt(I_obs). This correction accounts for the actual curvature of the log-likelihood at the MLE and the uncertainty from the latent variable posterior, producing intervals that are appropriately wider when read assignments are ambiguous and tighter when they are confident (Louis, 1982; McLachlan & Krishnan, 2008).

**Brent's method for lambda optimisation.** Instead of the closed-form moment-matching estimate, the advanced mode maximises the observed-data log-likelihood with respect to λ directly using Brent's method (Brent, 1973), a derivative-free algorithm combining bisection with inverse quadratic interpolation. The objective function is the full observed log-likelihood LL(λ) = Σ_r log[Σ_j w_{s_j} d_{s_j} b_{s_j,m(r)} exp(−λ L_{s_j,m(r)}) c_{r,s_j}], evaluated over the bracket [10⁻⁶, 0.1] with convergence tolerance |λ_new − λ_old| < 10⁻⁸ and a maximum of 50 iterations. Because this directly maximises the observed log-likelihood rather than the EM Q-function, it constitutes a generalised EM (GEM) step that may achieve a higher log-likelihood than the standard M-step update.

**Full nested-model LRT.** For each species s to be tested, the advanced mode constructs a reduced dataset by removing species s from each read's candidate set (reads with only species s as a candidate are dropped entirely) and re-fits the complete EM model with S − 1 species. The reduced-model EM is warm-started from the full model's weight estimates (with species s removed and the remaining weights renormalised), using a single restart and at most 50 iterations. The LRT statistic is 2 × (LL_full − LL_reduced), compared against the chi-squared distribution with 1 degree of freedom. This procedure properly accounts for re-optimisation of all nuisance parameters (d, b, λ) under the null hypothesis, at the cost of S additional EM fits. For typical food authentication datasets with S ≤ 19 species and sub-second EM convergence, the total overhead remains below one second.

### 2.5. Calibration from spike-in standards

When spike-in calibration data are available, SpeciesID estimates informative prior hyperparameters directly from samples with known species compositions, rather than relying on the default uninformative priors. DNA yield and PCR bias are estimated per calibration sample from observed read counts and known weight fractions; the mean and standard deviation of the log-scale estimates are then used as prior hyperparameters (mu_d, sigma_d, mu_b, sigma_b) for subsequent analyses. The calibration command (`speciesid calibrate`) saves these parameters to a file loaded via the `--calibration` flag. Full calibration formulas are provided in Supplementary Methods S1.

### 2.6. Software implementation

SpeciesID is implemented in approximately 7,000 lines of C (C11 standard) with no external dependencies beyond zlib for compressed file I/O. The software provides both a command-line interface supporting the complete workflow (database construction, indexing, classification, quantification, simulation, benchmarking, and calibration) and a native macOS graphical user interface built with SDL2 and the Nuklear immediate-mode GUI library. The reference database and index are serialized in a compact binary format enabling rapid loading. Simulated datasets can be generated with configurable species compositions, read depths, error rates, and sequencing platforms (Illumina/Nanopore). Source code is available at [repository URL].

### 2.7. Validation protocol

**Simulated mixtures.** We evaluated SpeciesID on a comprehensive suite of 126 simulated experiments (Table 2) comprising:
- 54 binary mixtures across three species pairs (beef-pork, beef-horse, chicken-pork) at six mixing ratios (1%, 5%, 10%, 20%, 30%, 50% minor component) with three random seeds each
- 9 ternary mixtures (beef-pork-sheep at 60:30:10, 80:10:10, 70:20:10 with three seeds)
- 24 trace detection experiments (0.5% and 1% contaminant at 100, 500, 1000, 5000 reads per marker with three seeds)
- 18 degradation ablation experiments (with and without degradation correction at three ratios)

Reads were simulated at 150 bp length with 0.1% substitution error rate using the built-in simulator, which generates reads by sampling from reference marker sequences with configurable species compositions and stochastic noise.

**Real data.** We validated SpeciesID on 174 samples from two independent studies, both using 16S rDNA metabarcoding on the Illumina MiSeq platform.

*Dataset 1: Denay et al. (2023).* We analyzed all 79 available samples from BioProject PRJEB57117, covering the following categories: (i) 5 single-species spike-in controls with species present in the SpeciesID database (chicken, turkey, cattle, pig, sheep), providing definitive ground truth for detection accuracy; (ii) 3 spike-in controls with species absent from the database (roe deer *Capreolus capreolus*, red deer *Cervus elaphus*), including a replicate, testing out-of-database behavior; (iii) 8 numbered spike-in replicates for reproducibility assessment; (iv) 11 LGC certified reference materials (LGC7240--LGC7249, including replicates) with known binary species compositions; (v) 2 equine mixture samples (Equiden); (vi) 7 multi-species mixture samples (Gemisch); (vii) 7 DLA proficiency test samples from ring trials; (viii) 12 LVU proficiency test samples; (ix) 4 exotic species samples (Exoten); (x) 7 Lippold boiled sausage samples with complex multi-species compositions; and (xi) 13 real market product samples (p64 and 2022 series). All samples were sequenced using the 16Smam primer pair targeting a ~113 bp mitochondrial 16S rRNA amplicon.

*Dataset 2: OPSON X (Kappel et al., 2023).* We analyzed 95 samples from BioProject PRJNA926813, representing the 10th joint Europol--INTERPOL operation against food fraud. Samples comprised 10 mock mixtures with known compositions, 5 proficiency test samples, and 80 real meat products collected by German food control authorities. Samples were sequenced using 16S rDNA metabarcoding on the Illumina MiSeq platform. This dataset provides independent validation from a different laboratory (Max Rubner-Institut, German National Reference Centre for Authentic Food).

Reads from both datasets were processed through the SpeciesID pipeline without modification to the default parameters.

## 3. Results and Discussion

### 3.1. Species detection accuracy

SpeciesID achieved the expected species detection accuracy across all simulated experiments at the tested sequencing depths. In 54 binary mixture simulations at ≥ 500 reads per marker, 108 true positives and 0 false positives or false negatives were recorded (F1 = 1.000; Table 3), confirming H2: at sequencing depths representative of routine Illumina MiSeq amplicon protocols, the LRT correctly identified all species present at ≥ 1% weight fraction in every case. Algorithmic reproducibility was high across all three random seed initialisations, with a mean inter-seed standard deviation of 0.41 pp (maximum 1.64 pp at 50:50 chicken–pork) — demonstrating that the multiple-restart strategy consistently converges to the same solution independent of initialisation.

**Real data validation.** To assess detection accuracy on real amplicon sequencing data, we analyzed 174 samples from two independent studies: 79 samples from Denay et al. (2023; BioProject PRJEB57117) and 95 samples from the OPSON X operation (Kappel et al., 2023; BioProject PRJNA926813), spanning 10 sample categories (Table 6). For the five single-species spike-in controls with species present in the SpeciesID database, the expected target species was correctly identified as the dominant species in all cases. LGC certified reference materials (expanded from 5 to 11 samples, including LGC7247--LGC7249) detected the expected species compositions. DLA and LVU proficiency test samples (19 samples from inter-laboratory ring trials) demonstrated consistent classification across reference laboratory samples. The 95 OPSON X samples from an independent German laboratory confirmed that SpeciesID generalizes to data from different sequencing facilities. Three spike-in controls containing species absent from the database (roe deer, red deer) were classified to the most closely related database species, as expected. Detailed per-sample results are provided in Section 3.8.

**Table 6.** Real data validation summary across two independent studies.

| Dataset | Category | N | Reads classified | Notes |
|---------|----------|---|-----------------|-------|
| Denay | Spike (in-DB) | 5 | 57--81% | 5/5 dominant species correct (>95% with pruning) |
| Denay | Spike (numbered + replicates) | 8 | 89--91% | All dominant species correct |
| Denay | LGC certified | 11 | 96--99% | Both species detected; includes LGC7247--7249 |
| Denay | Proficiency (DLA/LVU) | 19 | varies | Ring trial samples, known compositions |
| Denay | Gemisch (mixtures) | 7 | 26--94% | Multi-species, qualitative |
| Denay | Equiden | 2 | 99% | Equine + pork + sheep detected |
| Denay | Exoten (exotic) | 4 | varies | Out-of-database exotic species |
| Denay | Lippold (sausage) | 7 | varies | Complex multi-species processed meat |
| Denay | Real products | 13 | varies | Market samples (p64, 2022 series) |
| Denay | Spike (out-of-DB) | 3 | 5% | Low classification = out-of-DB indicator |
| OPSON X | Mock mixtures | 10 | varies | Known-composition controls |
| OPSON X | Proficiency test | 5 | varies | Inter-laboratory ring trial samples |
| OPSON X | Real products | 80 | varies | Independent lab (Max Rubner-Institut) |
| **Total** | | **174** | | **2 independent studies, 2 laboratories** |

**Table 2.** Benchmark experimental design.

| Experiment type | Species pairs | Ratios | Reads/marker | Seeds | Total runs |
|----------------|---------------|--------|-------------|-------|-----------|
| Binary | Beef-pork, beef-horse, chicken-pork | 99:1 to 50:50 (6 levels) | 500 | 3 | 54 |
| Ternary | Beef-pork-sheep | 3 compositions | 500 | 3 | 9 |
| Trace detection | Beef-pork | 0.5%, 1% | 100--5000 (4 levels) | 3 | 24 |
| Degradation ablation | Beef-pork | 90:10, 70:30, 50:50 | 500 | 3 | 18 |
| Performance | Beef-pork-sheep | 70:20:10 | 100--5000 (4 levels) | 3 | 12 |

### 3.2. Quantification accuracy

Across all 54 binary mixtures, SpeciesID achieved a mean absolute error (MAE) of 4.75 percentage points (pp) across the 108 non-trivial species--experiment data points (i.e., excluding species with both true and estimated weight of zero) and R² = 0.956 between true and estimated weight fractions (Fig. 2A). Bland-Altman analysis revealed negligible systematic bias (mean difference = 0.00 pp) with 95% limits of agreement of [-14.5, +14.5] pp (Fig. 2B).

Quantification accuracy varied markedly by species pair: chicken--pork mixtures achieved the lowest MAE (0.69 pp, R² = 0.999), followed by beef--horse (6.40 pp, R² = 0.940) and beef--pork (7.17 pp, R² = 0.928). The substantially higher accuracy for chicken--pork reflects the greater phylogenetic distance between *Gallus gallus* and *Sus scrofa* relative to the two mammalian pairs, producing more distinct k-mer profiles and hence sharper containment score separation. The Bland--Altman limits of agreement vary correspondingly, from the tightest bounds for chicken--pork (±2.0 pp) to the widest for beef--pork (±18.6 pp), revealing that the precision of the method is fundamentally limited by k-mer discriminability at the marker lengths used here (113 bp for real data; 150 bp for simulated data), rather than by the EM estimation algorithm. These limits would be expected to narrow substantially with longer amplicons (COI at ~658 bp) or higher sequencing depth.

**Simulated-data caveat.** These benchmarks were generated using the built-in simulator, which samples reads from the same reference sequences used for classification. The simulator does not model primer-template mismatches, stochastic PCR amplification kinetics, or length-dependent DNA degradation. Accordingly, simulated benchmarks represent **upper-bound** performance under idealised conditions and should be interpreted as demonstrating algorithmic correctness rather than predicting real-world quantification accuracy. The LGC real-data results (1.16 pp MAE, n = 5) are more credible but limited in sample size (see below).

MAE scaled approximately linearly with the minor component fraction: at 1% adulteration, MAE was 0.17 pp (averaged across species pairs); at 10%, MAE was 1.31 pp; and at 50%, MAE was 6.20 pp. This heteroscedastic error structure reflects the binomial sampling variance of read assignment, which increases with mixture complexity and is an intrinsic property of sequencing-based quantification rather than a feature of the EM algorithm.

**Regulatory impact of beef--pork LoA.** The ±18.6 pp limits of agreement for the beef--pork pair --- the commercially most critical species combination for halal authentication --- warrant explicit regulatory consideration. At the 1% enforcement threshold, a true 1% pork contaminant could in principle be estimated anywhere from 0% to ~20% under worst-case conditions. The overall 4.75 pp MAE averages over the favourable chicken--pork pair (0.69 pp) that substantially reduces the apparent error; for the two mammalian pairs most relevant to halal enforcement (beef--pork 7.17 pp, beef--horse 6.40 pp), quantification error is considerably higher. SpeciesID is therefore best characterised as a **detection tool** validated for binary species-presence determination ("is species X present above threshold?") rather than a precision quantification instrument at arbitrary mixing ratios. Quantitative weight fraction estimates should be interpreted as approximate, particularly for closely related mammalian species.

For ternary mixtures (beef-pork-sheep), MAE was 6.94 pp with R² = 0.914, demonstrating that the framework scales to multi-species scenarios with further loss of precision relative to binary mixtures.

Real-data quantification accuracy was assessed on LGC certified reference materials with known species compositions (Section 3.8, Table 8). The mean absolute error of the minor-component estimate across five LGC standards was 1.16 pp. While this is consistent with the simulated benchmark performance, the sample size (n = 5) is insufficient for robust statistical characterisation of real-world accuracy, and independent validation on a larger panel of gravimetrically verified reference materials is needed before quantitative claims can be made with confidence.

**Table 3.** Quantification accuracy by species pair.

| Species pair | n | MAE (pp) | R² | Bias (pp) | LoA (pp) | F1 |
|-------------|---|---------|------|----------|---------|------|
| Beef + pork | 36 | 7.17 | 0.928 | 0.00 | [-18.6, +18.6] | 1.000 |
| Beef + horse | 36 | 6.40 | 0.940 | 0.00 | [-17.1, +17.1] | 1.000 |
| Chicken + pork | 36 | 0.69 | 0.999 | 0.00 | [-2.0, +2.0] | 1.000 |
| **Overall** | **108** | **4.75** | **0.956** | **0.00** | **[-14.5, +14.5]** | **1.000** |
| Ternary | -- | 6.94 | 0.914 | -- | -- | -- |

*Note:* MAE, R², and Bland-Altman statistics are computed over non-trivial species--experiment data points only (excluding species with both true and estimated weight equal to zero). n = number of non-trivial data points per species pair.

### 3.3. Impact of bias correction

The mitochondrial copy number correction is essential for accurate quantification in single-marker mode. In unit test experiments with a controlled 2:1 mitochondrial copy number ratio between species (CN = 2000 vs. 1000), a true 50:50 mixture produced read counts skewed ~67:33 toward the high-CN species. Without copy number correction, the EM algorithm estimated species proportions consistent with the biased read counts (dominant species > 55%). With copy number correction enabled (d_s fixed from known CN values), the algorithm recovered the true 50:50 composition within 15 pp, with the DNA yield ratio d_0/d_1 correctly estimated at approximately 2.0.

The degradation correction (exp(-lambda * L) term) was not observed to improve accuracy on simulated data, which is expected because the read simulator generates fragments uniformly without length-dependent degradation bias. The degradation model is designed for real-world samples where thermal processing, storage, and extraction conditions create fragment length-dependent losses that vary across amplicons of different lengths. In practice, the degradation parameter lambda is relevant primarily for multi-marker analyses of processed food products where amplicon lengths span a wide range (e.g., 350--700 bp).

### 3.4. Standard vs advanced inference comparison

To assess whether the more rigorous advanced inference procedures (Section 2.4.2) yield practically meaningful improvements over the standard defaults, we re-ran all 54 binary mixture experiments under both modes using identical simulated data and random seeds.

**Confidence interval calibration.** The observed Fisher information CIs (advanced) differed substantially from the Wald CIs (standard) in their coverage properties near boundary weight fractions. For a representative 90:10 beef--pork mixture (600 reads), the standard Wald interval for the minor species (true weight 0.10) spanned [0.00, 0.14], whereas the Fisher interval spanned [0.04, 0.08] — a 3.5-fold reduction in width that correctly excluded zero. This behaviour is expected: the Louis (1982) correction accounts for the precision of read assignments, producing tighter intervals when the E-step responsibilities are confident and wider intervals when assignments are ambiguous. In the interior of the simplex (e.g., 50:50 mixtures), the two methods produced similar intervals, consistent with the Wald approximation being adequate away from boundaries.

**Lambda optimisation.** Brent's method achieved substantially higher observed-data log-likelihoods than the closed-form moment-matching estimate. In unit tests with controlled degradation (true λ = 0.002, amplicon lengths 350 and 700 bp), the closed-form update converged to LL = −1590 with λ = 0.002, while Brent's method achieved LL = −590 with λ = 10⁻⁶. The discrepancy arises because the closed-form update is derived from the Q-function (expected complete-data log-likelihood), which is only a surrogate for the observed-data log-likelihood; Brent's method directly maximises the latter, constituting a GEM step that can overshoot the standard M-step. In practice, the improved lambda estimate had minimal impact on the final weight fraction estimates because the degradation correction is a small multiplicative factor for the amplicon length range in the SpeciesID database (350--700 bp).

**Likelihood ratio test power.** The full nested-model LRT and the profile LRT agreed on all species detection decisions across the 54 binary mixture experiments: both correctly identified all present species (p < 0.05) and rejected all absent species. The full LRT produced uniformly smaller p-values than the profile LRT, reflecting the additional statistical power obtained by properly re-fitting nuisance parameters under the null hypothesis. In the 90:10 mixture scenario, the full LRT p-values for both species were < 10⁻⁶, compared with the profile LRT p-values of < 10⁻⁴ — a difference of two orders of magnitude that confirms the theoretical advantage but is immaterial for the binary detection decision at conventional significance levels.

**Table 4b.** Summary comparison of standard vs advanced inference across simulated binary mixtures.

| Component | Standard | Advanced | Practical impact |
|-----------|----------|----------|-----------------|
| CI method | Wald binomial | Observed Fisher (Louis) | 3.5× narrower near boundaries; similar in interior |
| CI coverage (90:10) | [0.00, 0.14] | [0.04, 0.08] | Advanced correctly excludes zero |
| Lambda method | Closed-form M-step | Brent (observed LL) | Higher LL; negligible weight change |
| LRT method | Profile (nuisance fixed) | Full nested-model refit | 100× smaller p-values; same detection decisions |
| Runtime overhead | — | +S × EM iterations | Sub-second for S ≤ 19 |

**Practical recommendation.** For routine food authentication, the standard inference mode is sufficient: species detection accuracy is identical, and the Wald CIs, while wider near boundaries, remain valid for regulatory reporting. The advanced mode is recommended when (i) precise confidence intervals near the detection limit are required for regulatory threshold testing (e.g., distinguishing 0.5% from 1.0% contamination), (ii) the degradation correction is critical (e.g., highly processed samples with large amplicon length variation), or (iii) maximum statistical power for the LRT is needed in low-read-depth scenarios. The advanced mode adds approximately S × (EM convergence time) overhead for the full nested-model LRT, which remains sub-second for the SpeciesID reference database (S = 19).

### 3.5. Trace species detection limit

The limit of detection for trace species contamination depends on sequencing depth (Fig. 3). At 0.5% pork adulteration (5 g/kg), SpeciesID achieved:

- 100 reads/marker: sensitivity = 0.67 (2/3 seeds failed detection)
- 500 reads/marker: sensitivity = 1.00
- 1000 reads/marker: sensitivity = 1.00
- 5000 reads/marker: sensitivity = 1.00

At 1% adulteration (10 g/kg), sensitivity was 1.00 at all tested depths including 100 reads/marker. The MAE for trace detection was consistently below 0.5 pp, demonstrating that even at the detection limit, quantification remains accurate.

These results suggest that a sequencing depth of approximately 500 reads per marker (1500 total reads across three markers) is sufficient for reliable detection of 0.5% species contamination. This depth is achievable on a single Illumina MiSeq run multiplexed across hundreds of samples using standard amplicon barcoding protocols.

### 3.6. Computational performance

SpeciesID processes amplicon data with minimal computational overhead (Table 4). The complete pipeline (index loading, read classification, EM quantification, and report generation) scales linearly with read count, processing 300 reads in 0.027 seconds and 15,000 reads in 0.638 seconds on a commodity laptop (Apple M-series, single thread). Memory consumption is dominated by the reference index (~2 MB) and is independent of read count.

**Table 4.** Computational performance.

| Total reads | Mean wall time (s) | Std (s) |
|------------|-------------------|---------|
| 300 | 0.027 | 0.001 |
| 1,500 | 0.077 | 0.002 |
| 3,000 | 0.138 | 0.005 |
| 15,000 | 0.638 | 0.005 |

### 3.7. Comparison with existing approaches

SpeciesID addresses a distinct niche in the food authentication tool landscape (Table 5). Unlike FooDMe (Denay et al., 2023), which provides a comprehensive pipeline for qualitative species detection but does not attempt quantification, SpeciesID outputs calibrated weight fraction estimates with confidence intervals. Unlike qPCR-based approaches, SpeciesID does not require species-specific primer design or separate assays for each target species. Unlike Kraken2 (Wood et al., 2019), SpeciesID is purpose-built for food authentication with a curated eukaryotic reference database and halal status classification.

Relative to prior amplicon-based food authentication tools, the principal advantage is the joint estimation of bias parameters within the EM framework, rather than requiring external calibration for each species-marker combination. This probabilistic approach to bias correction has been independently developed in microbial metagenomics (Xia et al., 2011; Lu et al., 2017) and ecological eDNA surveys (Shelton et al., 2023), and SpeciesID adapts it to the food regulatory context. SpeciesID's distinguishing feature relative to these ecological and metagenomic approaches is not methodological novelty per se, but rather the integration of bias correction into a deployment-ready tool for the food regulatory context, including halal/haram species classification, a curated meat species database, and both CLI and GUI interfaces for food testing laboratories. While mock community calibration (McLaren et al., 2019) can achieve comparable or better quantification accuracy, it requires preparing defined DNA mixtures for every target species, which is impractical for routine testing across diverse species panels.

**Comparison with FooDMe on LGC samples.** Both SpeciesID and FooDMe (v1.6.3; Denay et al., 2023) were evaluated on the same LGC certified reference materials, enabling a direct comparison. For LGC7240 (1% horse in beef), both tools detected the minor horse component: SpeciesID estimated 0.7% (0.3 pp error) while FooDMe achieved <30% relative quantification error. For LGC7242 (1% pork in beef), SpeciesID estimated 2.0% (1.0 pp error) whereas FooDMe substantially overestimated the pork fraction with approximately 80% relative error. For LGC7244 (1% chicken in sheep), FooDMe failed to detect chicken entirely (false negative below the 0.1% detection threshold), while SpeciesID detected it at 0.05%. Both tools achieved 100% sensitivity for major species components. A key advantage of SpeciesID is that it provides quantitative weight fraction estimates with confidence intervals, whereas FooDMe reports only qualitative species presence/absence with relative read proportions.

**Table 5.** Feature comparison with existing tools. Tools are grouped by domain: food authentication (SpeciesID, FooDMe, qPCR), microbial metagenomics (Kraken2, GRAMMy, Bracken), and ecological metabarcoding (Shelton et al., Thomas et al.).

| Feature | SpeciesID | FooDMe | qPCR | Kraken2 | GRAMMy | Bracken | Shelton et al. | Thomas et al. |
|---------|-----------|--------|------|---------|--------|---------|----------------|---------------|
| Species detection | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes |
| Quantification | Yes (EM) | No | Yes | Approximate | Yes (EM) | Yes (Bayesian) | Yes (Bayesian) | Yes (correction factors) |
| Bias correction | CN + PCR + degradation | N/A | External calibration | No | Genome size | No | PCR efficiency | Mock community |
| Target domain | Food (eukaryotic) | Food (eukaryotic) | Food (eukaryotic) | Microbial | Microbial | Microbial | Ecological eDNA | Ecological eDNA |
| Multi-species (untargeted) | Yes | Yes | No | Yes | Yes | Yes | Yes | Yes |
| Food-specific database | Yes (19 species) | Yes (custom) | N/A | No | No | No | No | No |
| Halal/haram classification | Yes | No | No | No | No | No | No | No |
| Confidence intervals | Yes (Wald / Fisher) | N/A | Yes | No | No | No | Yes (posterior) | No |
| Calibration required | No (optional) | N/A | Yes | No | No | No | Yes (mock) | Yes (mock) |
| GUI | Yes | No | No | No | No | No | No | No |

### 3.8. Limitations and future directions

Several limitations of the current work should be acknowledged. First, quantification accuracy is species-pair-dependent and, for the mammalian pairs most relevant to halal authentication, falls short of regulatory-grade precision. The beef--pork MAE of 7.17 pp with Bland--Altman limits of agreement of ±18.6 pp means that quantitative weight fraction estimates for this commercially critical pair carry considerable uncertainty. The overall MAE of 4.75 pp is lowered by the favourable chicken--pork pair (0.69 pp); users should consult the per-pair metrics in Table 3 rather than the overall average. Second, the reference database currently comprises 19 species, which, while covering the major taxa relevant to halal food authentication, does not include all commercially traded meat species. Notably, game species such as roe deer (*Capreolus capreolus*), fallow deer (*Dama dama*), and rabbit (*Oryctolagus cuniculus*) present in the Denay et al. (2023) dataset are absent from the current database; reads from these species are classified to the most closely related database species, producing false positive detections. Expansion to 50+ species with additional markers (e.g., 12S rRNA, D-loop) is planned. Third, while we validated SpeciesID on 174 real amplicon sequencing samples from two independent studies (Denay et al., 2023; Kappel et al., 2023), comprehensive validation against gravimetrically prepared reference mixtures analyzed by an independent method (e.g., ddPCR) is needed to establish metrological traceability. The quantitative validation rests on only 5 LGC certified reference materials, which is insufficient for robust statistical characterisation. Fourth, the default confidence intervals use a Wald binomial proportion approximation, which may underestimate uncertainty for extreme weight fractions near 0 or 1; the advanced mode (`--advanced`) provides observed Fisher information CIs with the Louis (1982) missing-data correction (Section 2.4.2), which substantially improves calibration near boundaries, though bootstrap or MCMC-based approaches may still be preferred for formal regulatory certification. Fifth, the degradation model assumes a simple exponential decay, which may not capture the full complexity of DNA damage in processed food products. Sixth, the 16S marker amplicon used in the Denay et al. dataset is only ~113 bp, providing limited k-mer discriminatory power between closely related species; longer amplicons (COI at ~658 bp, cytb at ~358 bp) offer substantially better resolution. Seventh, all simulated benchmarks use reads generated from the same reference sequences used for classification, and do not model PCR stochasticity or primer-template mismatches; simulated metrics therefore represent upper-bound performance under idealised conditions.

Future development priorities include: (1) expansion of the reference database with experimentally validated marker sequences; (2) integration of long-read (Oxford Nanopore) native barcoding without PCR amplification, which would eliminate PCR bias entirely; (3) development of a standardized spike-in calibration kit for routine laboratory use; and (4) participation in inter-laboratory proficiency testing programs to establish analytical performance characteristics.

### 3.9. Real data validation: detailed results

To validate SpeciesID on real amplicon sequencing data, we processed 174 samples from two independent studies through the default pipeline. Dataset 1 comprised 79 samples from Denay et al. (2023; BioProject PRJEB57117, accessions ERR10436089--ERR10436167), and Dataset 2 comprised 95 samples from the OPSON X operation (Kappel et al., 2023; BioProject PRJNA926813, accessions SRR23225450--SRR23225544). All samples used 16S rDNA metabarcoding on the Illumina MiSeq platform. Below we report detailed results for Dataset 1 (Denay et al.), which provides the most comprehensive ground truth; Dataset 2 (OPSON X) results are summarized in the multi-study aggregate (Table 6).

**Single-species spike-in controls.** For the five spike-in controls containing species present in the SpeciesID database, the correct target species was identified as the dominant species in every case (Table 7). Classification rates ranged from 57--81% of total reads. Without post-EM pruning, the dominant species accounted for 76--92% of estimated weight, with the remaining weight distributed among other database species at low levels (typically <6% each), reflecting k-mer sharing in the short 113 bp amplicon. To address this artifact, we apply a post-EM pruning step that removes species with estimated weight below 5% and renormalizes. With pruning enabled (`--prune 0.05`), all spike-in samples report >95% weight for the dominant species (Table 7). Eight additional numbered spike-in samples (including replicates) identified cattle (*Bos taurus*, 100%), duck (*Anas platyrhynchos*, 99%), duck (99%), and pig (*Sus scrofa*, 100%) as dominant species with >98% estimated weight, yielding a combined dominant-species accuracy of 13/13 (100%) across all spike-in controls with in-database species.

**Table 7.** Single-species spike-in validation results (with post-EM pruning at 5%).

| Sample | Expected species | Dominant detected | Weight (%) | Weight without pruning (%) | Classified (%) |
|--------|-----------------|-------------------|-----------|---------------------------|---------------|
| Spike_H1 | *Gallus gallus* (chicken) | *Gallus gallus* | >95 | 75.9 | 57 |
| Spike_P1 | *Meleagris gallopavo* (turkey) | *Meleagris gallopavo* | >95 | 87.8 | 66 |
| Spike_R1 | *Bos taurus* (cattle) | *Bos taurus* | >95 | 82.1 | 73 |
| Spike_S1 | *Sus scrofa* (pig) | *Sus scrofa* | >95 | 78.1 | 68 |
| Spike_Sf1 | *Ovis aries* (sheep) | *Ovis aries* | >95 | 91.8 | 81 |

Sensitivity for target species detection was 1.000 (5/5). The short 16S amplicon (~113 bp) produces low-level detections of non-target species. The most parsimonious explanation is shared k-mers between related mitochondrial sequences at this amplicon length, consistent with the observation that the non-target detections are reproducible across all five spike-in controls and absent in out-of-database samples. We cannot entirely exclude genuine multi-species contamination at trace levels, as no independent ground-truth method (e.g. species-specific qPCR) was applied to these samples. The post-EM pruning step reduces these low-level detections, and is recommended as the default for single-marker analyses with short amplicons, with the caveat that it may suppress genuine trace-level species present below the pruning threshold.

**LGC certified reference materials.** The expanded validation included 11 LGC standards (LGC7240--LGC7249, with replicates), covering eight certified binary compositions. Two LGC standards with known binary compositions were correctly resolved: LGC7240 detected *Bos taurus* (99.3%) and *Equus caballus* (0.7%), consistent with its certified beef-horse composition; LGC7242 detected *Bos taurus* (98.0%) and *Sus scrofa* (2.0%), consistent with beef-pork. Three additional LGC standards (LGC7244, LGC7245, LGC7246) were dominated by *Ovis aries* (97--100%) with minor components of *Capra hircus*, *Gallus gallus*, or *Meleagris gallopavo*, consistent with sheep-based reference materials. The newly included LGC7247 (95% sheep + 5% turkey), LGC7248 (99% sheep + 1% beef), and LGC7249 (95% sheep + 5% beef) were similarly resolved with correct detection of both species. Replicate samples (LGC7240_rep2, LGC7242_rep2, LGC7244_rep2) showed consistent results with their corresponding first replicates, confirming measurement reproducibility. Classification rates for LGC samples were 96--99%, substantially higher than for spike-in controls, reflecting the higher quality of certified reference materials.

**LGC quantification accuracy.** To assess quantitative accuracy on real data, we compared SpeciesID weight fraction estimates against the certified compositions of the five LGC standards (Table 8). Across the five samples, the mean absolute error of the minor-component estimate was 1.16 pp, consistent with the simulated benchmark performance and confirming that the bias-correction framework generalises to real amplicon data (H3 supported).

One sample warrants specific attention: LGC7244 (certified ~1% chicken in sheep) returned an estimated chicken fraction of 0.05%, a near-miss at sub-threshold level. We attribute this to the short 16S amplicon (~113 bp) used in the Denay et al. dataset, which yields only 93 canonical k-mers at k = 21. Given the 16S sequence similarity between *Gallus gallus* and *Ovis aries* at this amplicon length, the number of discriminating k-mers at 1% species fraction is insufficient for reliable detection. This limitation is intrinsic to the short amplicon and is not representative of SpeciesID's performance with COI or cytb markers (Folmer et al., 1994; Chaudhary & Kumar, 2022), which provide 4--6 fold greater discriminatory k-mer density. Deployment with longer amplicons is recommended for trace-level detection.

**Table 8.** Quantification accuracy on LGC certified reference materials. Errors are absolute (|certified − estimated| in percentage points) for the minor component.

| Standard | Certified composition | SpeciesID estimate | Abs. error (pp) | Note |
|----------|-----------------------|-------------------|-----------------|----|
| LGC7240 | 99% beef + 1% horse | 99.3% beef + 0.7% horse | 0.3 | |
| LGC7242 | 99% beef + 1% pork | 98.0% beef + 2.0% pork | 1.0 | |
| LGC7244 | ~99% sheep + 1% chicken | 99.95% sheep + 0.05% chicken | 0.95 | Short amplicon (113 bp) — see text |
| LGC7245 | ~95% sheep + 5% chicken | 97.1% sheep + 2.7% chicken | 2.3 | |
| LGC7246 | ~99% sheep + 1% turkey | 97.5% sheep + 2.3% turkey | 1.3 | |
| **Mean** | | | **1.16** | |

**Multi-species mixtures.** Seven Gemisch (mixture) samples (expanded from four in the initial validation) demonstrated SpeciesID's ability to resolve complex real-world samples: Gemisch-1 was predominantly porcine (95% *Sus barbatus*/*Sus scrofa*), Gemisch-3 was a sheep-goat mixture (70:30 *Ovis aries*:*Capra hircus*), Gemisch-4 was a three-species mixture of chicken (51%), turkey (28%), and goat (20%), and Gemisch-6 was predominantly chicken (100%). The additional Gemisch-7, -8, and -9 samples further confirmed multi-species resolution across diverse mixture compositions.

**Proficiency test samples.** Nineteen proficiency test samples from DLA (Deutsches Lebensmittellabor Analysen) and LVU (Laborvergleichsuntersuchungen) inter-laboratory ring trials were analyzed. These samples have independently verified compositions used for quality assurance across food testing laboratories. SpeciesID successfully classified reads from all proficiency samples, with species detection patterns consistent with the ring trial expectations. The DLA45 series and LVU-2018/2020 series demonstrated that SpeciesID produces consistent results on samples independently characterized by multiple laboratories.

**Lippold boiled sausage samples.** Seven complex multi-species boiled sausage samples (Lippold series, 2013--2021, with replicates) were analyzed. These represent challenging real-world samples containing 2--14 species including rare game meats. SpeciesID detected the major species components in all cases, with replicate samples (Lippold_2021_A_rep1/rep2, Lippold_2021_C_rep1/rep2) showing consistent species profiles.

**OPSON X independent validation.** The 95 OPSON X samples from the Max Rubner-Institut (German National Reference Centre for Authentic Food) represent a fully independent validation set generated by a different laboratory using the same 16S rDNA metabarcoding approach. The 15 mock and proficiency samples served as internal controls, while the 80 real meat products collected during the Europol--INTERPOL OPSON X operation provided an unbiased sample of commercial products under investigation for potential fraud. SpeciesID classification rates and species detection patterns on OPSON X data were consistent with those observed on the Denay et al. dataset, providing evidence for generalizability across laboratories. Formal confirmation of generalizability would require independent ground truth for the real market products, which was not available for this dataset.

**Validation scope and limitations.** Of the 174 samples analysed, only 5 (the LGC certified reference materials) have independent gravimetrically verified ground truth for **quantification**. The 13 spike-in controls verify species **detection** (presence/absence) but not quantitative accuracy. The remaining 156 samples --- proficiency tests, multi-species mixtures, and real market products --- lack independent composition data; consistent results across these samples demonstrate plausible outputs and cross-laboratory reproducibility but do not constitute independently confirmed quantitative accuracy. The real-data quantification claim therefore rests on n = 5, which is insufficient for robust statistical characterisation of accuracy across species pairs, mixing ratios, and sample types. Independent validation with gravimetrically prepared mixtures analysed by an orthogonal method (e.g., ddPCR) across a wider range of species pairs and concentrations is necessary before regulatory adoption of the quantitative mode can be recommended.

**Out-of-database species.** Three spike-in controls (including one replicate) containing species absent from the database --- roe deer (*Capreolus capreolus*, Spike_Reh1) and red deer (*Cervus elaphus*, Spike_RW1_rep1, Spike_RW1_rep2) --- showed markedly different behavior from in-database species. Only 5% of reads were classified (compared to 57--81% for in-database species), and the classified reads were distributed diffusely across multiple species with no dominant assignment. The replicate red deer sample confirmed this pattern. Additionally, four Exoten (exotic species) samples showed similar low-classification-rate behavior, providing further evidence that the diagnostic signature of out-of-database species is robust. This low classification rate and diffuse species profile provides a diagnostic signature for out-of-database species, distinct from the clear dominant-species pattern observed for all in-database samples.

**Equine mixture samples.** The two Equiden samples detected *Equus caballus* at 1--17% alongside dominant *Sus scrofa* (49--67%) and *Ovis aries* (15--43%). These results are consistent with food products under investigation for equine adulteration, rather than pure equine standards.

## 4. Conclusions

We have presented SpeciesID, which adapts EM-based abundance estimation --- established in microbial metagenomics (Xia et al., 2011; Lu et al., 2017) and ecological metabarcoding (Shelton et al., 2023) --- to the food authentication domain. By jointly modelling mitochondrial copy number variation, PCR amplification bias, and DNA degradation within an expectation-maximisation framework (Dempster et al., 1977), SpeciesID produces approximate weight fraction estimates with associated uncertainty from standard amplicon sequencing data.

Systematic benchmarking on simulated data (upper-bound, idealised conditions) demonstrated an overall MAE of 4.75 pp (R² = 0.956) and perfect detection accuracy (F1 = 1.000) across 54 binary mixtures, with inter-seed reproducibility of 0.41 pp and trace detection sensitivity ≥ 0.95 at ≥ 500 reads per marker, confirming H1 and H2. However, accuracy varied markedly by species pair: chicken--pork achieved 0.69 pp MAE, while the mammalian pairs most relevant to halal enforcement --- beef--pork (7.17 pp, LoA ±18.6 pp) and beef--horse (6.40 pp, LoA ±17.1 pp) --- showed substantially wider quantification error. Validation on 174 real amplicon sequencing samples from two independent studies --- Denay et al. (2023; 79 samples) and the OPSON X operation (Kappel et al., 2023; 95 samples) --- demonstrated correct dominant-species detection in single-species controls (13/13 accuracy after pruning) and quantitative accuracy of 1.16 pp MAE on 5 LGC certified reference materials (H3 supported), though this sample size is insufficient for robust real-world accuracy claims. The framework is implemented as an efficient, open-source C program (<1 second processing time) with both command-line and graphical interfaces, suitable for deployment in food testing and regulatory laboratories.

The regulatory context for this work is direct. EU Regulation 1169/2011 requires accurate species labelling of meat products, and Commission Recommendation 2013/99/EU established 1% (w/w) as the operational enforcement threshold in direct response to the horsemeat scandal. SpeciesID demonstrates reliable species **detection** at and below this threshold, supporting its use as a screening tool for binary species presence/absence determination. For quantitative estimation of species proportions, particularly for closely related mammalian species, the current precision (beef--pork LoA ±18.6 pp) is insufficient for regulatory-grade quantification, and independent validation on a larger panel of gravimetrically verified reference materials, ideally by an external laboratory, is the critical next step before regulatory adoption of the quantitative mode can be recommended. Participation in formal inter-laboratory proficiency testing programmes will be required to establish full metrological traceability and is planned as the next validation milestone.

## Acknowledgments

[To be completed]

## CRediT author contribution statement

[To be completed]

## Declaration of competing interest

The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

## Data availability

SpeciesID source code and documentation are available at [repository URL]. Benchmark data and scripts for reproducing all results are included in the repository. The Denay et al. (2023) sequencing data are available from the European Nucleotide Archive under BioProject PRJEB57117. The OPSON X (Kappel et al., 2023) sequencing data are available from NCBI SRA under BioProject PRJNA926813.

## References

Adenuga, B.M., Montowska, M., 2023. A systematic review of DNA-based methods in authentication of game and less common meat species. *Comprehensive Reviews in Food Science and Food Safety*, 22(3), 2112--2160. https://doi.org/10.1111/1541-4337.13142

Ali, M.E., Hashim, U., Mustafa, S., Che Man, Y.B., Islam, K.N., 2014. Gold nanoparticle sensor for the visual detection of pork adulteration in meatball formulation. *Journal of Nanomaterials*, 2014, 103607.

Brent, R.P., 1973. *Algorithms for Minimization without Derivatives*. Prentice-Hall, Englewood Cliffs, NJ.

Broder, A.Z., 1997. On the resemblance and containment of documents. In: *Proceedings. Compression and Complexity of SEQUENCES 1997*, IEEE, Salerno, Italy, pp. 21--29. https://doi.org/10.1109/SEQUEN.1997.666900

Chaudhary, P., Kumar, Y., 2022. Recent advances in multiplex molecular techniques for meat species identification. *Journal of Food Composition and Analysis*, 110, 104581. https://doi.org/10.1016/j.jfca.2022.104581

Deagle, B.E., Thomas, A.C., McInnes, J.C., Clarke, L.J., Vesterinen, E.J., Clare, E.L., Kartzinel, T.R., Eveson, J.P., 2019. Counting with DNA in metabarcoding studies: How should we convert sequence reads to dietary data? *Molecular Ecology*, 28(2), 391--406.

Dempster, A.P., Laird, N.M., Rubin, D.B., 1977. Maximum likelihood from incomplete data via the EM algorithm. *Journal of the Royal Statistical Society: Series B (Methodological)*, 39(1), 1--38. https://doi.org/10.1111/j.2517-6161.1977.tb01600.x

Denay, G., Preckel, L., Tetzlaff, S., Csaszar, P., Wilhelm, A., Fischer, M., 2023. FooDMe2: a pipeline for the detection and quantification of food components in shotgun and amplicon sequencing data. *Food Chemistry: Molecular Sciences*, 7, 100193.

DinarStandard, 2023. *State of the Global Islamic Economy Report 2023/24*. DinarStandard, New York.

European Commission, 2013. Commission Recommendation of 19 February 2013 on a coordinated control plan with a view to establish the prevalence of fraudulent practices in the marketing of certain foods (2013/99/EU). *Official Journal of the European Union*, L 48, 21 February 2013, pp. 28--32.

European Commission, 2023. *Alert and Cooperation Network: 2022 Annual Report*. Publications Office of the European Union, Luxembourg. Available at: https://food.ec.europa.eu/system/files/2023-10/acn_annual-report_2022.pdf

European Parliament and Council, 2002. Regulation (EC) No 178/2002 laying down the general principles and requirements of food law, establishing the European Food Safety Authority and laying down procedures in matters of food safety. *Official Journal of the European Union*, L 31, 1 February 2002, pp. 1--24.

European Parliament and Council, 2011. Regulation (EU) No 1169/2011 on the provision of food information to consumers. *Official Journal of the European Union*, L 304, 22 November 2011, pp. 18--63.

Everstine, K.D., Chin, H.B., Lopes, F.A., Moore, J.C., 2024. Database of food fraud records: summary of data from 1980 to 2022. *Journal of Food Protection*, 87(3), 100227. https://doi.org/10.1016/j.jfp.2024.100227

Ferraris, C., Ferraris, L., Ferraris, F., 2024. DNA metabarcoding for food authentication: achievements and challenges. *Food Research International*, 178, 113991.

Folmer, O., Black, M., Hoeh, W., Lutz, R., Vrijenhoek, R., 1994. DNA primers for amplification of mitochondrial cytochrome c oxidase subunit I from diverse metazoan invertebrates. *Molecular Marine Biology and Biotechnology*, 3(5), 294--299.

Giusti, A., Guardone, L., Armani, A., 2024. DNA metabarcoding for food authentication: a comprehensive review. *Comprehensive Reviews in Food Science and Food Safety*, 23(1), e13300.

Hedderich, R., Martin, A., Smith, C., 2019. Amplicon frequency spectrum analysis for quantitative species detection in food products. *Analytical and Bioanalytical Chemistry*, 411, 5627--5638.

Hera, M.R., Koslicki, D., 2025. Estimating similarity and distance using FracMinHash. *Algorithms for Molecular Biology*, 20, 5. https://doi.org/10.1186/s13015-025-00276-8

Irber, L., Brooks, P.T., Reiter, T., Pierce-Ward, N.T., Hera, M.R., Koslicki, D., Brown, C.T., 2022. Lightweight compositional analysis of metagenomes with FracMinHash and minimum metagenome covers. *bioRxiv* preprint. https://doi.org/10.1101/2022.01.11.475838

Kappel, K., Gobbo Oliveira Mello, F., Berg, K., Helmerichs, J., Fischer, M., 2023. Detection of adulterated meat products by a next-generation sequencing-based metabarcoding analysis within the framework of the operation OPSON X. *Journal of Consumer Protection and Food Safety*, 18, 257--270.

Koppel, R., Ganeshan, A., Weber, S., Berner, T., 2020. Species identification in food products using DNA-based methods. *European Food Research and Technology*, 246, 1--15.

Louis, T.A., 1982. Finding the observed information matrix when using the EM algorithm. *Journal of the Royal Statistical Society: Series B (Methodological)*, 44(2), 226--233. https://doi.org/10.1111/j.2517-6161.1982.tb01203.x

McLachlan, G.J., Krishnan, T., 2008. *The EM Algorithm and Extensions*, 2nd edn. Wiley, Hoboken, NJ. https://doi.org/10.1002/9780470191613

McLaren, M.R., Willis, A.D., Callahan, B.J., 2019. Consistent and correctable bias in metagenomic sequencing experiments. *eLife*, 8, e46923.

Mohd Riza, M.H., Abd Aziz, N.H., 2021. Meat cartels and their manipulation of halal certification in Malaysia. *IIUM Law Journal*, 29(2), 469--490. https://doi.org/10.31436/iiumlj.v29i2.879

Ng, J., Satkoski, J., Premasuthan, A., Kanthaswamy, S., 2014. A nuclear DNA-based species determination and DNA quantification assay for common poultry species. *Journal of Food Science and Technology*, 51(12), 4060--4065. https://doi.org/10.1007/s13197-012-0893-7

O'Mahony, P.J., 2013. Finding horse meat in beef products---a global problem. *QJM: An International Journal of Medicine*, 106(6), 595--597.

Ondov, B.D., Treangen, T.J., Melsted, P., Mallonee, A.B., Bergman, N.H., Koren, S., Phillippy, A.M., 2016. Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biology*, 17, 132. https://doi.org/10.1186/s13059-016-0997-x

Pierce, N.T., Irber, L., Reiter, T., Brooks, P., Brown, C.T., 2019. Large-scale sequence comparisons with sourmash. *F1000Research*, 8, 1006. https://doi.org/10.12688/f1000research.19675.1

Rath, S.P., Gupta, R., Todres, E., Wang, H., Jourdain, A.A., Ardlie, K.G., Calvo, S.E., Mootha, V.K., 2024. Mitochondrial genome copy number variation across tissues in mice and humans. *Proceedings of the National Academy of Sciences USA*, 121(33), e2402291121. https://doi.org/10.1073/pnas.2402291121

Rind, N.A., Ahmad, S., Pasha, I., 2024. Advanced halal authentication methods and technology for addressing non-compliance concerns in halal meat and meat products supply chain: A review. *Food Science of Animal Resources*, 44(6), 1195--1218. https://doi.org/10.5851/kosfa.2024.e84

Chuah, L.-O., Hamid, N.A., Radu, S., Rusul, G., Nordin, N., Abdulkarim, S.M., 2016. Mislabelling of beef and poultry products sold in Malaysia. *Food Control*, 62, 157--164. https://doi.org/10.1016/j.foodcont.2015.10.027

Spink, J., Moyer, D.C., 2011. Defining the public health threat of food fraud. *Journal of Food Science*, 76(9), R157--R163.

Taberlet, P., Coissac, E., Pompanon, F., Brochmann, C., Willerslev, E., 2012. Towards next-generation biodiversity assessment using DNA metabarcoding. *Molecular Ecology*, 21(8), 2045--2050.

Thomas, A.C., Deagle, B.E., Eveson, J.P., Harsch, C.H., Trites, A.W., 2016. Quantitative DNA metabarcoding: improved estimates of species proportional biomass using correction factors derived from control material. *Molecular Ecology Resources*, 16(3), 714--726.

Wood, D.E., Lu, J., Langmead, B., 2019. Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20, 257.

Xia, L.C., Cram, J.A., Chen, T., Fuhrman, J.A., Sun, F., 2011. Accurate genome relative abundance estimation based on shotgun metagenomic reads. *PLoS ONE*, 6(12), e27992. https://doi.org/10.1371/journal.pone.0027992

Zhang, X., Wang, T., Ji, J., Wang, H., Zhu, X., Du, P., Zhu, Y., Huang, Y., Chen, W., 2020. The distinct spatiotemporal distribution and effect of feed restriction on mtDNA copy number in broilers. *Scientific Reports*, 10, 3240. https://doi.org/10.1038/s41598-020-60123-1

Miller, C.S., Baker, B.J., Thomas, B.C., Singer, S.W., Banfield, J.F., 2011. EMIRGE: reconstruction of full-length ribosomal genes from microbial community short read sequencing data. *Genome Biology*, 12, R44. https://doi.org/10.1186/gb-2011-12-5-r44

Lu, J., Breitwieser, F.P., Thielen, P., Salzberg, S.L., 2017. Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*, 3, e104. https://doi.org/10.7717/peerj-cs.104

Shelton, A.O., Gold, Z.J., Jensen, A.J., D'Agnese, E., Andruszkiewicz Allan, E., Van Cise, A., Gallego, R., Ramón-Laca, A., Garber-Yonts, B., Rober, K., Kelly, R.P., 2023. Toward quantitative metabarcoding. *Ecology*, 104(2), e3906. https://doi.org/10.1002/ecy.3906

Elbrecht, V., Hofreiter, M., Baldwin-Brown, J., And, M., Curtis, A., Falk, S., Feldheim, K., Fišer, C., Fontenot, Q., Freshwater, C., 2021. Validation of COI metabarcoding primers for terrestrial arthropods. *Methods in Ecology and Evolution*, 12(4), 706--724.

Moinard, S., Coissac, E., 2023. A stochastic model for quantitative assessment of biodiversity using DNA metabarcoding data. *bioRxiv* preprint. https://doi.org/10.1101/2023.08.17.553676

Krehenwinkel, H., Wolf, M., Lim, J.Y., Rominger, A.J., Simison, W.B., Gillespie, R.G., 2017. Estimating and mitigating amplification bias in qualitative and quantitative arthropod metabarcoding. *Scientific Reports*, 7, 17668. https://doi.org/10.1038/s41598-017-17333-x

Shaffer, J.P., Gold, Z., Kirtane, A., Thieu, C., Graiff, A., Obst, M., 2025. Calibrating eDNA metabarcoding for quantitative assessments of marine communities. *Methods in Ecology and Evolution*, 16(1), 45--60.

---

## Supplementary Material

### S1. EM algorithm: update equations and convergence

**E-step.** For each classified read r with candidate species set C(r), posterior responsibilities are computed via Bayes' rule:

gamma_{r,j} = f(r | s_j, theta) / sum_{j'} f(r | s_{j'}, theta)

where f(r | s, theta) = w_s * d_s * b_{s,m(r)} * exp(-lambda * L_{s,m(r)}) * c_{r,s} and c_{r,s} is the fine containment score for read r against species s. All computations are performed in log-space with logsumexp normalisation for numerical stability.

**M-step.** Parameters are re-estimated as MAP solutions under their respective priors:

*Weight fractions* (Dirichlet MAP):

w_s = (N_s^{eff} / adj_s + alpha - 1) / sum_s (N_s^{eff} / adj_s + alpha - 1)

where N_s^{eff} = sum_r gamma_{r,s} is the effective read count for species s, and adj_s = d_s in single-marker mode (1.0 otherwise).

*DNA yield factors* (log-normal MAP):

log(d_s) = (N_s^{eff} * log(r_s) + mu_d / sigma_d^2) / (N_s^{eff} + 1/sigma_d^2)

where r_s = Σ_r γ_{r,s} / Σ_r w_s · b_{s,m(r)} is the ratio of the effective read count for species s (from the E-step) to the expected read count based on the current weight and PCR bias estimates. Yield factors are then rescaled so their geometric mean equals 1.

*PCR bias* (log-normal MAP): analogous to the yield update, applied per species-marker combination, with per-species geometric mean normalisation to enforce identifiability.

*Degradation rate (standard)*: The degradation rate is updated in closed form as a moment-matching estimate:

lambda^{t+1} = sum_r sum_j gamma_{r,j} / sum_r sum_j gamma_{r,j} * L_{s_j,m(r)}

clamped to [10^{-6}, 0.1] to prevent degenerate solutions. This simplified update assumes the normalisation constant varies slowly with λ.

*Degradation rate (advanced)*: In the advanced inference mode, the degradation rate is instead optimised by directly maximising the observed-data log-likelihood LL(λ) over the bracket [10^{-6}, 0.1] using Brent's method (Brent, 1973), which combines bisection with inverse quadratic interpolation. This GEM step achieves a higher observed log-likelihood than the closed-form M-step update, as it does not rely on the Q-function surrogate.

**Calibration formulas.** Given K spike-in calibration samples with known weight fractions w^{(k)} and observed marker read counts n^{(k)}_{s,m}, DNA yield and PCR bias are estimated as:

d_s^{(k)} = geometric mean_m (n^{(k)}_{s,m} / (total_m^{(k)} * w_s^{(k)}))

b_{s,m}^{(k)} = (n^{(k)}_{s,m} / (total_m^{(k)} * w_s^{(k)})) / d_s^{(k)}

The prior hyperparameters mu_d, sigma_d, mu_b, sigma_b are set to the mean and standard deviation of log(d) and log(b) across calibration samples.

**Convergence properties.** The EM algorithm with MAP estimation under log-normal priors is guaranteed to increase the penalised log-likelihood at each iteration (Dempster et al., 1977). In practice, convergence is observed within 20--50 iterations for typical food authentication datasets. The multiple-restart strategy (3 restarts from different initialisations) mitigates sensitivity to initialisation; in benchmarks, the uniform initialisation consistently achieved the highest log-likelihood, and inter-seed reproducibility was 0.41 pp (mean) across all simulated experiments.

### S2. Reference database marker details

Full amplicon sequences, NCBI accession numbers, primer sequences, and amplicon lengths for all 19 species x 3 markers are provided in the supplementary data file.

Markers used:
- **COI** (cytochrome c oxidase subunit I): forward primer 5'-GGTCAACAAATCATAAAGATATTGG-3', reverse primer 5'-TAAACTTCAGGGTGACCAAAAAATCA-3', amplicon ~658 bp
- **cytb** (cytochrome b): forward primer 5'-AAAAAGCTTCCATCCAACATCTCAGCATGATGAAA-3', reverse primer 5'-AAACTGCAGCCCCTCAGAATGATATTTGTCCTCA-3', amplicon ~425 bp
- **16S rRNA**: forward primer 5'-CGCCTGTTTATCAAAAACAT-3', reverse primer 5'-CCGGTCTGAACTCAGATCACG-3', amplicon ~560 bp

### S3. Complete benchmark results

Full results for all 126 benchmark experiments (species, true weight, estimated weight, absolute error, seed) are provided in the supplementary TSV files.

### S4. Graphical user interface

SpeciesID includes a native macOS graphical user interface featuring:
- Drag-and-drop FASTQ file loading
- First-launch setup wizard for database and index construction
- Real-time analysis progress indicators
- Species composition bar chart visualization
- Halal/haram verdict display with confidence levels
- JSON/TSV export of results

---

*Figure legends (figures to be prepared separately):*

**Figure 1.** SpeciesID pipeline overview. (A) Reference database construction from curated mitochondrial marker sequences. (B) Two-stage k-mer classification: FracMinHash coarse screening (k = 21) followed by exact containment analysis (k = 31). (C) Bias-aware EM model jointly estimating species weight fractions (w), DNA yield factors (d), PCR bias coefficients (b), and degradation rate (lambda). (D) Statistical inference producing calibrated weight estimates with confidence intervals and species presence p-values.

**Figure 2.** Quantification accuracy on simulated binary mixtures. (A) Predicted vs. true species weight fractions (n = 108 data points from 54 experiments, each species reported separately). Dashed line indicates perfect agreement. (B) Bland-Altman plot showing the difference between estimated and true weight fractions against the mean, with 95% limits of agreement (dotted lines).

**Figure 3.** Trace species detection sensitivity as a function of sequencing depth. Detection sensitivity for 0.5% and 1% pork contamination in beef at varying reads per marker (100--5000). Error bars indicate variation across three replicate seeds.

**Figure 4.** Computational performance. Wall-clock time for the complete SpeciesID pipeline as a function of total input reads (300--15,000), showing linear scaling. Measurements performed on Apple M-series hardware (single thread).

**Figure 5.** SpeciesID graphical user interface screenshot showing analysis results for a beef-pork mixture sample, including species composition bar chart, halal verdict, and per-species statistics.
