

#maybe do general PRS plots with the whole dataset using elastic?
    #Association and goodness-of-fit metrics
        #https://www.nature.com/articles/s41596-020-0353-1


"""
Notes about Genomic control and confounding effects from Dr. Speed's slides

- Confounding due to population structure in GWAS
    Compared with observational epidemiology, GWAS have few opportunities for confounding bias. The main problem is population structure, which refers to 
        - mating patterns within a pop -> subpops (more relatedness within than between)
        - allele frequency differences across subpops 
        - environmental exposures may also vary across subpops. 

    Phenotypes can also vary across subpops, because 
        - the causal alleles vary in frequency and/or
        - they vary with environmental factors correlated with pop structure, and/or 
        - ascertainment bias: recruitment of phenotypic groups differs across subpops 

    These can lead to significant associations not with CVs but with SNPs whose allele frequencies correlate with trait across subpops

- Adjustments for pop structure confounding: genomic control (GC) + PCA
    GC was an early approach to adjusting for confounding, based on the idea that pop structure can lead to many significant SNPs genome-wide 
        - all association test stats (with χ21 null distribution) are divided by the ratio of empirical to null medians (called a genomic inflation factor, GIF) provided GIF > 1 
        - assumes sparsity: true causals are rare, so there are few non-null test stats, so median test statistic is close to null value. 
    
    However, the omnigenic nature of many complex traits means that the assumption is false and GC is overly conservative. 
    
    The first few eigenvectors (or principal components) of XX T often reflect pop structure   
        - Included as covariates in GWAS regression models, they can absorb pop structure effects on the trait. 
    
    Now Mixed Model Association Analysis (MMAA) is the preferred approach to adjusting for pop structure effects in tests of association (see slide 40 from Dr. Speed's slides; teaching_slides.pdf)

- Given we have performed a fine analysis of population structure, removing many PCA outliers (we only have 1 ancestry) and then calculating PCAs on very clean data and considering these PCAs in the models, we have covered this point, so no need for genomic control.
"""

"""
Quantile plot of transformed phenotype against PRS quantiles

Quantile plots corresponding to the effect of a PRS on a normally distributed target trait should reflect the S-shape of the probit function (Fig. 5a). This is because the trait values are more spread out between quantiles at the tails of a normal distribution. Thus, plotting quantiles of PRS versus (absolute) effect on trait shows increasingly larger jumps up/down the y-axis from the median to the extreme upper/lower quantiles. When unequal strata are plotted, with the smallest strata at the tails, then this effect appears stronger. When the target outcome is disease status and prevalence or OR are plotted on the y-axis, then the shape is expected to be different: here, the shape is asymmetrical, showing a marked inflection at the upper end (Fig. 5b), since cases are enriched at the upper end only. Thus, inflections of risk at the tails of the PRS distribution82,83 should be interpreted according to these statistical expectations and not as interesting in themselves.

https://www.nature.com/articles/s41596-020-0353-1
"""


def manhattan_plot(covariate_dataset, response_variable):
        

    dict_change_names={
        "family_id": "FID",
        "AGRF code": "IID",
        "Age": "age",
        "Week 1 Body Mass": "week_1_weight",
        "Week 1 Beep Test": "week_1_beep",
        "Week 1 Distance (m)": "week_1_distance",
        "Week 1 Pred VO2max": "week_1_vo2"
    }

    print_text("specify the covariates", header=3)
    selected_covariates = pheno_subset.columns[~pheno_subset.columns.isin(["family_id", "AGRF code", response_variable])]
    print(selected_covariates)
    #check we have correct covariates
    total_list_covariates = ["Age", "sex_code", "Week 1 Body Mass", "Week 1 Beep Test", "Week 1 Distance (m)", "Week 1 Pred VO2max"] + [f"PCA{i}" for i in range(1,21)]
    if(sum([1 for cov in selected_covariates if cov not in total_list_covariates])!=0):
        raise ValueError("ERROR: FALSE! WE HAVE COVARIATES THAT ARE NOT IN THE TOTAL LIST")


    pheno_subset_transform = pd.read_csv( \
        "./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_transform.tsv", \
        sep="\t", \
        header=True, \
        index=False, \
        na_rep="NA" \
    )

    #save the response variable after changing the names
    pheno_subset_transform[["family_id", "AGRF code", response_variable]].rename(columns=dict_change_names).to_csv( \
        "./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_transform_response.tsv", \
        sep="\t", \
        header=True, \
        index=False, \
        na_rep="NA" \
    )

    #save the covariates that are factors
    if ("sex_code" in selected_covariates):

        #save sex_code
        pheno_subset_transform[["family_id", "AGRF code", "sex_code"]].rename(columns=dict_change_names).to_csv( \
            "./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_transform_covars_factors.tsv", \
            sep="\t", \
            header=True, \
            index=False, \
            na_rep="NA" \
        )

    #save the covariates that are continuous
    selected_covariates_cont = [cov for cov in selected_covariates if cov != "sex_code"]
    pheno_subset_transform[["family_id", "AGRF code"] + selected_covariates_cont].rename(columns=dict_change_names).to_csv( \
        "./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_transform_covars_cont.tsv", \
        sep="\t", \
        header=True, \
        index=False, \
        na_rep="NA" \
    )

    #create a file with the samples of the selected set
    pheno_subset_transform[["family_id", "AGRF code"]].to_csv( \
        "./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_samples_in.tsv", \
        sep="\t", \
        header=False, \
        index=False, \
        na_rep="NA" \
    )

    #use that file to select the corresponding sample from the plink fileset
    run_bash(" \
        plink \
            --bfile ./data/plink_filesets/" + covariate_dataset + "/" + response_variable + "_filesets/" + response_variable + "_subset_missing_clean_maf_hwe_sample_snp_missing \
            --keep ./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_samples_in.tsv \
            --make-bed \
            --out ./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_plink_fileset \
    ")
        #--keep accepts a space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, and removes all unlisted samples from the current analysis. --remove does the same for all listed samples.


    run_bash(" \
        ldak6.1.linux \
            --linear ./results/final_results/" + covariate_dataset + "/" + response_variable + "/full_dataset/" + response_variable + "_full_dataset_linear \
            --bfile ./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_plink_fileset \
            --pheno ./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_transform_response.tsv \
            --covar ./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_transform_covars_cont.tsv \
            --factors ./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_transform_covars_factors.tsv \
            --permute NO \
    ")

    #--linear
        #https://dougspeed.com/single-predictor-analysis/


    #clumping?

    print("load assoc results to pandas")
    assoc_results = pd.read_csv( \
        "./results/final_results/" + covariate_dataset + "/" + response_variable + "/full_dataset/" + response_variable + "_full_dataset_linear.assoc", \
        sep="\t", \
        header=0, \
        low_memory=False)
    print(assoc_results)
        #CHECK COLUMNS
        #compare wald P with the .pvalue file
    

    print("check we have the correct dtypes for dash_bio.ManhattanPlot (see code below)")
    print(assoc_results["Chromosome"].dtype == "int64")
    print(assoc_results["Basepair"].dtype == "int64")
    print(assoc_results["Wald_P"].dtype == "float64")
    print(assoc_results["Predictor"].dtype == "O")

    dict_titles = { \
        "distance_change": "Change in distance (m)", \
        "beep_test": "Change in beep test", \
        "vo2max": "Change in predicted VO2max", \
        "weight": "Change in body mass", \
    }

    # import libraries
    from scipy.stats import uniform
    from scipy.stats import randint
    import matplotlib.pyplot as plt


    # -log_10(pvalue)
    assoc_results['minuslog10pvalue'] = -np.log10(assoc_results["Wald_P"])
    assoc_results["Chromosome"] = assoc_results["Chromosome"].astype('category')
    assoc_results = assoc_results.sort_values('Chromosome')

    # How to plot gene vs. -log10(pvalue) and colour it by chromosome?
    assoc_results['ind'] = range(len(assoc_results))
    df_grouped = assoc_results.groupby(('Chromosome'))

    # manhattan plot
    fig = plt.figure(figsize=(22, 8)) # Set the figure size
    ax = fig.add_subplot(111)
    colors = ['darkred','darkgreen','darkblue', 'gold']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='ind', y='minuslog10pvalue',color=colors[num % len(colors)], ax=ax)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)

    genome_wide_threshold = -np.log10(0.05/assoc_results.shape[0])
    ax.axhline(y=genome_wide_threshold, color='red', linestyle='--', linewidth=1.5, label=f'Genome-wide threshold ({genome_wide_threshold:.2f})')
    
    # set axis limits
    ax.set_xlim([0, len(assoc_results)])
    if (assoc_results['minuslog10pvalue'].max() > genome_wide_threshold):
        ax.set_ylim([0, assoc_results['minuslog10pvalue'].max() + 2])
    else:
        ax.set_ylim([0, genome_wide_threshold + 2])

    # x axis label
    ax.set_xlabel('Chromosome')

    # Save the plot as a static image
    plt.savefig("./results/final_results/" + covariate_dataset + "/" + response_variable + "/full_dataset/" + response_variable + ".png", dpi=300)

    #code from:
        #https://python-graph-gallery.com/manhattan-plot-with-matplotlib/



    #check manhatan plots
        #there is a strange gap in VO2 max for one of the first chromosomes in the prelim results



