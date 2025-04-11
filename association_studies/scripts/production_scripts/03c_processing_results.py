
###process prredictive power and make manhatann plots


#we need to check inflation? maybe just bonferroni and show no snp is signifncat, so PRSs are useful

#check overlap across iterations for the samples with the highest PRS values in the training or test set?


##WHEN INTERPRETING THE RESULTS OF THE PRS, LOOK SLIDES FROM DOUG
    #https://dougspeed.com/short/
###WHEN DONE, YOU CAN CHECK THE SECTION OF INTERPRETATION FROM ORRELLY
    #Interpretation and presentation of results


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


    ##INFLATION FACTOR?

    #check manhatan plots
        #there is a strange gap in VO2 max for one of the first chromosomes in the prelim results



#FOR BAT ANALYSES
    #this would an additional step in this project that would be outside of the paper
    #take the 1000kb gene windows for all coding genes, liftover to hg38. If the USCS tool accepts genomic ranges, just use them as input, if not, split in two datasets the start and the end of the gene windows
    #for each phenotype (VO2, beep....), calculate the average (better than median because want influence of outliers within gene like in iHS, if a SNPs is veery important in a gene that should influence the info about the whole gene) p-value for the association of SNPs inside each gene
    #then, calculate 1000 random sets of genes, within each set, calculate the median association of all genes inside the set and compare with the BAT set to obtain a distribution and empirical p-value (is association lower in BAT? LOOF BAT PAPER SCRIPTS FOR THIS). Here I want median because i do not want a gene outliser change things, I want the overall impact of BAT genes in general, not just a few genes.
