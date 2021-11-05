
PROP = ../properties/supplement.properties
CACHE = ../supplement_cache
NORM = --normalization RobustLinearNormalization
PARALLEL = --multiprocessing
#PARALLEL = --lsf-directory ../computation --memory 4
PARALLEL_LONG = --lsf-directory ../computation --memory 4
PARALLEL_LONG8 = --lsf-directory ../computation --memory 8
#PARALLEL_LONG = 

all: figures tables

.SECONDARY:

#
# Rules for subsampling cells for building models. We sample
# $(SUBSAMPLE_SIZE) cells, a number that corresponds to 10% of the
# cells in the experiment.
#
# The experiment has 306144 control cells and 148649 cells treated
# with the ground-truth set. The total is 454794 cells.
#
# We compute 20 subsamples so that the methods that compute models
# based on subsamples can be run 20 times based on different
# subsamples.
#

SUBSAMPLE_SIZE = 45479
SUBSAMPLE = ../outputs/1.$(SUBSAMPLE_SIZE).controls.subsample

subsample: $(SUBSAMPLE)

$(SUBSAMPLE):
	python -m cpa.profiling.subsample --multiprocessing -pv $(NORM) \
	$(PROP) $(CACHE) $@ $(SUBSAMPLE_SIZE)

nsubsample: \
	../outputs/1.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/2.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/3.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/4.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/5.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/6.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/7.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/8.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/9.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/10.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/11.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/12.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/13.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/14.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/15.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/16.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/17.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/18.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/19.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/20.$(SUBSAMPLE_SIZE).controls.subsample \
	../outputs/1.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/2.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/3.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/4.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/5.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/6.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/7.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/8.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/9.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/10.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/11.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/12.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/13.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/14.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/15.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/16.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/17.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/18.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/19.$(SUBSAMPLE_SIZE).noncontrols.subsample \
	../outputs/20.$(SUBSAMPLE_SIZE).noncontrols.subsample

%.$(SUBSAMPLE_SIZE).controls.subsample:
	python -m cpa.profiling.subsample --multiprocessing -pv $(NORM) \
	-f controls \
	$(PROP) $(CACHE) $@ $(SUBSAMPLE_SIZE)

%.$(SUBSAMPLE_SIZE).noncontrols.subsample:
	python -m cpa.profiling.subsample --multiprocessing -pv $(NORM) \
	-f noncontrols \
	$(PROP) $(CACHE) $@ $(SUBSAMPLE_SIZE)

%.$(SUBSAMPLE_SIZE).both.subsample:
	python -m cpa.profiling.subsample --multiprocessing -pv $(NORM) \
	$(PROP) $(CACHE) $@ $(SUBSAMPLE_SIZE)

#
# Rules for computing preprocessors (factor analysis, PCA, and
# factor-analysis--based feature selection) based on subsampled cells.
#

../outputs/1.30614.controls.50.fasel.preprocessor: $(SUBSAMPLE)
	python -m cpa.profiling.factor_analysis \
	--variable-selection-only $< 50 $@

../outputs/1.30614.controls.50.pca.preprocessor: $(SUBSAMPLE)
	python -m cpa.profiling.pca $< 50 $@

#
# Run core profiling methods to produce per-well profiles
#

confusion: confusion-mean confusion-ks confusion-svm confusion-gm
confusion-mean: \
	../outputs/mean.treatment.confusion.pdf \
	../outputs/fa/50.fa.mean.treatment.confusion.pdf \
	../outputs/50.pca.mean.treatment.confusion.pdf \
	../outputs/50.fasel.mean.treatment.confusion.pdf
confusion-ks: \
	../outputs/ksstatistic.treatment.confusion.pdf \
	../outputs/50.fa.ksstatistic.treatment.confusion.pdf \
	../outputs/50.pca.ksstatistic.treatment.confusion.pdf \
	../outputs/50.fasel.ksstatistic.treatment.confusion.pdf
confusion-svm: \
	../outputs/svm.treatment.confusion.pdf \
	../outputs/50.fa.svm.treatment.confusion.pdf \
	../outputs/50.pca.svm.treatment.confusion.pdf \
	../outputs/50.fasel.svm.treatment.confusion.pdf
confusion-gm: \
	../outputs/gmm.treatment.confusion.pdf \
	../outputs/50.fa.gmm.treatment.confusion.pdf \
	../outputs/50.pca.gmm.treatment.confusion.pdf \
	../outputs/50.fasel.gmm.treatment.confusion.pdf
confusion-fa: \
	../outputs/fa/$(SUBSAMPLE_SIZE).controls.50.fa.mean.treatment.confusion.pdf

../outputs/mean.image.profiles.txt:
	python -m cpa.profiling.profile_mean $(PARALLEL)  \
	-o $@ -f noncontrols $(PROP) $(CACHE) Image

../outputs/mean.well.profiles.txt:
	python -m cpa.profiling.profile_mean $(PARALLEL)  \
	-o $@ -f noncontrols $(PROP) $(CACHE) Well

../outputs/mean+std.well.profiles.txt:
	python -m cpa.profiling.profile_mean $(PARALLEL)  \
	-o $@ -f noncontrols --method mean+std $(PROP) $(CACHE) Well

../outputs/mode.well.profiles.txt:
	python -m cpa.profiling.profile_mean $(PARALLEL_LONG) \
	-o $@ -f noncontrols --method mode $(PROP) $(CACHE) Well

../outputs/median.well.profiles.txt:
	python -m cpa.profiling.profile_mean $(PARALLEL)  \
	-o $@ -f noncontrols --method median $(PROP) $(CACHE) Well

../outputs/median+mad.well.profiles.txt:
	python -m cpa.profiling.profile_mean $(PARALLEL)  \
	-o $@ -f noncontrols --method median+mad $(PROP) $(CACHE) Well

../outputs/ksstatistic.well.profiles.txt:
	python -m cpa.profiling.profile_ksstatistic $(PARALLEL_LONG) \
	-o $@ -f noncontrols $(PROP) $(CACHE) Well controls

../outputs/svmnormalvector.well.profiles.txt:
	python -m cpa.profiling.profile_svmnormalvector $(PARALLEL) \
	-o $@ -f noncontrols $(PROP) $(CACHE) Well controls

../outputs/svmnormalvector_rfe.well.profiles.txt:
	python -m cpa.profiling.profile_svmnormalvector --rfe $(PARALLEL_LONG) \
	-o $@ -f noncontrols $(PROP) $(CACHE) Well controls

# ... and preprocessed

dimred: \
	../outputs/mean.treatment.confusion \
	../outputs/fa/$(SUBSAMPLE_SIZE).controls.50.fa.mean.treatment.confusion \
	../outputs/50.pca.mean.treatment.confusion \
	../outputs/50.fasel.mean.treatment.confusion \
	../outputs/ksstatistic.treatment.confusion \
	../outputs/50.fa.ksstatistic.treatment.confusion \
	../outputs/50.pca.ksstatistic.treatment.confusion \
	../outputs/50.fasel.ksstatistic.treatment.confusion \
	../outputs/svmnormalvector.treatment.confusion \
	../outputs/50.fa.svmnormalvector.treatment.confusion \
	../outputs/50.pca.svmnormalvector.treatment.confusion \
	../outputs/50.fasel.svmnormalvector.treatment.confusion \
	../outputs/$(SUBSAMPLE_SIZE).both.25.gmm.mean.treatment.confusion \
	../outputs/50.fa.gmm.treatment.confusion \
	../outputs/50.pca.gmm.treatment.confusion \
	../outputs/50.fasel.gmm.treatment.confusion

%.mean.well.profiles.txt: %.preprocessor
	python -m cpa.profiling.profile_mean $(PARALLEL)  \
	--preprocess $< \
	-o $@ -f noncontrols $(PROP) $(CACHE) Well

%.ksstatistic.well.profiles.txt: %.preprocessor
	python -m cpa.profiling.profile_ksstatistic $(PARALLEL_LONG8)  \
	--preprocess $< \
	-o $@ -f noncontrols $(PROP) $(CACHE) Well controls

%.svmnormalvector.well.profiles.txt: %.preprocessor
	python -m cpa.profiling.profile_svmnormalvector $(PARALLEL) \
	--preprocess $< \
	-o $@ -f noncontrols $(PROP) $(CACHE) Well controls

%.gmm.well.profiles.txt: %.preprocessor
	python -m cpa.profiling.profile_gmm $(PARALLEL) \
	--preprocess $< -o $@ --components 25 \
	$(PROP) $(CACHE) $(SUBSAMPLE) Well

# Combine per-well profiles into per-treatment and per-MOA profiles

%.treatment.profiles.txt: %.well.profiles.txt
	python -m cpa.profiling.median_profiles -o $@ $(PROP) $< \
	CompoundConcentration

%.moa.profiles.txt: %.well.profiles.txt
	python -m cpa.profiling.median_profiles -o $@ $(PROP) $< \
	MOA

#
#
#

FIGURES = \
	../figures/fa_varyfactors.pdf \
	../figures/confusion-fa.pdf \
	fig_confusion \
	../figures/distance-concentration.pdf \
	fig_confusion-population-statistics \
	../figures/gmm_varycomponents.pdf \
	fig_confusion-fa-others

figures: $(FIGURES)

clean_figures:
	rm -f $(FIGURES)

# Table 1: Accuracies for classifying compound treatments into mechanisms of action
tab_accuracies: \
	../outputs/mean.treatment.confusion \
	../outputs/ksstatistic.treatment.confusion \
	../outputs/svmnormalvector.treatment.confusion \
	../outputs/svmnormalvector_rfe.treatment.confusion \
	../outputs/$(SUBSAMPLE_SIZE).both.25.gmm.mean.treatment.confusion \
	../outputs/fa/$(SUBSAMPLE_SIZE).controls.50.fa.mean.treatment.confusion

# Figure 2: Distributions of classification accuracies for 20 runs of the factor analysis method for each possible choice of the number of factors from 2 to 100.
../figures/fa_varyfactors.pdf:
	./make_fa_varyfactors.sh
	python fa_varyfactors_plot.py $@

# Figure 3: Confusion matrix for the factor-analysis method.
../figures/confusion-fa.pdf: ../outputs/fa/$(SUBSAMPLE_SIZE).controls.50.fa.mean.treatment.confusion
	python pretty-confusion-matrix.py -f $< $@

# Supplementary Fig. 2: Histograms of distances between compound-concentrations that have the same MOA and between compound- concentrations that have different MOA.
../figures/distance-concentration.pdf: ../outputs/mean.treatment.profiles.txt
	python -m cpa.profiling.inter_intra_distances -o $@ $(PROP) $< MOA

# Supplementary Fig. 4: Confusion matrices for the five profiling methods.
fig_confusion: \
	../figures/confusion_a.pdf \
	../figures/confusion_b.pdf \
	../figures/confusion_c.pdf \
	../figures/confusion_d.pdf \
	../figures/confusion_e.pdf \
	../figures/confusion_f.pdf

../figures/confusion_a.pdf: ../outputs/mean.treatment.confusion
	python pretty-confusion-matrix.py -t 'A: Means' $< $@

../figures/confusion_b.pdf: ../outputs/ksstatistic.treatment.confusion
	python pretty-confusion-matrix.py -t 'B: KS statistics' $< $@

../figures/confusion_c.pdf: ../outputs/svmnormalvector.treatment.confusion
	python pretty-confusion-matrix.py -t 'C: SVM normal vector' $< $@

../figures/confusion_d.pdf: ../outputs/svmnormalvector_rfe.treatment.confusion
	python pretty-confusion-matrix.py -t 'D: SVM normal vector with recursive feature elimination' $< $@

../figures/confusion_e.pdf: ../outputs/$(SUBSAMPLE_SIZE).both.25.gmm.mean.treatment.confusion
	python pretty-confusion-matrix.py -f -t 'E: Gaussian mixture (mean of 20 models)' $< $@

../figures/confusion_f.pdf: ../outputs/fa/$(SUBSAMPLE_SIZE).controls.50.fa.mean.treatment.confusion
	python pretty-confusion-matrix.py -f -t 'F: Factor analysis (mean of 20 models)' $< $@

# Supplementary Fig. 7: Confusion matrices for each method, using factor analysis as a preprocessing step.
fig_confusion-fa-others: \
	../figures/confusion-fa-others_a.pdf \
	../figures/confusion-fa-others_b.pdf \
	../figures/confusion-fa-others_c.pdf \
	../figures/confusion-fa-others_d.pdf

../figures/confusion-fa-others_a.pdf: ../outputs/fa/$(SUBSAMPLE_SIZE).controls.50.fa.mean.treatment.confusion
	python pretty-confusion-matrix.py -f -t 'A: Means (mean of 20 models)' $< $@

../figures/confusion-fa-others_b.pdf: ../outputs/50.fa.ksstatistic.treatment.confusion
	python pretty-confusion-matrix.py -t 'B: KS statistics' $< $@

../figures/confusion-fa-others_c.pdf: ../outputs/50.fa.svmnormalvector.treatment.confusion
	python pretty-confusion-matrix.py -t 'C: SVM normal vector' $< $@

../figures/confusion-fa-others_d.pdf: ../outputs/50.fa.gmm.treatment.confusion
	python pretty-confusion-matrix.py -f -t 'D: Gaussian mixture' $< $@

#
# Plot confusion matrices
#

%.confusion: %.profiles.txt
	python -m cpa.profiling.leave_one_out -H Compound $(PROP) $< MOA > $@

%.sva.confusion: %.profiles.txt
	python -m cpa.profiling.leave_one_out -s -H Compound $(PROP) $< MOA > $@

%.gmm.mean.treatment.confusion.pdf: %.gmm.mean.treatment.confusion
	python pretty-confusion-matrix.py -f $< $@

%.confusion.pdf: %.confusion
	python pretty-confusion-matrix.py $< $@

# Supplementary Fig. 3: Confusion matrices for other population statistics.
fig_confusion-population-statistics: \
	../figures/confusion-population-statistics_a.pdf \
	../figures/confusion-population-statistics_b.pdf \
	../figures/confusion-population-statistics_c.pdf \
	../figures/confusion-population-statistics_d.pdf

../figures/confusion-population-statistics_a.pdf: ../outputs/mean+std.treatment.confusion
	python pretty-confusion-matrix.py -t 'A: Means' $< $@

../figures/confusion-population-statistics_b.pdf: ../outputs/mode.treatment.confusion
	python pretty-confusion-matrix.py -t 'B: Mode' $< $@

../figures/confusion-population-statistics_c.pdf: ../outputs/median.treatment.confusion
	python pretty-confusion-matrix.py -t 'C: Median' $< $@

../figures/confusion-population-statistics_d.pdf: ../outputs/median+mad.treatment.confusion
	python pretty-confusion-matrix.py -t 'D: Median + MAD' $< $@


# Supplementary Fig. 5: Distributions of classification accuracies for 20 runs of the Gaussian mixture method 
../figures/gmm_varycomponents.pdf:
	./make_gmm_varycomponents.sh
	python gmm_varycomponents_plot.py $@

# Supplementary Fig. 6: The GM method performs similarly whether the model is built from a subsample of control cells, a subsample of non-control (treated) cells, or a mixture of both.
figure_s6new: \
	../figures/gmm_fromcontrols.pdf \
	../figures/gmm_fromnoncontrols.pdf \
	../figures/gmm_fromboth.pdf

../figures/gmm_fromcontrols.pdf: ../outputs/$(SUBSAMPLE_SIZE).controls.25.gmm.mean.treatment.confusion
	python pretty-confusion-matrix.py -f -t 'A: Controls' $< $@
../figures/gmm_fromnoncontrols.pdf: ../outputs/$(SUBSAMPLE_SIZE).noncontrols.25.gmm.mean.treatment.confusion
	python pretty-confusion-matrix.py -f -t 'B: Non-controls'  $< $@
../figures/gmm_fromboth.pdf: ../outputs/$(SUBSAMPLE_SIZE).both.25.gmm.mean.treatment.confusion
	python pretty-confusion-matrix.py -f -t 'C: Both'  $< $@

../outputs/$(SUBSAMPLE_SIZE).controls.25.gmm.mean.treatment.confusion: 
	bash make_gmm.sh controls
	python mean_confusion_gmm.py $@ controls
../outputs/$(SUBSAMPLE_SIZE).noncontrols.25.gmm.mean.treatment.confusion: 
	bash make_gmm.sh noncontrols
	python mean_confusion_gmm.py $@ noncontrols
../outputs/$(SUBSAMPLE_SIZE).both.25.gmm.mean.treatment.confusion: 
	bash make_gmm.sh both
	python mean_confusion_gmm.py $@ both

# Supplementary Fig. 8: The FA method performs similarly whether the model is built from a subsample of control cells, a subsample of non-control (treated) cells, or a mixture of both.
figure_s8new: \
	../figures/fa_fromcontrols.pdf \
	../figures/fa_fromnoncontrols.pdf \
	../figures/fa_fromboth.pdf

../figures/fa_fromcontrols.pdf: ../outputs/fa/$(SUBSAMPLE_SIZE).controls.50.fa.mean.treatment.confusion
	python pretty-confusion-matrix.py -f -t 'A: Controls' $< $@
../figures/fa_fromnoncontrols.pdf: ../outputs/fa/$(SUBSAMPLE_SIZE).noncontrols.50.fa.mean.treatment.confusion
	python pretty-confusion-matrix.py -f -t 'B: Non-controls'  $< $@
../figures/fa_fromboth.pdf: ../outputs/fa/$(SUBSAMPLE_SIZE).both.50.fa.mean.treatment.confusion
	python pretty-confusion-matrix.py -f -t 'C: Both'  $< $@

../outputs/fa/$(SUBSAMPLE_SIZE).controls.50.fa.mean.treatment.confusion:
	bash make_fa.sh controls
	python mean_confusion_fa.py $@ controls
../outputs/fa/$(SUBSAMPLE_SIZE).noncontrols.50.fa.mean.treatment.confusion:
	bash make_fa.sh noncontrols
	python mean_confusion_fa.py $@ noncontrols
../outputs/fa/$(SUBSAMPLE_SIZE).both.50.fa.mean.treatment.confusion:
	bash make_fa.sh both
	python mean_confusion_fa.py $@ both

# Supplementary Fig. 8: Confusion matrices from the data in the Loo paper

fig_loo-confusion: \
	../figures/loo-confusion_a.pdf \
	../figures/loo-confusion_b.pdf \
	../figures/loo-confusion_c.pdf \
	../figures/loo-confusion_d.pdf

../outputs/loo_dna_anilin_sc35.confusion: ../inputs/Loo\ confusion\ matrices\ -\ DNA-Anilin-SC35.csv
	python loo_confusion.py "$<" $@

../outputs/loo_dna_mt_actin.confusion: ../inputs/Loo\ confusion\ matrices\ -\ DNA-MT-Actin.csv
	python loo_confusion.py "$<" $@

../outputs/loo_dna_p53_cFOS.confusion: ../inputs/Loo\ confusion\ matrices\ -\ DNA-p53-cFOS.csv
	python loo_confusion.py "$<" $@

../outputs/loo_dna_pp38_pERK.confusion: ../inputs/Loo\ confusion\ matrices\ -\ DNA-pp38-pERK.csv
	python loo_confusion.py "$<" $@

../figures/loo-confusion_a.pdf: ../outputs/loo_dna_anilin_sc35.confusion
	python pretty-confusion-matrix.py -t 'A: DNA–anilin–SC35' $< $@

../figures/loo-confusion_b.pdf: ../outputs/loo_dna_p53_cFOS.confusion
	python pretty-confusion-matrix.py -t 'B: DNA–p53–cFOS' $< $@

../figures/loo-confusion_c.pdf: ../outputs/loo_dna_mt_actin.confusion
	python pretty-confusion-matrix.py -t 'C: DNA–MT–Actin' $< $@

../figures/loo-confusion_d.pdf: ../outputs/loo_dna_pp38_pERK.confusion
	python pretty-confusion-matrix.py -t 'D: DNA–pp38–pERK' $< $@


#
# Tables
#

TABLES = \
	tab_accuracies \
	../tables/compounds_table.tex \
	#../tables/features.docx \
	../tables/features.tex \
	../tables/ground_truth_table.tex \
	../tables/mean_moa_variables.tex \
	tab_dimred \
	../tables/fa_loadings.tex \
	tab_moa-batches \
	../outputs/aucs_pvalues_table.tex

tables: $(TABLES)

clean_tables:
	rm -f $(TABLES) $(AUCS_PVALUES)

# Supplementary Tab. 1: Compounds table
# (depends on database not included in Supplementary Data)
../tables/compounds_table.tex: 
	python make_compounds_table.py $@

# Supplementary Tab. 2: Features table
../tables/features.docx:
	python make_features_table.py $@
../tables/features.tex:
	python make_features_table.py $@

# Supplementary Tab. 3: Ground-truth set
# (depends on database not included in Supplementary Data)
../tables/ground_truth_table.tex:
	python make_ground_truth_table.py $@

# Supplementary Tab. 4: Scoring by eye
# Compiled manually from observation by eye, not reproducible

# Supplementary Tab. 5: The 15 most important image features for distinguishing each mechanism of action using the means method.
../tables/mean_moa_variables.tex: ../outputs/mean.moa.profiles.txt
	python -m cpa.profiling.rank_variables --latex -o $@ $(PROP) $<

# Supplementary Tab. 6: Features most commonly selected by SVMRFE
# Generated in an interactive session, not reproducible

# Supplementary Tab. 7: Accuracies for all combinations of dimensionality-reducing preprocessing method and profiling method
tab_dimred: \
	../outputs/mean.treatment.confusion \
	../outputs/ksstatistic.treatment.confusion \
	../outputs/svmnormalvector.treatment.confusion \
	../outputs/$(SUBSAMPLE_SIZE).both.25.gmm.mean.treatment.confusion \
	../outputs/fa/$(SUBSAMPLE_SIZE).controls.50.fa.mean.treatment.confusion \
	../outputs/50.fa.ksstatistic.treatment.confusion \
	../outputs/50.fa.svmnormalvector.treatment.confusion \
	../outputs/50.fa.gmm.treatment.confusion \
	../outputs/50.pca.mean.treatment.confusion \
	../outputs/50.pca.ksstatistic.treatment.confusion \
	../outputs/50.pca.svmnormalvector.treatment.confusion \
	../outputs/50.pca.gmm.treatment.confusion \
	../outputs/50.fasel.mean.treatment.confusion \
	../outputs/50.fasel.ksstatistic.treatment.confusion \
	../outputs/50.fasel.svmnormalvector.treatment.confusion \
	../outputs/50.fasel.gmm.treatment.confusion

# Supplementary Tab. 8: The 15 features most heavily loaded onto each factor in the 50-factor model
../tables/fa_loadings.tex: ../outputs/1.30614.controls.50.fa.preprocessor
	python -m cpa.profiling.factor_loadings --latex -o $@ $<

# Supplementary Tab. 9: Distribution of mechanisms of action across batches
tab_moa-batches: ../outputs/batch_effects.png

../outputs/batch_effects.png:
	python batch_effects.py -o $@ $(PROP)

# Supplementary Tab. 10: AUCs and p-values
../outputs/aucs_pvalues_table.tex: $(AUCS_PVALUES)
	python make_aucs_pvalues_table.py $@

AUCS_PVALUES = \
 ../outputs/mean.treatment.aucs_pvalues \
 ../outputs/ksstatistic.treatment.aucs_pvalues \
 ../outputs/svmnormalvector.treatment.aucs_pvalues \
 ../outputs/gmm.treatment.aucs_pvalues \
 ../outputs/factoranalysis_mean.treatment.aucs_pvalues \

%.aucs_pvalues: %.profiles.txt
	python calculate_aucs.py $(PROP) $< MOA > $@

# Plot factor loadings
%.fa_loadings.txt: %.fa.preprocessor
	python -m cpa.profiling.factor_loadings -o $@ $<

%.fa_loadings.tex: %.fa.preprocessor
	python -m cpa.profiling.factor_loadings --latex -o $@ $<

# Classify DMSO wells

../outputs/mean.well.dmsoprofiles.txt:
	python -m cpa.profiling.profile_mean $(PARALLEL)  \
	-o $@ -f controls $(PROP) $(CACHE) Well

%.mean.well.dmsoprofiles.txt: %.preprocessor
	python -m cpa.profiling.profile_mean $(PARALLEL)  \
	--preprocess $< \
	-o $@ -f controls $(PROP) $(CACHE) Well

../outputs/mean.classify_dmso.txt: ../outputs/mean.well.dmsoprofiles.txt ../outputs/mean.treatment.profiles.txt
	python classify_dmso.py $@ $(PROP) $^

../outputs/fa.classify_dmso.txt: ../outputs/1.30614.controls.50.fa.mean.well.dmsoprofiles.txt ../outputs/1.30614.controls.50.fa.mean.well.profiles.txt
	python classify_dmso.py $@ $(PROP) $^



#
# Not included in the paper
#


# Make list of misclassified treatments
%.misclassified.txt: %.treatment.profiles.txt
	python misclassified.py -H Compound $(PROP) $< MOA > $@

# Plot distance matrix
%.distances.png: %.profiles.txt
	python -m cpa.profiling.plot_distances -o $@ $(PROP) $< MOA

# Plot profiles

profile_figures: \
 ../gmm.treatment.profiles.pdf \
 ../factoranalysis_mean.treatment.profiles.pdf

%.treatment.profiles.pdf: %.treatment.profiles.txt
	python -m cpa.profiling.plot_profiles -o $@ $(PROP) $< MOA

%.moa.profiles.pdf: %.moa.profiles.txt
	python -m cpa.profiling.plot_profiles -o $@ $(PROP) $< MOA

%.inter_intra.pdf: %.profiles.txt
	python -m cpa.profiling.inter_intra_distances -o $@ $(PROP) $< MOA

%.variables.txt: %.profiles.txt
	python -m cpa.profiling.rank_variables -o $@ $(PROP) $<

%.variables.tex: %.profiles.txt
	python -m cpa.profiling.rank_variables --latex -o $@ $(PROP) $<

../tables/fa_moa_variables.tex: ../outputs/1.30614.controls.50.fa.mean.moa.variables.tex
	cp $< $@

%.kaiser.txt: %.profiles.txt	
	python -m cpa.profiling.kaiser -o $@ $<

# Compare inter- and inter-MOA distances with DNA only.

PROP_DNAONLY = ../properties/dnaonly.properties
CACHE_DNAONLY = ../cache_dnaonly

$(CACHE_DNAONLY):
	python -m cpa.profiling.cache $(PROP_DNAONLY) $@ "Image_Metadata_Compound = 'DMSO'"

../outputs/dnaonly.well.profiles.txt: $(CACHE_DNAONLY)
	python -m cpa.profiling.profile_mean $(PARALLEL)  \
	-o $@ -f noncontrols $(PROP_DNAONLY) $(CACHE_DNAONLY) Well

../figures/distance-concentration-dnaonly.pdf: ../outputs/dnaonly.treatment.profiles.txt
	python -m cpa.profiling.inter_intra_distances -o $@ $(PROP_DNAONLY) $< MOA
