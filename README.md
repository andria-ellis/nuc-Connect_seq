**Project Overview:** This GitHub repo is for the paper Nuc-Connect-Seq paper, which seeks to identify and sequence activated neuron's in mouse hypothalamus, and use their transcriptomes to find important marker genes or neuropeptides for mouse restraint stress response.

Analysis for this project was completed in the following 7 code blocks:

**1_combine_align_all_cells.R:** Combines and aligns all cell data.

**2_cell_type_clustering.R:** Performs clustering analysis on different cell types.

**3_analyze_buck_data_only.R:** Analyzes buck data in isolation.

**4_neuron_type_labelling.R:** Labels neuron types based on clustering results.

**5_glu_subtype_analysis.R:** Analyzes subtypes of glutamatergic neurons.

**6_gaba_subtype_analysis.R:** Analyzes subtypes of GABAergic neurons.

**7_regression_testing.R:** Conducts regression testing to validate findings.

These code blocks should be run in order, with intermediate files being created along the way. Each code block is commented on top with the summary, the input files needed, the output files created, and the figures created. 

**Installation instructions:** The following packages will need to be installed to run these code blocks. 

  install.packages('dplyr', 'plyr', 'ggplot2', 'monocle3', 'reshape2', 'Garnett', 'Seurat')

