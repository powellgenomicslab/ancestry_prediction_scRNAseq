Results!
=========
After running those jobs, you should be done! 


You should have the following results directories:

.. code-block:: bash



- There will be an ``Ancestry_PCAs.png`` figure generated for each individual in each pool.

  - This figure shows the 1000G individual locations in PC space compared to the individuals in each pool. For example:

  .. figure:: _figures/
      :align: right
      :figwidth: 300px
    

- There will be an ``ancestry_assignments.tsv`` file generated for each individual in each pool and one that has all the individuals joined together in the base output directory.

    - This file has the annotations and probabilities for each pool. For example:



Additional Results with Reference-based Exectution
----------------------------------------------------

If you have reference SNP genotypes (*i.e.* microarray or whole exome or genome sequencing-called SNPs) and decided to estimate SNP-based ancestry, you will have additional results that compare the reference and single-cell ancestry predictions.


These are the main comparison files with the two most informative ones highlighted:

.. code-block:: bash
  :emphasize-lines: 6,7

    snp_ancestry_predictions.tsv
    reference_ancestry_numbers.png
    predicted_ancestry_numbers_correct.png
    predicted_ancestry_numbers_correct_identified.png
    statistics_heatmap.png
    predicted_numbers_statistics_heatmap_combined.png
    assignments_probabilities_w_ref.png


The ``predicted_numbers_statistics_heatmap_combined.png`` figure shows the number of individuals classified to each ancestry with the single cell data and colored by if they correctly or incorrectly match the reference-annotated SNP genotype data
and the heatmap below shows some statistical metrics for quantifying the single cell derived ancestry predictions:

  .. figure:: _figures/
      :align: right
      :figwidth: 300px


The ``assignments_probabilities_w_ref.png`` figure show the probability of each sample to be classified to each of the different ancestries. 
This includes the reference-based predictions (microarray or whole exome or genome sequencing data) compared to the single cell based predictions

  .. figure:: _figures/
      :align: right
      :figwidth: 300px