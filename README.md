# Supplementary Materials

Here, I provide the scripts used in the computational pipeline for detecting and classifying transposable elements (TEs). An outline of the pipeline can be seen in the figure below 

![FlowChart_pipeline](https://github.com/uio-cels/Repeats/blob/master/figures/FlowChart_pipeline.png)

The scripts used for creating species-specific libraries of TEs are called 'repeats\_master\_pipeline.slurm' and 'repeats\_worker\_script.slurm' and is run like this:

```
sbatch repeats_master_pipeline.slurm genome.fasta
```

The pipeline runs around a week on medium sized genomes (0.4 - 1.5 Gb).

In addition, the data analysis steps for making the figures can be found in 'notebooks'. Parts of the data frames that were analyzed can be found in 'data' (unfortunately some are too big to be shared entirely). 




