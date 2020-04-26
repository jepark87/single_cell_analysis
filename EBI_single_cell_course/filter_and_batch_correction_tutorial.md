# Preprocessing single-cell data and batch correction

## Objectives

By the end of this short tutorial should you should understand:

 * how to check the quality of single-cell experiment and clean it
 * how to preprocess dataset for analysis (normalisation, log-transform, scaling)
 * how to perform demensionality reduction (PCA, umap) and clustering
 * how to perform batch correction

## Input files required

You will most likely have multiple files from different experiments (like control vs treatment, multiple samples collected from different individuals, or run on different single-cell methods), which needs to be concatenated

 * matrix.mtx, barcodes.tsv, genes.tsv files from alignment tool (star-solo, cellranger, salmon alevin etc)
 
For the sake of time, we will start with toy example (*.h5ad format) which have concatenated multiple single-cell experiment, containing data from several donors and two different 10X single-cell chemistries.

 * adata.h5ad file

# Introduction

In the previous session, you have learned how to generate matrix files from fastq files. This session will guide you to how to check the quality of these matrix files and analyse them. This includes: checking the quality of the generated matrix files, preprecessing the data for analysis, correcting batch effect and visualising them by clustering. This is most basic workflow that you will run in any single-cell data analysis.

There are many packages, which have collection of functions for running this workflow. Such as [Seurat](https://satijalab.org/seurat/), [Scanpy](https://icb-scanpy.readthedocs-hosted.com/en/stable/), [Monocle3](https://cole-trapnell-lab.github.io/monocle3/), [Scater](http://bioconductor.org/packages/release/bioc/html/scater.html). They offer streamlined workflow and tutorials that can be easily followed.

Today, we will follow the [scanpy workflow](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html). Scanpy is python package based on object called [anndata](https://icb-anndata.readthedocs-hosted.com/en/stable/anndata.AnnData.html), which is data storage format made easy to store metadata and processed data from each stages.

We'll be using a setup based on Galaxy to illustrate the process for teaching purposes, but all steps are peformed using commodity tools you could run from the command line if you're more comforable there.

## 1. Example data

The data you will process today is small toy example taken from [human thymus dataset](https://zenodo.org/record/3711134#.XqPKJVNKhrB). This study is good example for batch correction as the data has collected from multiple donors and using two different 10X chemistries (3' v2 chemistry and 5' chemistry).

Down-sampled reads and some associated annotation are provided for you in what Galaxy calls a 'history'. Access the shared history:

![Getting to the shared data](goto_histories.png)

![Select the specific history](specific_history.png)


If you click through to the history itself, you'll find the 'toy_example_thymus.h5ad' file at the bottom. This file contains concatenated matrix from 5 individual 10X lanes.

Batch | Donor | 10X chemistry
--------------|-------|--------------
1|F21|3'
2|F29|3'
3|F29|5'
4|F30|3'
5|F30|5'

All datasets have been mapped with [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger), which is alignment tool provided by 10X. The resulting count matrix have been loaded using [scanpy.read](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.read.html#scanpy.read) function and concatenated into single anndata object (h5ad format).

Metadata contains following contents:

 * donor: human donor ID
 * method: 10X chemistry used. 5GEX (10X 5') and 3GEX (10X 3' v2)
 * batch: donor x method

For general QC, we have calculated some basic features, which are as following:

 * n_genes: number of total genes detected per cell
 * n_loggenes: log10(n_genes)
 * n_counts: number of total UMI counts per cell
 * n_logcounts: log10(n_counts)
 * mito: fraction of total UMI counts assigned to mitochondrial genes to total UMI counts
 * doublet_scores: doublet score calculated by [scrublet](https://github.com/AllonKleinLab/scrublet) [(Wolock 2019)](https://www.cell.com/cell-systems/pdfExtended/S2405-4712(18)30474-5)
 
The dataset also contains cell type annotation, which was taken from thymus cell atlas, to aid visualistion of batch effect

 * celltype

## 2. The workflow

### 1. Checking the quality of data

So, we have data that has just come out from mapping pipeline. Now we need to check whether the run has been successful. There are many steps that can go wrong. Your cells might have died before loading, which would yield low UMI counts from cells and higher background. Sometimes cell counting might have failed, resulting in over-loading. In this case, you will have too many cells and doublets. Thus, it's important to perform general QC for each run. Otherwise, your end result will be deteriorated by these low-quality signals from failed run.

We will start by visualising the distribution of number of genes detected per cell, number of total UMI counts per cell, and fraction of mitochondiral mRNAs.

![Step1](Step1.png)



 1. Search 'plot' from Galaxy search toolbar
 2. From 'Scanpy' package, click 'Plot with scanpy'
 3. Link h5ad input with plot function
 4. Go to __Details__, set __Method used for plotting__ as __Generic: Scatter plot along observations or variable axes, using 'pl.scatter'__. 
 5. set __x- and y- coordinates__ with __'n_logcounts', 'n_loggenes'__ or __'n_loggenes', 'mito'__ or __'n_loggenes', 'doublet_scores'__'''
 


![Result1](Result1.png)

> Q: Which threshold would you use to filter bad cells? 

The plots show distribution of basic QC measures and their relationship. As expected, n_counts and n_genes are proportional. At the same time, you could notice droplets with disproportionally less n_genes compared to n_counts, which also have higher mito fraction. This is often indicative of empty droplet which was mistaken from initial filter and you would want to remove them.

In this dataset, most cells with good gene number coverge (high n_genes) contain less than 20% of mitochondrial mRNA. If you see mitochondrial fraction much higher than this (like 50~80%), this is probably resulting from dying cells.

You can notice high doublet scores for some droplets. You would generally expect ~4% doublet rate if you aimed to recover 5000 cells. Generally, doublets have higher n_genes, as they contain genes from both cell types. But some cells do have more n_genes and other cell in general. So you should be careful about setting doublet filter based on n_genes. (Also, some cells have significanly less n_genes compared to others. So be careful for lower bound filters too!)

![Result2](Result2_violin.png)

Scatter plot is good to compare relationship between QC measures. However, it is often difficult to compare distributions across multiple samples. For this purpose, let's draw some violin plots. You can simply change the step 4 and step 5 from above to:



 4. Go to __Details__, set __Method used for plotting__ as __Generic: Violin plot, using 'pl.violin'__.
 5. Set __Keys for accessing variables__ as __Subset of variables in 'adata.var_names' or fields of '.obs'__ and provide comma separte list of __n_loggenes,n_logcounts,mito__. Also, set __The key of the observation grouping to consider__ as __batch__.



> Q: Which sample has the lowest coverage? Are there any difference between 10X chemistries?

### 2. Filter cells based on QC metrics.

After visual inspection, We have come up with following conclusion: _Let's remove any droplets with less than 500 n_genes, 1000 n_counts, or higher than 0.2 mito fraction_. This can be simply achieved by built in function in scanpy called __filter_cells__. We will also remove genes that are not expressed in this dataset, by taking genes which are detected in at least 3 cells. There's another function for this, which is called __filter_genes__.


 1. Go to 'Tools' and within 'Scanpy' package, click 'filter_cells' or 'filter_genes'
 2. Set up consecutive filters. Set __parameters to select genes (cells) to keep__ in __Details__ as suggested above.
 3. Click each  




To generate gene-level quantifications based on transcriptome quantification, Alevin and similar tools require a conversion between transcript and gene identifiers, which we can derive from the annotations available in genome resources such as Ensembl. The transcripts in this list need to match the ones we used in the sequences reference. The study we're working with here added spike-ins to their samples, so when we created the index above in preparation for this tutorial, we combined sequence information for the biological sequences and the spike-ins. So now we need to generate a transcript to gene mapping with both types of sequence.

In your example data you will see 'Homo_sapiens.GRCh38.95.gtf.gz', which is the human reference annotations as retrieved from Ensembl in GTF format, and a transcript-level description of the ERCC spike-ins (also in GTF format). We will use these to generate the transcript/ gene mapping, first by combining the two types of annotation and second by passing that information to a tool that extracts just the identifiers we need:

Find the 'concatenate' tool in Galaxy under 'Text Manipulation':

![concatenating annotation files](concatenate.png)

Select biological and technical sequences from your history, and use them as inputs. Run the tool and generate your combined annotations.
 
> EXERCISE: Which of the 'attributes' in the last column of the GTF files contains the transcript and gene identifiers?

Now we have combined annotations, we can parse the GTF file using the [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html) package in R. We've made a convenience tool in Galaxy to wrap that functionality for you- search for 'GTF2GeneList':

![GTF2GeneList](gtf2genelist.png)

This tool is designed to do a variety of things. For our purposes parameterise like:

 * "Feature type..." : transcript
 * "Field to place first..." : (use the transcript ID field you identified above)
 * "Suppress header" : yes
 * "Comma separated list of field names" : (transcript ID field),(gene ID field)
 * "Append version to transcript identifiers" : yes 
 * "Flag mitochondrial features" : no
 * "Filter a FASTA": no

> EXERCISE: compare the transcript identifiers you'll find in the FASTA file you were given at the start with the transcript ID values in the GTF, to see why we say 'yes' to the 'transcript identifiers' question

### 3. Running quantification

Now we have everything we need to run Alevin:

 * A Salmon index for our transcriptome
 * Barcode/ UMI reads
 * cDNA reads
 * transcript/ gene mapping

We can now run Alevin. Locate the tool by searching with the 'search tools' box. 

![The Alevin tool](alevin.png)

> EXERCISE: Using [Alevin's documentation](https://salmon.readthedocs.io/en/latest/alevin.html), what you know of the protocol used in this experiment (see above), and what you identified by examining the read files, set the options appropriately and run Alevin to produce a matrix market format (MTX) format output. Try to determine the correct library options yourself, but do select the built-in index (building one yourself from this history will take longer than you have!), activate MTX output, select 'dumpFeatures' and do click 'Retrieve all output files'. Many options can be left unset.

Because we're only using a million or so reads Alevin will run quickly and finish within a few minutes. Alevin produces many file outputs not all of which we'll use. You can refer to the [Alevin documentation](https://salmon.readthedocs.io/en/latest/alevin.html) if you're curious what they all are, but we're most interested in the matrix itself (quants_mat.mtx.gz), the row (cell/ barcode) identifiers (quants_mat_rows.txt) and the column (gene) labels (quants_mat_cols.txt). 

> EXERCISE: Once you've run Alevin, look through the files and see if you can find: 1) the mapping rate 2) how many cells are present in the output.

### 4. Basic QC with barcode rank plots

Congratulations- you've made an expression matrix! We could almost stop here. But it's sensible to do some basic QC, and one of the things we can do is look at a barcode rank plot.

The question we're looking to answer here, is: "do we have mostly a have single cell per droplet"? That's what experimenters are normally aiming for, but it's not entirely straightforward to get exactly one cell per droplet. Sometimes almost no cells make it into droplets, other times we have too many droplets with two or more cells. But at a minimum, we should easily be able to distiguish droplets with cells from those without. 

Locate the barcode rank plot tool by searching for it in the search box. Select 'No' to input MTX, and select the 'raw_cb_frequencies.txt' file you should have in your history from running Alevin. If you do not, then you didn't select 'dumpFeatures' when you ran Alevin- so go back and try again. Set a title if you wish, but leave other options at defaults. 

![droplet barcode plot tool](droplet_barcode_tool.png)

You'll end up with a plot like: 

![barcode plot from raw barcode counts](barcodes_raw.png)

This is our own formulation of the barcode plot based on a [discussion](https://github.com/COMBINE-lab/salmon/issues/362#issuecomment-490160480) we had with community members. The left hand plot is the main one, showing the counts for individual cell barcodes ranked from high to low. We expect a sharp drop-off between cell-containing droplets and ones that are empty or contain only cell debris. The right hand plot is a density from the first one, and the thresholds are generated either using [dropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html) or by the method described in that discussion. We could use any of these thresholds to select cells, assuming that anything with fewer counts is not a valid cell. By default, Alevin does something similar, and we can learn something about that by plotting just the barcodes Alevin retains. Go back and re-run the droplet barcode plot, this time selecting MTX input (quants_mat.mtx.gz). You will need to select the option to assume cells are by row (more on that later). 

![droplet barcode plot tool](droplet_barcode_tool2.png)

This will use the actual sum of cell-wise counts produced in Alevin's outputs to make the plot. The output will be like:

![barcode plot from processed barcode counts](barcodes_processed.png)

You should see a completely vertical drop-off where Alevin has trunctated the distribution.

In experiments with relatively simple characteristics, this 'knee detection' method works relatively well. But some populations present difficulties due to sub-populations of small cells that cannot be distinguished from empty droplets based purely on counts by barcode. Some libraries produce multiple 'knees' for multiple sub-populations. The [emptyDrops](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y) method has become a popular way of dealing with this. emptyDrops still retains barcodes with very high counts, but also adds in barcodes that can be statistically distinguished from the ambient profiles, even if total counts are similar.

## 5. Apply emptyDrops to remove empty cells

To use emptyDrops effectively, we need to go back and re-run Alevin, stopping it from applying it's own thresholds. If you look under 'optional commands' you will see 'keepCBFraction'. Set this to 1 to retain all cell barcodes. Then set freqThreshold to 3: this will only remove cell barcodes with a frequency of less than 3, a low bar to pass but useful way of avoiding processing a bunch of almost certainly empty barcodes. Trigger the Alevin re-run.

> EXERCISE: How many cells are in the output now?

Alevin outputs MTX format, which we can pass to the dropletUtils package and run emptyDrops. Unfortunately the matrix is in the wrong orientation for tools expecting files like those produced by 10X software (which dropletUtils does). We need to 'transform' the matrix such that cells are in columns and genes are in rows. We've provided you with a tool for doing this, search for 'salmonKallistoMtxTo10x':

![matrix transformation tool](transform_tool.png)

Run this tool, supplying the inputs specified in the help text (make sure you don't mix files from different Alevin runs)). The output will be a matrix in the correct orientation to pass to the next step. 

emptyDrops works with a specific form of R object called a SingleCellExperiment. We need to convert our transformed MTX files into that form, using the DropletUtils Read10x tool:

![reading matrix files into a SingleCellExperiment](read_10x.png)

Now we have the data in the right format, we can run emptyDrops. Search for it with the tools search box:

![the emptyDrops tool](emptydrops_tool.png)

Be sure to select 'yes' when asked 'Should barcodes estimated to have no cells be removed from the output object?'.

> EXERCISE: how many cell barcodes remain after the emptyDrops treatment? Why might that be? (hint: is this a real/ complete set of data?). Go back and tweak parameters, re-running the tool. 

Assuming you've completed the last exercise, you have an expression matrix ready to go, in the SingleCellExperiment format of R. The other trainers will mostly be using a tool called Scanpy. They won't be using these dummy data, but if you wanted to pass this matrix to that tool, you would need to convert to a format called annData, which is a variant of a file format called hdf5. To help you with this we've provided you with a tool called sceasy:

![converting formats with sceasy](sceasy.png)

If you use the options as selected in the image above, you will have an annData object with raw data, ready to be passed to downstream analysis with Scanpy.

# 6. Running the whole workflow

In real life we do not run through analysis steps piecemeal as above, we run workflows in an automated manner. A complete workflow composed of the above steps is available to you via 'Shared Data' -> 'Workflows', labelled 'Droplet Quantification and preprocessing'. 

![get the whole workflow](workflow1.png)

Click the dropdown next to the workflow to import it. Then go back, select the 'Workflow' tab and click on the workflow. This will open the workflow editor and you can overview what we did above:

![the workflow editor](workflow2.png)

If you click the play icon you get an interface allowing you to execute the whole workflow. Try it out:

![the workflow runner](workflow3.png)


# 7. Summary

You've reached the end of this session. We have:

 * Taken raw read data and annotations and necessary input files for quantification.
 * Run Alevin in two different parameterisations, both allowing Alevin to make its own calls on what constitutes empty droplets, and applying emptyDrops instead.
 * Deployed barcode rank plots as a way of quickly assessing the signal present in droplet datasets. Further assessments of quality are the subject of subsequent sections of this course.
 * Applied the necessary conversion to pass these data to downstream processes. 
 * Run a workflow of all steps at once.

# References

 - Lun ATL, Riesenfeld S, Andrews T, et al. EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. Genome Biol. 2019;20(1):63.
 - Melsted P, Ntranos V, Pachter L. The Barcode, UMI, Set format and BUStools. Bioinformatics. 2019;
 - Srivastava A, Malik L, Smith T, Sudbery I, Patro R. Alevin efficiently estimates accurate gene abundances from dscRNA-seq data. Genome Biol. 2019;20(1):65.
 - Zheng GX, Terry JM, Belgrader P, et al. Massively parallel digital transcriptional profiling of single cells. Nat Commun. 2017;8:14049.
