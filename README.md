#Genepop2SPAM.GCL
**An R paqckage for converting Genepop files to SPAM baseline and mixture files.**

Alaska Dept. of Fish & Game

Gene Conservation Laboratory


_This function is currently in a testing stage. Please send all non-working examples to charles.thompson@alaska.gov for review if possible._

**About:**

The following table represents how Genepop stores ordered sample alleles differently between the SPAM mixture file for a SNP.

|  | Genepop | SPAM |
| --- | --- | --- |
| A/a | 0101 | 20 |
| A/a | 0102 | 11 |
| a/a | 0202 | 02 |

For mixture files, this function changes the ids of the allele listed in the Genepop file to the counts of the allele in a number in a vector as long as the number of alleles.

**Installation:**

To use this tool, install the most recent version of [R]  (http://www.r-project.org/). After loading into R install the additional package, "stringr". In the console enter:

```
>install.packages("stringr")
>library(devtools)
>install_github("hansthompson", "Genepop2SPAM")
>library(Genepop2SPAM)
>?Genepop2SPAM.GCL
```

**Considerations For Each Estimation and Simulation:**

1. This function is set up to work only with diploid alleles for SNPs and microsatellites or mtDNA haplotypes with ninor less SNPs at this time.

2. All the files for the baseline and mixtures should be read in at the same time to create SPAM files. This is because for each locus it will need to know how many different alleles are present across all files to properly id each one into the correct position (see the About section table). Alleles with ids that are not in a consecutive order for the same locus in Genepop format (for an ordered version of sample allele names in a microsatellite locus ie 001001, 002002, 004004 …) are fine. However, this function will give them new corresponding id positions.  mtDNA haplotypes will be recognized as single alleles and given corresponding ids as well. 

3. Each population in the Genepop file should be separated by POP, Pop, or pop. Any other population separator will not be recognized. 

4. The locus names listed at the top of each Genepop file should be in the same order and have the same number of loci across all files.

5. Within the Genepop file there should only be one leading zero in front of the allele id for each data point.

6. Samples with no data for any loci should be represented as 0's for the same character length of the other samples in the locus.

**Application:**

Here is an example of how to use the function:

* Bring the function, DevGen2SPAM.GCL, into the R environment by entering into the console,
```
library(Genepop2SPAM)
```
* Organize your Genepop files

a. Find the working directory and place the genepop file/s in this directory.
```
>getwd()
```

something similar to "C:/Users/user/Documents" should appear and be where you should place the genepop files.

b. If you would like to run the function from a different folder outside this example, you can change the working directory.
```
>setwd(yourDirectory)
```

c. Or enter the full file path into the function.

* Run the function without labeling the full file paths. If you have multiple mixtures you would like to test against the baseline, add all of them at the same time to the function. To do this, wrap the file names in the concatenation function c with quotes surrounding each file and a comma separating each quotation:
```
>c("file1", "file2", "file3")```
#So it should run as so (without the fullpath for this example):
>Genepop2SPAM.GCL(baselineFile = "file1.gen", mixtureFiles = c("file2.gen", "file3.gen", "file4.gen")
```
Note that baselineFile is singular and mixtureFiles can be any number of files. 

This should write your baseline file and three mixture files using the same file names that were provided. If you would like to put new names for the files, additional arguments baselineOutPath and baselineOutPath can be used to match the number of files that were used as inputs. It is important that
