## hmgu-itg-challenge
### Solutions
---

### Support/resource management/Shell

> 1. A user has several versions of R installed in their path. Each version of R has a number of locally installed libraries. The user is confused, and would like to know which library is installed for each version of R. Can you write a command to help them out?

```sh
which -a R | while read Rbin; do 
  echo "R at: $Rbin" 
  "$Rbin" -e '.libPaths(); installed.packages()[, c("Package", "Version")]'
  echo "----------------------"
done
```

> 2. A common problem with shared filesystems is a disk quota overflow. This can be due to 1) a large combined size on disk or 2) too many files present on disk. We would like to help users who encounter this problem to locate the problematic files/directories. Write a command to sort all subdirectories of level n (n determined by the user) by their human-readable size. Write another command to sort all subdirectories of level n according to the number of files they contain.
```sh
du -ah --max-depth=$n | sort -h
find . -mindepth $n -maxdepth $n -type d | xargs du -sh | sort -h
```

```sh
find . -mindepth $n -maxdepth $n -type d | while read d; do 
  echo "$(find "$d" -type f | wc -l) $d"
done | sort -n
```

> 4. A user is running commands like this one cat file1 <(cut -d " " -f 1-15,17,18 file2) > file3. What does this command do? It runs fine on the command line, but then the user includes it into a file with other commands, saves it and runs chmod +x on it. However, that line of code throws the following error : syntax error near unexpected token '('. What has the user forgotten?

The script could be missing the shebang line

> 6. Programmatic download
>- You have to download all autosomal files from this location: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ onto your server. You connect to the server via SSH. Using only the command line, how do you perform this download?
>- You are at a conference abroad and you quickly realise that your connection is unstable. You get disconnected constantly, which interrupts the download. How do you ensure the download survives these disconnections?

```sh
scp username@host:/path/remote/file /path/local/file
```
```sh
rsync -av --partial source_dir username@destinationhost:/destination_dir
```
```sh
screen -S mysession #Keeping sessions on
```

> 8. An analysis you need to run on the cluster requires a particular python library, but you do not have administrator rights. IT is on holiday. What do you do?

Local installation
```sh
pip install --user library_name

pip install --target=/home/python3/site-packages library_name

export PYTHONPATH="/path/to/local-site-packages:${PYTHONPATH}"
```
A virtual environment through the ```virtualenv``` package can also be created

---

### Bioinformatics

> 2. From an existing VCF with an arbitrary number of samples, how do you produce a VCF file without any samples using ```bcftools ```? 

```sh
bcftools view input.vcf.gz -G -Oz -o grouped_sample_out.vcf.gz
```

> 4. How do you convert a gzipped VCF to the bimbam format? 

```sh
plink --vcf vcf_file.gz --recode bimbam --out vcf_bimbam_out
```

> 5. A user sends you a small number of chromosome and positions in build 38 that they want to know the rsID of.
>- What is missing from their request? What kind of unexpected output can they expect?
>- Given this file, honour their request using the Ensembl REST API.
>- Do the same, but offline, using the dbSNP rs.150 VCF file.
>- What would change if these positions were in build 37?
>- If the user sends you 7,000 such chromosome positions, how would the above methods perform? Do you know of any alternatives?

- It is unclear whether to retrieve rsIDs for SNPs only, or both SNPs and indels.
- The reference database(s) to use (e.g., dbSNP, Ensembl, ClinVar) are not specified.
- The request does not indicate whether the chromosome positions are 1-based or 0-based, which could lead to mismatches.
Unexpected outputs can include incorrect IDs, or missing IDs.

```sh
curl -X GET "https://rest.ensembl.org/overlap/region/human/1:1553203-1553203?feature=variation" -H "Content-Type: application/json"
```
Example output: The first coordinate- 1:1553203 corresponds to ```rs6603790```
```json
[{"consequence_type":"intron_variant","source":"dbSNP","strand":1,"clinical_significance":[],"assembly_name":"GRCh38","feature_type":"variation","alleles":["C","T"],"start":1553203,"seq_region_name":"1","id":"rs6603790","end":1553203}]% 
```
The retrieved ```rs6603790``` when using dbSNP corresponds to chr1:1553103-1553303 for GRCh38
and chr1:1488483-1488683 on GRCh37.

- ```bedtools intersect``` with a dbSNP VCF
- Running a SQL query on whole vcf and ```join``` with input data.
- Querying ```biomaRt``` package on R 
- running the VCF through Ensembl VEP

> 6. How would you change the chromosome numbers in the file above to chromosome names (e.g. "chr1" instead of "1")?
>- How would you change the names back to the original? Would your solution work if an additional column containing text of arbitrary length and content is appended at the left of the file?
>- These positions are extracted from a VCF. Convert this file to the BED format.

```sh
awk 'BEGIN{OFS=FS=" "} {$2="chr"$2; print}' rand.chrpos.txt > rand.chrpos_nochr.txt
```
Removing `"chr"` prefix:
```sh
awk '{$2=substr($2,4); print}' rand.chrpos_nochr.txt > rand.chrpos_chr.txt
```
```sh
awk '!/^#/ {print $1, $2-1, $2-1}' OFS='\t' rand.chrpos.txt > rand.chrpos.txt.bed
```

> 9. We want to round a column of numbers to n decimal places, with values with 5 as their rightmost significant digit rounded up. Use the language of your choice.

```sh
awk -v n=2 '{printf "%.*f\n", n, $1 + (0.5 / (10^n))}' input.txt
```

> 10. Is this HRC-imputed file missing any chromosomes? Try to find out in seconds if you can.

```sh
gzcat hrc.positions.txt.bgz|  awk '                                                                       
  BEGIN { for (i=1; i<=22; i++) expected[i] = i; expected["X"]="X"; expected["Y"]="Y"; expected["MT"]="MT"; }
  { seen[$1] = 1; }
  END { for (chr in expected) if (!(chr in seen)) print chr; }
```

> 12. How would you convert a VCF file to the Plink binary format? How would you do the reverse, and what kind of problems do you anticipate?

```sh
plink2 --vcf test.vcf.gz --make-bed --out test
```

```sh
plink --bfile test \
      --recode vcf \
      --out testVCF
```
Plink should be given a reference allele file to specify reference allele.


> 13. Write a snippet to reformat a PED file so as to output a file with the following header sample_name genotype_SNP1 genotype_SNP2 ... where genotypes are coded VCF-style (e.g A/C, the order of the alleles in the output is not important).

```sh
awk 'BEGIN {OFS="\t"} NR==1 {printf "sample_name"; for (i=1; i<=(NF-6)/2; i++) printf "\tgenotype_SNP" i; print ""} {printf "%s", $2; for (i=7; i<NF; i+=2) printf "\t%s/%s", $i, $(i+1); print ""}' input.ped > output.tsv

```

> 15. This file contains eQTL overlap data for SNPs that arise as signals in GWAS for several phenotypes. Reformat this file to have one line per SNP/phenotype pair, and two additional columns formatted as such : GENE1(tissue1, tissue2),GENE2(tissue1, tissue3), and GENE1(2),GENE2(2). Each line should contain the SNP/phenotype pair, all genes found overlapping and their respective tissues, and all genes found overlapping with the number of tissues.

[Reformated_output_file](output/Bioinformatics_Q15.txt)

---

### Statistical genetics

> 3. A common practice when performing genetic association studies is to perform an ethnicity check as a quality control step. Can you explain how this is done?
> - You are dealing with whole-genome sequencing data in 2,326 Bulgarian samples. How would you perform such an ethnicity check, and which projection dataset would you use?
  
- 1000G European (EUR)
- GnomAD Non-Finnish European (NFE)
- HGDP Central/Eastern European populations

Running PCA with `smartpca` (EIGENSOFT)
```sh
smartpca -i bulgarian_data.ped -a reference_data.ped -o pca_output
```
Principal components can be generated for visualization.

> 4. You perform a single-point association with blood lipids and find a variant with MAF=0.7% associated at p=1e-14. Do you expect the effect size to be large or small? What would be your next steps investigating this signal?

Expectation is that a **larger effect size** will generate a **smaller p-value**.
Next steps could be:
- Effect Size estimmation
- Conditional analysis as strong associations can be due to LD with a known signal.
- Replication in other datasets
- Meta-analyses

> 6. An analyst studies a population of remote villages in Eastern Europe. They are interested in a particular variant, and compare the frequency in their villages (3.5%) to the EUR population frequency in the 1000 Genomes (0.03%). They conclude that the variant has increased in frequency in their villages. Do you agree, and if not, what would your advice be?

We need to assess whether the observed difference is statistically significant, taking population size into account. Additionally, remote villages may have distinct ancestry compared to the EUR reference population, which could influence the results.

---

