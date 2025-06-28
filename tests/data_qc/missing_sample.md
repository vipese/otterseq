# Missing samples

### Step 1: Understanding the Warnings

PLINK provided two key warnings:

1. **Warning: Variant 7 (post-sort) triallelic; setting rarest alleles missing.**
2. **Warning: Variant 10 (post-sort) quadallelic; setting rarest alleles missing.**

This indicates that for these two variants, PLINK found more than two alleles across
the samples, which is not typical for biallelic SNPs. As a result, PLINK set the rarest alleles to missing.

### Step 2: Identifying the Affected Variants

Let’s determine which variants are triallelic or quadallelic:

1. **Variant 7**: 
   - After sorting, PLINK determined this variant is triallelic.
   - Affected by multiple alleles at this position.

2. **Variant 10**:
   - PLINK determined this variant is quadallelic.
   - Affected by four different alleles at this position.

### Step 3: Calculate Missing Data Due to Triallelic/Quadallelic Variants

To understand why specific individuals were removed, let's analyze the genotypes at Variants 7 and 10.

Given the `.ped` file:
```
FAM001 1 0 0 1 2 AA AT GG CC AT AA AA AT GG AA
FAM001 2 0 0 1 2 GG TT GC CC AA GG GG TT TT GC
FAM001 3 0 0 1 2 AG TT GG CT AA AG AG AT GG GG
FAM001 4 0 0 1 2 AA TT GG CC AA AA AA TT GG GG
FAM001 5 0 0 1 2 GG TT GG CC AT GG GG TT GG TT
FAM001 6 0 0 1 2 AG TT GG CC AA AG AG TT GG GG
FAM001 7 0 0 1 2 GG TT GC CC AA GG GG TT GC TT
FAM001 8 0 0 1 2 AA TT GG CC AT AA AA TT GG GG
FAM001 9 0 0 1 2 GG TT GG CC AA GG GG TT GG TT
FAM001 10 0 0 1 2 AA TT GG CT AT AA AA TT GG GG
```

#### Variant 7 Genotypes:

- **FAM001 1**: AA
- **FAM001 2**: GC
- **FAM001 3**: AG
- **FAM001 4**: AA
- **FAM001 5**: GG
- **FAM001 6**: AG
- **FAM001 7**: GC
- **FAM001 8**: AA
- **FAM001 9**: GG
- **FAM001 10**: AA

This variant has three alleles: A, G, C. PLINK flags this as triallelic, so the
rarest alleles (C) are marked as missing.

#### Variant 10 Genotypes:

- **FAM001 1**: AA
- **FAM001 2**: GC
- **FAM001 3**: GG
- **FAM001 4**: GG
- **FAM001 5**: TT
- **FAM001 6**: GG
- **FAM001 7**: TT
- **FAM001 8**: GG
- **FAM001 9**: TT
- **FAM001 10**: GG

This variant has four alleles: A, G, C, T. PLINK flags this as quadallelic, so
the rarest alleles (in this case, potentially A or C) are marked as missing.

### Step 4: Calculate Missing Data for Each Individual

Now, let’s count the number of genotypes marked as missing due to these
triallelic and quadallelic variants:

1. **FAM001 1**:
   - Missing in Variant 7: No (AA, common allele)
   - Missing in Variant 10: Yes (AA, rare allele)
   - **Total missing genotypes**: 1/10 = 10% missing (Exceeds 5%)

2. **FAM001 2**:
   - Missing in Variant 7: Yes (GC, rare allele)
   - Missing in Variant 10: Yes (GC, rare allele)
   - **Total missing genotypes**: 2/10 = 20% missing (Exceeds 5%)

3. **FAM001 3**:
   - Missing in Variant 7: No (AG, common alleles)
   - Missing in Variant 10: No (GG, common alleles)
   - **Total missing genotypes**: 0/10 = 0% missing (Below 5%)

4. **FAM001 4**:
   - Missing in Variant 7: No (AA, common allele)
   - Missing in Variant 10: No (GG, common alleles)
   - **Total missing genotypes**: 0/10 = 0% missing (Below 5%)

5. **FAM001 5**:
   - Missing in Variant 7: No (GG, common allele)
   - Missing in Variant 10: No (TT, common alleles)
   - **Total missing genotypes**: 0/10 = 0% missing (Below 5%)

6. **FAM001 6**:
   - Missing in Variant 7: No (AG, common alleles)
   - Missing in Variant 10: No (GG, common alleles)
   - **Total missing genotypes**: 0/10 = 0% missing (Below 5%)

7. **FAM001 7**:
   - Missing in Variant 7: Yes (GC, rare allele)
   - Missing in Variant 10: Yes (TT, rare allele)
   - **Total missing genotypes**: 2/10 = 20% missing (Exceeds 5%)

8. **FAM001 8**:
   - Missing in Variant 7: No (AA, common allele)
   - Missing in Variant 10: No (GG, common alleles)
   - **Total missing genotypes**: 0/10 = 0% missing (Below 5%)

9. **FAM001 9**:
   - Missing in Variant 7: No (GG, common allele)
   - Missing in Variant 10: No (TT, common alleles)
   - **Total missing genotypes**: 0/10 = 0% missing (Below 5%)

10. **FAM001 10**:
    - Missing in Variant 7: No (AA, common allele)
    - Missing in Variant 10: No (GG, common alleles)
    - **Total missing genotypes**: 0/10 = 0% missing (Below 5%)

### Step 5: Conclusion

- **FAM001 1, FAM001 2, and FAM001 7** were removed because their missing
genotype rates exceeded the `--mind 0.05` threshold. Specifically:
  - **FAM001 1**: 10% missing genotypes.
  - **FAM001 2**: 20% missing genotypes.
  - **FAM001 7**: 20% missing genotypes.

This is a direct result of the rare alleles in triallelic and quadallelic
variants being set as missing, which increased the missing data for these individuals
beyond the allowed threshold, causing their removal by PLINK.