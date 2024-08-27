# Missing variants explanation

### Genotype Data After Sample Removal

The original data (before sample removal):

``` text
FAM001 1 0 0 1 2 AA AT GG CC AT AA AA AT GG AA
FAM001 2 0 0 1 2 GG TT GC CC AA GG GG TT TT GC
FAM001 3 0 0 1 2 AG TT GG 00 AA 00 AG 00 GG GG
FAM001 4 0 0 1 2 AA TT GG CC AA AA AA TT GG GG
FAM001 5 0 0 1 2 GG TT GG CC AT GG GG TT GG TT
FAM001 6 0 0 1 2 AG TT GG CC AA AG AG TT GG GG
FAM001 7 0 0 1 2 GG TT GC CC AA GG GG TT GC TT
FAM001 8 0 0 1 2 AA TT GG CC AT AA AA TT GG GG
FAM001 9 0 0 1 2 GG TT GG CC AA GG GG TT GG TT
FAM001 10 0 0 1 2 AA TT GG CT AT AA AA TT GG GG
```

**After removing IID 1, IID 2, and IID 7,** the remaining genotypes are:

``` text
FAM001 3 0 0 1 2 AG TT GG 00 AA 00 AG 00 GG GG
FAM001 4 0 0 1 2 AA TT GG CC AA AA AA TT GG GG
FAM001 5 0 0 1 2 GG TT GG CC AT GG GG TT GG TT
FAM001 6 0 0 1 2 AG TT GG CC AA AG AG TT GG GG
FAM001 8 0 0 1 2 AA TT GG CC AT AA AA TT GG GG
FAM001 9 0 0 1 2 GG TT GG CC AA GG GG TT GG TT
FAM001 10 0 0 1 2 AA TT GG CT AT AA AA TT GG GG
```

### Allele Frequency Calculation After Sample Removal

Let's calculate the MAF for **Variant 3** and others after removing these samples:

1. **Variant 1:** 
   - Alleles: {AG, AA, GG, GG, AA, GG, AA}
   - G: 8, A: 6 
   - MAF: 6/14 = 0.43 **(Retained)**

2. **Variant 2:** 
   - Alleles: {TT, TT, TT, TT, TT, TT, TT}
   - T: 14, A: 0 
   - MAF: 0/14 = 0.00 **(Removed)**

3. **Variant 3:** 
   - Alleles: {GG, GG, GG, GG, GG, GG, GG}
   - G: 14, C: 0 
   - MAF: 0/14 = 0.00 **(Removed)**

4. **Variant 4:** 
   - Alleles: {00, CC, CC, CC, CC, CC, CT}
   - C: 12, T: 2
   - MAF: 2/14 = 0.14 **(Retained)**

5. **Variant 5:** 
   - Alleles: {AA, AA, AT, AA, AT, AA, AT}
   - A: 11, T: 3
   - MAF: 3/14 = 0.21 **(Retained)**

6. **Variant 6:** 
   - Alleles: {AA, AA, AG, AG, AA, GG, AA}
   - A: 12, G: 2
   - MAF: 2/14 = 0.14 **(Retained)**

7. **Variant 7:** 
   - Alleles: {AG, AA, GG, GG, AA, GG, GG}
   - G: 10, A: 4
   - MAF: 4/14 = 0.29 **(Retained)**

8. **Variant 8:** 
   - Alleles: {00, TT, TT, TT, TT, TT, TT}
   - T: 12, A: 0
   - MAF: 0/14 = 0.00 **(Removed)**

9. **Variant 9:** 
   - Alleles: {GG, GG, GG, GG, GG, GG, GG}
   - G: 14, C: 0
   - MAF: 0/14 = 0.00 **(Removed)**

10. **Variant 10:** 
   - Alleles: {GG, GG, GG, GG, GG, GG, GG}
   - G: 14, T: 0
   - MAF: 0/14 = 0.00 **(Removed)**

### Conclusion:
**Variant 3** was removed because the MAF dropped to 0.00 after the samples with
minor alleles (C in this case) were removed. This is a direct consequence of
filtering by `--mind 0.05`, which led to recalculating the MAFs, ultimately
triggering the `--maf 0.05` filter to remove the variant.