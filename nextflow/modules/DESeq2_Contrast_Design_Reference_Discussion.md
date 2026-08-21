# DESeq2 Contrast Design Reference — HPrEC1/HPrEC2 GrowthHormone vs Control

## Dataset Structure

- 2 cell lines (HPrEC1, HPrEC2) × 2 treatments (Control, GrowthHormone) × 3 replicates each
- Batch == Replicate (A=Rep1, B=Rep2, C=Rep3)
- PCA: PC1 (56.3%) = cell line effect, PC2 (23.2%) = batch effect
- Treatment effect not visible in PCA — buried under cell line + batch noise
- Batch correction is therefore essential before testing

---

## Design Definition

### Standard DESeq2 Design (intercept-based)

```
~ Batch + Comparisons
```

- Intercept = reference level (e.g. BatchA, HPrEC1.Control)
- All other coefficients are relative to the intercept
- Simple pairwise contrasts use: `c("Comparisons", "numerator", "denominator")`
- Cannot easily express multi-group or additive contrasts with `c()` syntax

### Parser Design (this pipeline)

```
~ Batch + Comparisons  (same formula)
```

- Contrasts are expressed as human-readable arithmetic strings e.g. `"HPrEC1.GrowthHormone - HPrEC1.Control"`
- Parser builds a numeric contrast vector over the full model matrix by averaging model matrix rows per group (`colMeans`) then evaluating the arithmetic expression — see Section 1 of `extract_deseq2_results()`
- `batch_vars` triggers zeroing of batch coefficients from the contrast vector so batch does not contribute to the biological comparison
- More general than `c()` — handles 3+ term, additive, and interaction contrasts that `c()` cannot express at all
- Standard `c()` contrast is still run in parallel as cross-validation for simple pairwise cases — see Section 3

### `"0+Comparisons"` vs `"~Batch+Comparisons"`

|-----------------------|----------------|--------------|-------------------------|
| Design                | Batch in model | Parser works | Recommended             |
|-----------------------|----------------|--------------|-------------------------|
| `0+Comparisons`       | ❌             | ✅           | ❌ batch ignored       |
| `0+Batch+Comparisons` | ✅             | ❌           | ❌ breaks group naming |
| `~Batch+Comparisons`  | ✅             | ✅           | ✅ correct choice      |
|-----------------------|----------------|--------------|-------------------------|

- `0+Comparisons` fits group means directly but excludes Batch — batch variation inflates noise and reduces power
- `0+Batch+Comparisons` would break parser because `design_vars` would unite Batch+Comparisons into combined group names like `"A_HPrEC1.Control"` which do not match the contrast strings
- `~Batch+Comparisons` (with intercept) is correct — Batch is in the model for dispersion estimation, contrast strings still reference Comparisons levels directly, `batch_vars` zeroes batch coefficients

---

## Contrast Dictionary

### (1) GH Effect Within Each Cell Line — primary question

> "What genes does Growth Hormone turn on or off?"
> Answered per cell line independently, batch corrected.

**Standard DESeq2:**
```r
c("Comparisons", "HPrEC1.GrowthHormone", "HPrEC1.Control")
c("Comparisons", "HPrEC2.GrowthHormone", "HPrEC2.Control")
```

**Parser (this pipeline):**
```yaml
contrasts:
  - "HPrEC1.GrowthHormone-HPrEC1.Control"
  - "HPrEC2.GrowthHormone-HPrEC2.Control"
```

> Note: For simple pairwise contrasts like these, standard and parser methods produce identical results — cross-validation in Section 3 confirms this. Parser is used for consistency across all contrasts.

---

### (2) Consistent GH Effect Across Both Cell Lines — preferred approach

> "Which GH effects are reproducible regardless of cell line?"

**Recommended — overlap of DEG lists (NOT a contrast):**
```r
consistent_genes <- intersect(
  degs_hprec1 %>% filter(padj < 0.1) %>% pull(gene_id),
  degs_hprec2 %>% filter(padj < 0.1) %>% pull(gene_id)
)
# Further filter to same-direction genes for high-confidence GH response
```

**Standard DESeq2:** cannot express this — no single `c()` contrast exists

**Parser (averaged contrast — NOT recommended, shown for reference only):**
```yaml
- "(HPrEC1.GrowthHormone+HPrEC2.GrowthHormone)/2-(HPrEC1.Control+HPrEC2.Control)/2"
```

**Why NOT use the averaged contrast?**
- A gene strongly up in HPrEC1 but strongly down in HPrEC2 cancels out to non-significant — even though GH clearly acts on it in both lines
- The overlap approach requires independent significance + same direction, which is a stricter and more biologically meaningful definition of "consistent"
- Averaged contrast can also let noise from both lines combine into false positive signal

---

### (3) Interaction — Does GH Effect Differ Between Cell Lines?

> "Which genes respond to GH more strongly in HPrEC2 vs HPrEC1?"
> This is a difference-of-differences, NOT a GH up/down list.

```
(HPrEC2.GrowthHormone - HPrEC2.Control) - (HPrEC1.GrowthHormone - HPrEC1.Control)
```

**Standard DESeq2:** cannot express with `c()` — requires numeric vector:
```r
# get coefficient names first
resultsNames(dds)
# then manually construct numeric vector e.g. c(1, -1, -1, 1) over relevant coefficients
```

**Parser (this pipeline):**
```yaml
- "(HPrEC2.GrowthHormone-HPrEC2.Control)-(HPrEC1.GrowthHormone-HPrEC1.Control)"
```

|-----------------|--------------------------------------------|
| Result          | Interpretation                             |
|-----------------|--------------------------------------------|
| Positive LFC    | Stronger GH response in HPrEC2 than HPrEC1 |
| Negative LFC    | Stronger GH response in HPrEC1 than HPrEC2 |
| Not significant | GH response is similar in both lines       |
|-----------------|--------------------------------------------|

---

### (4) Cell Line Differences — for reference only

> "Which genes are different between HPrEC1 and HPrEC2 regardless of GH?"
> These answer cell line biology, not GH biology. Only run if characterising the two lines is a secondary aim.

**Standard DESeq2:**
```r
c("Comparisons", "HPrEC2.GrowthHormone", "HPrEC1.GrowthHormone")
c("Comparisons", "HPrEC2.Control",        "HPrEC1.Control")
```

**Parser (this pipeline):**
```yaml
- "HPrEC2.GrowthHormone-HPrEC1.GrowthHormone"   # under GH treatment
- "HPrEC2.Control-HPrEC1.Control"               # under control
```

---

## Final Config (Recommended)

```yaml
deseq2:
  design: "~Batch+Comparisons"
  batch_vars: ["Batch"]
  contrasts:
    - "HPrEC1.GrowthHormone-HPrEC1.Control"
    - "HPrEC2.GrowthHormone-HPrEC2.Control"
```
