# Cell-count analysis summary

## Sample and population summary

- Samples: 10,500
- Populations: 5
- Each nonzero sample's population percentages sum to 100%.

## Responder comparison

Cohort: melanoma + miraclib + PBMC. Percentages were averaged within each subject before a two-sided Mann-Whitney U test; Benjamini-Hochberg correction was applied across populations.

- Independent subjects analyzed: 656
- Populations significant at FDR < 0.05: **None**

| population | n_responder_subjects | n_nonresponder_subjects | p_value | p_value_bh | significant_fdr_0_05 |
| --- | --- | --- | --- | --- | --- |
| b_cell | 331 | 325 | 0.345796 | 0.432245 | False |
| cd4_t_cell | 331 | 325 | 0.0124221 | 0.0621105 | False |
| cd8_t_cell | 331 | 325 | 0.622144 | 0.622144 | False |
| monocyte | 331 | 325 | 0.264492 | 0.432245 | False |
| nk_cell | 331 | 325 | 0.126745 | 0.316862 | False |

## Baseline cohort query

Filters: melanoma + PBMC + time_from_treatment_start = 0 + miraclib.

- Matching samples: 656

### Samples by project

| project | sample_count |
| --- | --- |
| prj1 | 384 |
| prj3 | 272 |

### Distinct subjects by response

| response | unique_subject_count |
| --- | --- |
| no | 325 |
| yes | 331 |

### Distinct subjects by sex

| sex | unique_subject_count |
| --- | --- |
| F | 312 |
| M | 344 |
