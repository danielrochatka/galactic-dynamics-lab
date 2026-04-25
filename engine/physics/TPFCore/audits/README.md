# TPFCore audits

This folder stores **audit evidence records** for contamination/sanity checks in the TPFCore spike work.

## Scope and intent

- These records are **evidence logs**, not theory claims.
- These records are **not** validation/proof statements for TPF.
- These records are intended for spike-branch transparency and traceability.

## Evidence standard

Audit notes must be grounded in executable artifacts, not narrative statements:

- cite executable code paths and function call chains,
- cite concrete files/functions/line ranges inspected,
- list terminal commands used,
- include observed outputs (when applicable).

## Non-evidence sources

The following are useful context but are **not proof of behavior** by themselves:

- comments,
- README prose,
- run_info text,
- PR summaries.

Behavior claims must still be tied back to executable code paths and observed command output.

## Required finding buckets

Each audit should separate findings into:

- **A. confirmed clean**
- **B. suspicious but harmless**
- **C. actual contamination / must fix**
- **D. unsupported / not found**

## Current records

- `static_4d_residual_contamination_audit.md` — contamination/sanity audit of `tpf_4d_static_residual_benchmark` dispatch and static residual path.

