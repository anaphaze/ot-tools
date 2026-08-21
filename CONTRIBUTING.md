# Contributing

Thanks for contributing to this repository, maintained by the **Gan Lab (Anaphase)** —
a structural cell biology lab studying the inner workings of the cell nucleus.
This guide covers how we work, with particular attention to how AI-assisted code
is documented so that our software remains transparent and citable in publications.

- Lab: https://github.com/anaphaze
- Website: https://www.anaphase.org

## General guidelines

- Open an issue before starting substantial work, so we can discuss the approach.
- Keep pull requests focused; describe what changed and why.
- Include or update tests for any code that produces analysis results.
- Match the existing code style of the file you are editing.

## AI-assisted code

Some code in this project is generated with the assistance of a large language
model (Claude, by Anthropic). This is permitted, but any AI-generated or
AI-substantially-rewritten source file **must carry a provenance header** at the
top of the file, and the contributor remains fully responsible for reviewing,
testing, and validating that code before it is committed.

The purpose of this policy is transparency and reproducibility. LLM output cannot
be reproduced from a prompt, so the code we commit — not the prompt — is the record.
The header makes it clear which code was AI-assisted, with which model, and when, so
that this can be disclosed accurately in any resulting paper.

### Required header

Prepend a header using the correct comment syntax for the file's language
(`#` for Python and shell, `//` or `/* */` for C/C++/JavaScript, etc.). For
Python:

```python
# -----------------------------------------------------------------------------
# Generated with the assistance of Claude (Anthropic).
# Model: <model-id and version>        Date: <YYYY-MM-DD>
# Reviewed, tested, and validated by the Gan Lab (Anaphase); the authors take
# full responsibility for the correctness of this code.
# Lab: https://github.com/anaphaze | https://www.anaphase.org
# -----------------------------------------------------------------------------
```

Notes:

- Fill in the exact model identifier/version you used and the date of generation.
  A model cannot always report its own version reliably, so set these values
  yourself rather than trusting an auto-filled guess.
- If only part of a file is AI-generated, name the affected functions or sections
  in the header rather than implying the whole file was generated.
- Do not list an AI system as an author of the code or of any resulting paper.

## Citing this software in papers

When AI-assisted code contributes to a publication:

- Include a disclosure statement (Methods or Acknowledgements) naming the tool,
  the model and version, and the date of use.
- Release the actual code, archived with a persistent identifier (e.g., a Zenodo DOI),
  rather than relying on prompts for reproducibility.
- Describe how the code was validated (tests, reference comparisons, manual review).
