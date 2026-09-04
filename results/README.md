# Result index

Published results: **0**.

`index.jsonl` contains no entries. Do not add placeholder research claims to the index.

For a future specifically authorized release, use one entry per independently verified result:

```json
{
  "result_id": "<public-result-id>",
  "title": "<precise-result-title>",
  "source_entry_id": "<original-source-entry-id>",
  "statement_scope": "<exact-proved-or-disproved-statement>",
  "source_references": ["<verified-public-source-link>"],
  "artifact_path": "<public-release-directory>",
  "verification_path": "<public-verification-summary>",
  "release_version": "<public-git-commit-or-tag>",
  "released_on": "<YYYY-MM-DD>"
}
```

This template is not a result or release authorization. Public verification summaries must be prepared for the authorized package and must not expose private repository locations, unpublished candidate lists, or private research history.
