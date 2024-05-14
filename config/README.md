# Workflow configuration

Configure this workflow by making two settings in `config.yml`:

- `pmc_baseline_date`: Set the date string for the baseline files, e.g., "2023-12-18".
To see the latest available baseline files, see https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/xml/ or
https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_noncomm/xml/.
- `subset`: Set the list of strings for PMC ID prefix subsets to be considered when running the workflow, e.g.
    - "001"
    - "002"
    - "003"
    - "004"
    - "005"
    - "006"
    - "007"
    - "008"
    - "009"
    - "010"
